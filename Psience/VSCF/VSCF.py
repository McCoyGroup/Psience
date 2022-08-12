import abc, dataclasses, typing, numpy as np
from McUtils.Zachary import Mesh

__all__ = [
    "SelfConsistentOptimizer",
    "GridSCF"
]

class SCFComponentWavefunctions(typing.Protocol):
    def expectation(self, operator)->'typing.Callable|np.ndarray':
        ...
    def overlap(self, wfn:'SCFComponentWavefunctions')->float:
        ...
    def __getitem__(self, item)->'SCFComponentWavefunctions':
        ...
class SCFSolver(typing.Protocol):
    def __call__(self, pot, **kwargs)->SCFComponentWavefunctions:
        ...

@dataclasses.dataclass
class SelfConsistentIterationData:
    wavefunctions:'list[SCFComponentWavefunctions]'
    iteration:int

class SelfConsistentOptimizer:
    __props__ = (
        'initial_point',
        'overlap_threshold',
        'max_iterations'
    )
    def __init__(self,
                 potential,
                 solvers_1D:'Iterable[SCFSolver]',
                 initial_point=None,
                 overlap_threshold=.999999,
                 max_iterations=100
                 ):
        self.pot = potential
        self.solv = solvers_1D
        if initial_point is None:
            initial_point = self.get_initial_point()
        self.center = initial_point
        self.threshold = overlap_threshold
        self.its = max_iterations
    @abc.abstractmethod
    def get_initial_point(self):
        raise NotImplementedError("...")
    @abc.abstractmethod
    def get_SCF_potential(self, axis:int, wavefunctions:'Iterable[SCFComponentWavefunctions]', target:'Iterable[int]'):
        raise NotImplementedError("...")
    @abc.abstractmethod
    def get_cut_potential(self, axis, center):
        raise NotImplementedError("...")
    def initialize(self):
        wfns = [
            s(self.get_cut_potential(i, self.center))
            for i,s in enumerate(self.solv)
        ]
        return SelfConsistentIterationData(wfns, 0)
    def step(self, iteration_data:SelfConsistentIterationData, target:'Iterable[int]'):
        new_wfns = [
            solver(
                self.get_SCF_potential(i, iteration_data.wavefunctions, target)
            )
            for i, solver in enumerate(self.solv)
        ]
        return SelfConsistentIterationData(new_wfns, iteration_data.iteration+1)
    def check_overlaps(self,
                       previous:SelfConsistentIterationData,
                       current:SelfConsistentIterationData,
                       threshold,
                       target
                       ):
        # ov = [
        #     old[(t,)].overlap(new[(t,)])[0][0]
        #     for old,new,t in zip(previous.wavefunctions, current.wavefunctions, target)
        # ]
        # print(ov)
        # for old,new,t in zip(previous.wavefunctions, current.wavefunctions, target):
        #     old[:1].plot(
        #         figure=new[:1].plot()
        #     ).show()

        return all(
            old[(t,)].overlap(new[(t,)])[0][0] > threshold
            for old,new,t in zip(previous.wavefunctions, current.wavefunctions, target)
        )
    def run(self, target=None):
        if target is None:
            target = [0] * len(self.solv)
        cur = self.initialize()
        # print(cur.wavefunctions[0].frequencies() * 219475)
        # fig = None
        # for old in cur.wavefunctions:
        #     fig = old[:1].plot(figure=fig)
        # fig.show()
        for _ in range(self.its):
            old = cur
            cur = self.step(old, target)
            if self.check_overlaps(old, cur, self.threshold, target):
                break
        else:
            raise ValueError("SCF didn't converge in {} steps".format(self.its))
        return cur

class GridSCF(SelfConsistentOptimizer):
    """
    A more concrete implementation of an SCF solver
    which uses a fixed grid of points and potential
    values at those points to construct the averaged potential
    """
    def __init__(self, mesh, potential_values, solvers, **opts):
        self.grid = Mesh(mesh)
        self.vals = potential_values
        super().__init__(None, solvers, **opts)
    def get_initial_point(self, vals=None):
        if vals is None:
            vals = self.vals
        pos = np.unravel_index(np.argmin(vals), vals.shape)
        return pos
    def get_cut_potential(self, axis, center):
        """
        :param axis:
        :type axis:
        :param center: index of the minimum point in the grid
        :type center:
        :return:
        :rtype:
        """
        center = list(center)
        center[axis] = slice(None, None, None)
        center = tuple(center) # numpy needs this
        return self.vals[center]
    def get_SCF_potential(self, axis, wavefunctions:'list[SCFComponentWavefunctions]', target:'Iterable[int]'):
        v = np.moveaxis(self.vals, axis, -1) # shift final remaining axis to the back
        # import McUtils.Plots as plt
        # plt.ContourPlot(*self.grid.transpose(2, 0, 1), v, plot_range=[[1, 2], [1, 2]]).show()
        for i,w in enumerate(wavefunctions):
            if i != axis:
                v = w[(target[i],)].expectation(v).squeeze()
        # plt.Plot(wavefunctions[axis].grid, v).show()
        return v
