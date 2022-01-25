
import numpy as np
from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import NullLogger, Checkpointer
from .Common import PerturbationTheoryException, _safe_dot

__all__ = [
    "PerturbationTheoryCorrections"
]



class PerturbationTheoryCorrections:
    """
    Represents a set of corrections from perturbation theory.
    Can be used to correct other operators in the basis of the original calculation.

    """
    def __init__(self,
                 states,
                 coupled_states,
                 total_basis,
                 energy_corrs,
                 wfn_corrections,
                 all_energy_corrections=None,
                 degenerate_states=None,
                 degenerate_transformation=None,
                 degenerate_energies=None,
                 logger=None
                 ):
        """
        :param states:
        :type states: BasisStateSpace
        :param coupled_states:
        :type coupled_states: BasisMultiStateSpace
        :param total_basis:
        :type total_basis: BasisMultiStateSpace
        :param energy_corrs:
        :type energy_corrs: np.ndarray
        :param wfn_corrections:
        :type wfn_corrections: Iterable[SparseArray]
        :param degenerate_states:
        :type degenerate_states: None | np.ndarray
        :param degenerate_transformation:
        :type degenerate_transformation: None | np.ndarray
        :param degenerate_energies:
        :type degenerate_energies: None | np.ndarray
        """
        self.states = states
        self.coupled_states = coupled_states
        self.total_basis = total_basis
        self.energy_corrs = energy_corrs
        self.all_energy_corrs = all_energy_corrections
        self.wfn_corrections = wfn_corrections
        self.degenerate_states = degenerate_states
        self.degenerate_transf = degenerate_transformation
        self.degenerate_energies = degenerate_energies
        self.logger = logger

    @classmethod
    def from_dicts(cls,
                   states,
                   corrections,
                   hamiltonians,
                   **opts
                   ):
        """
        :param states: a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
        :type states: dict
        :param corrections: the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
        :type corrections: dict
        """
        state_space = states['states']
        coupled_states = states['coupled_states']
        total_basis = states['total_states']
        energy_corrs = corrections['energies']
        all_energy_corrs = corrections['energy_corrections'] if 'energy_corrections' in corrections else None
        wfn_corrections = corrections['wavefunctions']
        if 'degenerate_states' in states:
            degenerate_states = states['degenerate_states']
        else:
            degenerate_states = None

        if 'degenerate_transformation' in corrections:
            degenerate_transf = corrections['degenerate_transformation']
        else:
            degenerate_transf = None

        if 'degenerate_energies' in corrections:
            degenerate_energies = corrections['degenerate_energies']
        else:
            degenerate_energies = None

        return cls(
            state_space,
            coupled_states,
            total_basis,
            energy_corrs,
            wfn_corrections,
            all_energy_corrections=all_energy_corrs,
            degenerate_states=degenerate_states,
            degenerate_transformation=degenerate_transf,
            degenerate_energies=degenerate_energies,
            **opts
        )

    @property
    def degenerate(self):
        return self.degenerate_transf is not None

    @property
    def energies(self):
        if self.degenerate:
            return self.degenerate_energies
        else:
            return np.sum(self.energy_corrs, axis=1)

    @property
    def order(self):
        return len(self.energy_corrs[0])

    def take_subspace(self, space):
        """
        Takes only those elements that are in space
        :param space:
        :type space:
        :return:
        :rtype:
        """

        new_states = self.states.find(space)
        # raise Exception(new_states)
        return type(self)(
            self.states.take_states(space),
            self.coupled_states.take_states(space),
            self.total_basis,
            self.energy_corrs[new_states],
            [w[new_states, :] for w in self.wfn_corrections],
            # not sure what to do with all this...
            degenerate_states=self.degenerate_states,
            degenerate_transformation=self.degenerate_transf,
            degenerate_energies=self.degenerate_energies,
            logger=self.logger
        )

    @staticmethod
    def _fmt_operator_rep(full_ops, operator_symbol, conversion, real_fmt="{:>.8e}", padding_fmt='{:>16}'):
        tag_line = None
        rep_lines = None

        op_dim = None

        for (a, b, c), subrep in full_ops:
            if isinstance(subrep, SparseArray):
                subrep = subrep.asarray()
            elif isinstance(subrep, (int, float, np.integer, np.floating)):
                if subrep == 0:
                    if op_dim is None:
                        raise ValueError("was lazy and haven't filled operator dim yet...")
                    subrep = np.zeros(op_dim)
                else:
                    raise ValueError("don't know what to do with representation '{}'".format(subrep))

            if op_dim is None:
                op_dim = subrep.shape

            if conversion is not None:
                subrep = subrep * conversion

            subrep_lines = [
                " ".join(padding_fmt.format(real_fmt.format(e)) for e in line)
                for line in subrep
            ]
            line_len = len(subrep_lines[0])

            if rep_lines is None:
                rep_lines = subrep_lines
            else:
                rep_lines = [x + " " + y for x,y in zip(rep_lines, subrep_lines)]

            tag_fmt = "{:<" + str(line_len) + "}"
            base_tag=tag_fmt.format("<{a}|{A}({c})|{b}>".format(A=operator_symbol, a=a, c=c, b=b))
            if tag_line is None:
                tag_line = base_tag
            else:
                tag_line += " " + base_tag

        # we want to return a line list so the logger can add any prefixes it needs
        rep_lines.insert(0, tag_line)
        return rep_lines

    def operator_representation(self, operator_expansion, order=None, subspace=None, contract=True,
                                logger_symbol="A",
                                logger_conversion=None
                                ):
        """
        Generates the representation of the operator in the basis of stored states

        :param operator_expansion: the expansion of the operator
        :type operator_expansion: Iterable[float] | Iterable[np.ndarray]
        :param order: the order of correction to go up to
        :type order: Iterable[float] | Iterable[np.ndarray]
        :param subspace: the subspace of terms in which the operator expansion is defined
        :type subspace: None | BasisStateSpace
        :return: the set of representation matrices for this operator
        :rtype: Iterable[np.ndarray]
        """

        mordor = self.order - 1
        if order is None:
            order = mordor
        if order > mordor:
            raise PerturbationTheoryException("{}: can't correct up to order {} when zero-order states were only corrected up to order {}".format(
                type(self).__name__,
                order,
                mordor
            ))
        order = order + 1 # so that we actually do get up to the request order after accounting for the zeros...
        if len(operator_expansion) < order:
            operator_expansion = list(operator_expansion) + [0]*(order - len(operator_expansion))

        # we stopped supporting indexing based on the total set of inds...
        if subspace is None:
            wfn_corrs = self.wfn_corrections[:order]
        else:
            # need to make the subspace good with the subspace in which the corrections are defined...
            subspace_sel = self.total_basis.find(subspace, check=True)
            wfn_corrs = []
            for k in range(order):
                wfn_corrs.append(self.wfn_corrections[k][:, subspace_sel])

        # generalizes the dot product so that we can use 0 as a special value...
        dot = _safe_dot
        logger = self.logger
        logger = None if logger is None or isinstance(logger, NullLogger) else logger

        # does the dirty work of acutally applying the rep...
        reps = [[] for _ in range(order)]
        full_ops = []
        for k in range(order):
            tags = []
            op = []
            # apply each thing up to requested order...
            for a in range(k+1): # if k == 2: a=0, a=1, a=2
                for b in range(k-a+1): # if k==2, a==0: b=0, b=1, b=2; a==1: b=0, b=1
                    c = k - (a + b) # a + b + c == k
                    rop = operator_expansion[c]
                    if isinstance(rop, (int, float, np.integer, np.floating)): # constant reps...
                        if rop != 0: # cheap easy check
                            subrep = rop * dot(wfn_corrs[a], wfn_corrs[b].T)
                            op.append(subrep)
                        else:
                            subrep = 0
                            op.append(0)
                    else:
                        subrep = dot(dot(wfn_corrs[a], rop), wfn_corrs[b].T)
                        op.append(subrep)

                    full_ops.append([
                        (a, b, c),
                        subrep
                    ])

            if contract:
                op = sum(op)
            reps[k] = op

        if logger is not None:
            logger.log_print(full_ops, logger_symbol, logger_conversion, message_prepper=self._fmt_operator_rep)

        return reps

    def get_overlap_matrices(self):
        """
        Returns the overlap matrices for the set of corrections
        at each order of correction

        :return:
        :rtype:
        """

        wat = []
        for k in range(2 + 1):
            ov = None
            for i in range(k + 1):
                c1 = self.wfn_corrections[i].asarray()
                c2 = self.wfn_corrections[k - i].asarray()
                if ov is None:
                    ov = np.dot(c1, c2.T)
                else:
                    ov += np.dot(c1, c2.T)
            wat.append(ov)

        return wat

    # def checkpoint_save(self, checkpoint:Checkpointer):
    #     """
    #     Writes correction arrays to checkpoint
    #
    #     :param checkpoint:
    #     :type checkpoint:
    #     :return:
    #     :rtype:
    #     """
    #     import gc, sys
    #
    #     # self.states = states
    #     # self.coupled_states = coupled_states
    #     # self.total_basis = total_basis
    #     # self.energy_corrs = energy_corrs
    #     # self.all_energy_corrs = all_energy_corrections
    #     # self.wfn_corrections = wfn_corrections
    #     # self.degenerate_states = degenerate_states
    #     # self.degenerate_transf = degenerate_transformation
    #     # self.degenerate_energies = degenerate_energies
    #     # self.logger = logger
    #     checkpoint["wfn_corrections"] = self.wfn_corrections
    #     print(">>>>", sys.getrefcount(self.wfn_corrections))
    #     self.wfn_corrections = None
    #     gc.collect()
    #
    # def checkpoint_reload(self, checkpoint:Checkpointer):
    #     self.wfn_corrections = checkpoint['wfn_corrections']
    #
    # def disk_backed(self):
    #     return self.chk_backer(self)
    #
    # class chk_backer:
    #     def __init__(self, parent):
    #         self.parent = parent
    #         self.chk = None
    #     def load_chk(self):
    #         import tempfile as tf
    #         target = tf.NamedTemporaryFile(suffix=".hdf5").name
    #         self.chk = Checkpointer.from_file(target)
    #     def unload_chk(self):
    #         import os
    #         os.remove(self.chk.checkpoint_file)
    #     def __enter__(self):
    #         self.load_chk()
    #         self.chk.__enter__()
    #         self.parent.checkpoint_save(self.chk)
    #     def __exit__(self, exc_type, exc_val, exc_tb):
    #         self.parent.checkpoint_reload(self.chk)
    #         self.chk.__exit__(exc_type, exc_val, exc_tb)
    #         self.unload_chk()

    def savez(self, file):
        raise NotImplementedError("old and wrong now")
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections
        )
        if self.degenerate_states is not None:
            keys['degenerate_states'] = self.degenerate_states
        if self.degenerate_transf is not None:
            keys['degenerate_transformation'] = self.degenerate_transf
        if self.degenerate_energies is not None:
            keys['degenerate_energies'] = self.degenerate_energies
        np.savez(file, **keys)
    @classmethod
    def loadz(cls, file):
        keys = np.load(file)
        return cls.from_dicts(
            {
                "states":keys['states'],
                 "coupled_states":keys['coupled_states'],
                 "total_states":keys['total_states'],
                 "degenerate_states":keys['degenerate_states'] if 'degenerate_states' in keys else None
             },
            {
                "energies":keys['energies'],
                "wavefunctions":keys['wavefunctions'],
                "degenerate_transformation": keys['degenerate_transformation'] if 'degenerate_transformation' in keys else None,
                "degenerate_energies": keys['degenerate_energies'] if 'degenerate_energies' in keys else None
            },
            keys['hamiltonians']
        )

    def to_state(self, serializer=None):
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections,
            degenerate_states=self.degenerate_states,
            degenerate_transformations=self.degenerate_transf,
            degenerate_energies=self.degenerate_energies
        )
        return keys
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls.from_dicts(
            {
                "states": serializer.deserialize(data['states']),
                "coupled_states": serializer.deserialize(data['coupled_states']),
                "total_states": serializer.deserialize(data['coupled_states']),
                "degenerate_states": serializer.deserialize(data['degenerate_states']),
            },
            {
                "energies": data['energies'],
                "wavefunctions": data['wavefunctions'],
                "degenerate_transformation": data['degenerate_transformations'],
                "degenerate_energies": data['degenerate_energies']
            },
            data['hamiltonians'] # we probably want to ditch this for memory reasons...
        )