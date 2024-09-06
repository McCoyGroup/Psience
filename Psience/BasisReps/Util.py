
import numpy as np

__all__ = [
    "StateMaker"
]

class StateMaker:
    """
    A tiny but useful class to make states based on their quanta
    of excitation
    """

    def __init__(self, ndim, mode='low-high'):
        self.ndim = ndim
        self.mode = mode

    @classmethod
    def parse_state(cls, state, mode='low-high'):
        nzp = np.nonzero(state)
        if len(nzp) > 0: nzp = nzp[0]
        if len(nzp) == 0: return "()"
        if mode == 'low-high':
            pos = len(state) - nzp
            state_quants = reversed(list(zip(pos, nzp)))
        elif mode == 'high-low':
            pos = nzp + 1
            state_quants = list(zip(pos, nzp))
        else:
            pos = nzp
            state_quants = list(zip(pos, nzp))
        return "".join("{}({})".format(p, state[i]) for p,i in state_quants)

    def make_state(self, *specs, mode=None):

        if mode is None:
            mode = self.mode

        state = [0] * self.ndim
        for s in specs:
            if isinstance(s, (int, np.integer)):
                i = s
                q = 1
            else:
                i,q = s
            if mode == 'low-high':
                state[-i] = q
            elif mode == 'high-low':
                state[i-1] = q
            elif mode == 'normal':
                state[i] = q
            else:
                raise ValueError("don't know what to do with filling mode '{}'".format(mode))
        return state

    def __call__(self, *specs, mode=None):
        return self.make_state(*specs, mode=mode)