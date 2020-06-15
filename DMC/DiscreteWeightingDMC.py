from .AbstractDMC import AbstractDMC

__all__ = ["DiscreteWeightingDMC"]

class DiscreteWeightingDMC(AbstractDMC):

    def branch(self): # this needs to be implemented properly to work here
        raise NotImplemented

    def update_weights(self, potentials, v_refs):
        raise NotImplemented