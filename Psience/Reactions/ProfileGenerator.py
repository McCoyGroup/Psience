import abc

from ..Molecools import Molecule
from ..Molecools.Evaluator import EnergyEvaluator

class ProfileGenerator:
    def __init__(self,
                 reactants:list[Molecule],
                 energy_evaluator:EnergyEvaluator
                 ):
        self.reactants = reactants
        self.eval = energy_evaluator

    @abc.abstractmethod
    def generate(self, **opts):
        ...

class OptimizingProfileGenerator(ProfileGenerator):

    def generate(self, initial_search_dir, **opts):
        ...


class NudgedElasticBand(ProfileGenerator):
    ...