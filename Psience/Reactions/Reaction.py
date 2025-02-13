

__all__ = [
    "Reaction"
]

from ..Molecools import Molecule

class Reaction:
    def __init__(self,
                 reactants:list[Molecule],
                 products:list[Molecule]
                 ):
        self.reactants = reactants
        self.products = products


