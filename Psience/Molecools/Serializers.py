
from McUtils.Scaffolding import Checkpointer
from McUtils.Devutils import Schema

MoleculeSchema = Schema(
    {
        'atoms': 'array',
        'coords': 'array'
    },
    {
        'bonds': 'array',
        'potential_derivatives': 'array',
        'dipole_derivatives': 'array',
        'normal_modes': {
            'matrix': 'array',
            'inverse': 'array',
            'freqs': 'array',
            'origin': 'array'
        }
    }
)

class MoleculePropertyCache:
    def __init__(self, mol, checkpoint=None):
        if checkpoint is None:
            checkpoint = {}
        self.checkpointer = Checkpointer.build_canonical(checkpoint)
        self.mol = mol

    def preserialize(self):
        self.checkpointer.update({
            'atoms':self.mol.atoms,
            'coords':self.mol.coords
        })

    def __getitem__(self, item):
        return self.checkpointer[item]
    def __setitem__(self, key, value):
        self.checkpointer[key] = value
    def __delitem__(self, key, value):
        self.checkpointer[key] = value
    def update(self, vals):
        self.checkpointer.update(vals)

    def cached_eval(self,
                    key, generator,
                    *,
                    condition=None,
                    args=(),
                    kwargs=None):
        return self.checkpointer.cached_eval(
            key,
            generator,
            condition=condition,
            args=args,
            kwargs=kwargs
        )



