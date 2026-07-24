
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
        """
        **LLM Docstring**

        Build a checkpoint-backed cache for a molecule's derived properties.

        Wraps `checkpoint` in a canonical `Checkpointer` (an empty dict is used if none is given) and keeps a reference to the owning `mol`.

        :param mol: the molecule this cache stores properties for
        :type mol: AbstractMolecule
        :param checkpoint: a checkpoint source (e.g. a file path or dict) to build the underlying `Checkpointer` from; if `None`, an empty dict is used
        :type checkpoint: dict | str | None
        :return: None
        :rtype: None
        """
        if checkpoint is None:
            checkpoint = {}
        self.checkpointer = Checkpointer.build_canonical(checkpoint)
        self.mol = mol

    def preserialize(self):
        """
        **LLM Docstring**

        Write the molecule's `atoms` and `coords` into the checkpointer so they are present before the rest of the cache is serialized.

        :return: None
        :rtype: None
        """
        self.checkpointer.update({
            'atoms':self.mol.atoms,
            'coords':self.mol.coords
        })

    def __getitem__(self, item):
        """
        **LLM Docstring**

        Look up a cached value by key, delegating to the underlying checkpointer.

        :param item: the key to look up
        :type item: str
        :return: the cached value stored under `item`
        :rtype: object
        """
        return self.checkpointer[item]
    def __setitem__(self, key, value):
        """
        **LLM Docstring**

        Store a value under `key` in the underlying checkpointer.

        :param key: the key to store the value under
        :type key: str
        :param value: the value to cache
        :type value: object
        :return: None
        :rtype: None
        """
        self.checkpointer[key] = value
    def __delitem__(self, key, value):
        """
        **LLM Docstring**

        Intended to delete a cached entry, but as written it actually **sets** `self.checkpointer[key] = value` (identical to `__setitem__`) rather than deleting anything. This looks like a bug: the extra `value` parameter is accepted but no deletion occurs.

        :param key: the key that would be deleted
        :type key: str
        :param value: unused for deletion; currently assigned to `self.checkpointer[key]` instead
        :type value: object
        :return: None
        :rtype: None
        """
        self.checkpointer[key] = value
    def update(self, vals):
        """
        **LLM Docstring**

        Bulk-update the checkpointer with multiple key/value pairs.

        :param vals: mapping of keys to values to merge into the checkpointer
        :type vals: dict
        :return: None
        :rtype: None
        """
        self.checkpointer.update(vals)

    def cached_eval(self,
                    key, generator,
                    *,
                    condition=None,
                    args=(),
                    kwargs=None):
        """
        **LLM Docstring**

        Evaluate and cache a value under `key`, delegating entirely to the checkpointer's own `cached_eval`.

        :param key: the cache key to look up or populate
        :type key: str
        :param generator: callable used to compute the value when it is not already cached
        :type generator: callable
        :param condition: optional predicate controlling whether the cached value should be recomputed
        :type condition: callable or None
        :param args: positional arguments passed to `generator` if it is called
        :type args: tuple
        :param kwargs: keyword arguments passed to `generator` if it is called
        :type kwargs: dict or None
        :return: the cached or newly computed value
        :rtype: object
        """
        return self.checkpointer.cached_eval(
            key,
            generator,
            condition=condition,
            args=args,
            kwargs=kwargs
        )



