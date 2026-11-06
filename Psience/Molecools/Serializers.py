
import warnings
from McUtils.Scaffolding import Checkpointer, ObjectStore
from McUtils.Devutils import Schema, chained_key

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
    def __init__(self, mol, checkpoint=None, identifier=None):
        """
        **LLM Docstring**

        Build a checkpoint-backed cache for a molecule's derived properties.

        Wraps `checkpoint` in a canonical `Checkpointer` (an empty dict is used if none is given) and keeps a reference to the owning `mol`.

        :param mol: the molecule this cache stores properties for
        :type mol: AbstractMolecule
        :param checkpoint: a checkpoint source (e.g. a file path or dict) to build the underlying `Checkpointer` from; if `None`, an empty dict is used. An `ObjectStore` stores `cached_eval` results as store entries instead; plain item access then uses an in-memory dict. `True` is rejected (it used to be a silent no-op that failed on first use; there is no implicit store location).
        :type checkpoint: dict | str | ObjectStore | bool | None
        :param identifier: what makes this molecule's cached values reusable: a string (or any JSON-able value) that you choose and that must change whenever the molecule does. The cache never derives one from the geometry. With an identifier, `cached_eval` entries are scoped by `(identifier, key, options)` and the molecule's derived properties (energy, dipole, normal modes, point group) are cached automatically; an `ObjectStore` cache requires one to store anything.
        :type identifier: str | object | None
        :return: None
        :rtype: None
        """
        self.store = None
        self.identifier = identifier
        self.source = checkpoint  # what `Molecule.modify` hands to the new molecule's cache
        self.changed = set()  # `modify` arguments / setters used since the identifier was assigned (see `valid`)
        self.persistent = checkpoint is not None and checkpoint is not False
        if checkpoint is True:
            raise ValueError(
                "checkpoint_file=True is not supported; pass a file path, a dict, or an explicit ObjectStore "
                "(for example ObjectStore('/path/to/store'))"
            )
        if isinstance(checkpoint, ObjectStore):
            self.store = checkpoint
            checkpoint = {}
        if checkpoint is None or checkpoint is False:
            checkpoint = {}
        self.checkpointer = Checkpointer.build_canonical(checkpoint)
        self.mol = mol

    def note_change(self, *names):
        """Records that something a cached value may depend on (`'coords'`, `'potential_derivatives'`, ...) was changed in place"""
        self.changed.update(names)

    def valid(self, key):
        """
        Whether the identifier still vouches for `key`: false once something that `key` depends on (the
        molecule's `cascade`; for keys it doesn't list, anything at all) has been changed by `Molecule.modify`,
        a setter or a property manager without a new identifier being given. Such values are computed, not cached.
        """
        return not self.mol.cascade.is_stale(key, self.changed)

    @property
    def active(self):
        """Whether derived properties are cached automatically: a persistent checkpoint and an identifier"""
        return self.persistent and self.identifier is not None

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
    def __delitem__(self, key):
        """
        **LLM Docstring**

        Delete a cached entry from the underlying checkpointer.

        :param key: the key to delete
        :type key: str
        :return: None
        :rtype: None
        """
        del self.checkpointer[key]
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
                    kwargs=None,
                    options=None):
        """
        **LLM Docstring**

        Evaluate and cache a value for this molecule. With an `identifier` the entry is scoped by `(identifier, key, options)`, so a different identifier or different options miss the cache; the identifier is yours to keep in step with the molecule (nothing is derived from the geometry). Without an identifier and without options, the checkpointer is used with the plain `key`, as it always was (an `ObjectStore` needs an identifier and raises `ValueError` without one). With an `ObjectStore` the value is a store entry (a value the store can't serialize is returned uncached, with a warning); otherwise it is saved in the checkpointer under `key@<scope>`.

        :param key: the property name
        :type key: str
        :param generator: callable used to compute the value when it is not already cached
        :type generator: callable
        :param condition: optional predicate on `key`; when it is false the cache is bypassed
        :type condition: callable or None
        :param args: positional arguments passed to `generator` if it is called
        :type args: tuple
        :param kwargs: keyword arguments passed to `generator` if it is called
        :type kwargs: dict or None
        :param options: anything besides the molecule that the value depends on
        :type options: object
        :return: the cached or newly computed value
        :rtype: object
        """
        if kwargs is None:
            kwargs = {}
        if self.store is None and self.identifier is None and options is None:
            return self.checkpointer.cached_eval(key, generator, condition=condition, args=args, kwargs=kwargs)
        if (condition is not None and not condition(key)) or not self.valid(key):
            return generator(*args, **kwargs)
        ident = self.identifier
        if self.store is not None and ident is None:
            raise ValueError("an ObjectStore property cache needs an identifier: pass `checkpoint_id=...` to the "
                             "Molecule (nothing is derived from the geometry)")
        if self.store is not None:
            entry = self.store.key_for(ident, key, options)
            replace = False
            if self.store.has(entry):
                try:
                    return self.store.get(entry)
                except Exception as err: # an entry that can't be read back is recomputed and replaced
                    warnings.warn("couldn't load cached {!r} from the object store ({}); recomputing".format(key, err))
                    replace = True
            value = generator(*args, **kwargs)
            try:
                self.store.put(value, key=entry, replace=replace,
                               meta={'property': key, 'identifier': ident if isinstance(ident, (str, int, float)) else repr(ident)})
            except (ValueError, TypeError, NotImplementedError) as err:
                warnings.warn("couldn't cache {!r} in the object store: {}".format(key, err))
            return value
        scoped = "{}@{}".format(key, chained_key(ident, key, options)[:32])
        return self.checkpointer.cached_eval(scoped, generator, args=args, kwargs=kwargs)
