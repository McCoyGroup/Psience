"""
The raw evaluation layer for electronic-structure models: MLIPs, semiempirical, and ab initio
methods, all behind one interface.

A backend evaluates properties (`'energy'`, `'charges'`, `'dipole'`, `'polarizability'`, plus
whatever a model exposes natively) and their analytic derivatives with respect to the nuclear
Cartesian coordinates, in the model's *native* units. It does no unit conversion, no
finite differences, no coordinate embedding, and no optimization: those live in
`Psience.Molecools.Evaluator`, whose evaluators wrap these backends.
"""

import abc
import contextlib
import importlib
import inspect
import os
import sys

import numpy as np

from McUtils.Data import AtomData
import McUtils.Numputils as nput

__all__ = [
    "ElectronicStructureBackend",
]


def import_or_explain(module, package, extra=None):
    """
    Import an optional dependency, raising an `ImportError` that says how to install it.
    """
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(
            f"{module!r} is required for this backend (install `{package}`)" + (f"; {extra}" if extra else "")
        ) from e


def resolve_torch_device(device=None):
    """
    The PyTorch device to run a model on: `device`, or `'cuda'` if available, else `'cpu'`.
    """
    if device is None:
        import torch
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    return device


def hf_authenticate(hf_token=None):
    """
    Make a HuggingFace token available to `huggingface_hub` for this process.

    `huggingface_hub.login` writes the token to `$HF_HOME/stored_tokens`, which races when many jobs share an `HF_HOME`;
    when a token is available we instead export it as `HF_TOKEN` (which `huggingface_hub` reads directly) and never touch the disk.
    An interactive `login()` is only used when there is no token at all.
    """
    if hf_token is None:
        hf_token = os.environ.get("HF_TOKEN")
    if hf_token is not None:
        os.environ["HF_TOKEN"] = hf_token
    else:
        from huggingface_hub import login
        login()


class quiet_mode:
    """
    Optionally silence stdout for the duration of a `with` block.
    """

    def __init__(self, quiet=True):
        self.quiet = quiet
        self._stdout = None
        self._devnull = None

    def __enter__(self):
        if self.quiet:
            self._stdout = sys.stdout
            self._devnull = open(os.devnull, 'w+').__enter__()
            sys.stdout = self._devnull

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.quiet:
            self._devnull.__exit__(exc_type, exc_val, exc_tb)
            sys.stdout = self._stdout


@contextlib.contextmanager
def environment_override(**variables):
    """
    Temporarily set (or, with `None`, leave alone) environment variables.
    """
    old = {}
    try:
        for k, v in variables.items():
            if v is None:
                continue
            old[k] = os.environ.get(k)
            os.environ[k] = v
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


class ElectronicStructureBackend(metaclass=abc.ABCMeta):
    """
    One electronic-structure model evaluated on Cartesian geometries of a fixed set of atoms.

    `evaluate(coords, order, properties)` returns, for each requested property, the list
    `[value, d value/dx, d^2 value/dx^2, ...]` up to `order`, where `x` are the flattened
    Cartesians (so the `k`-th term has shape `base_shape + (3 * natoms,) * k + property_shape`).
    Coordinates are in `distance_units`; property values are in `property_units[prop]`.
    Derivatives beyond `analytic_derivative_orders[prop]` are not available from a backend
    (the `Molecools` evaluators fill those in by finite differences).

    The total `charge` and spin `multiplicity` set at construction can be overridden per call.

    Subclasses implement `evaluate_properties` on a flat batch of structures.
    """

    distance_units = 'Angstroms'
    property_units = {}
    analytic_derivative_orders = {}
    property_shapes = {}  # per property, as a function of natoms: callables or tuples
    default_charge = 0
    default_multiplicity = None
    backend_registry = {}

    def __init__(self, atoms, charge=None, multiplicity=None, quiet=True):
        """
        :param atoms: atom symbols
        :type atoms: Sequence[str]
        :param charge: total charge (default `default_charge`)
        :type charge: int | None
        :param multiplicity: spin multiplicity (default `default_multiplicity`)
        :type multiplicity: int | None
        :param quiet: silence the underlying program's output
        :type quiet: bool
        """
        self.atoms = list(atoms)
        self.numbers = np.array([AtomData[a, "Number"] for a in self.atoms], dtype=int)
        self.charge = self.default_charge if charge is None else charge
        self.multiplicity = self.default_multiplicity if multiplicity is None else multiplicity
        self.quiet = quiet

    # region registry

    @classmethod
    def register(cls, name, backend=None):
        """
        Register a backend class under `name` (usable as a decorator).
        """
        if backend is None:
            def register(backend, name=name):
                return cls.register(name, backend)
            return register
        ElectronicStructureBackend.backend_registry[name] = backend
        return backend

    @classmethod
    def get_backends(cls):
        """
        The registered backends, by name.

        :rtype: dict[str, type]
        """
        return dict(ElectronicStructureBackend.backend_registry)

    @classmethod
    def handle_specialization(cls, tag):
        """
        Options implied by a `name:tag` specialization (e.g. `'mace:small'`).
        """
        raise NotImplementedError(f"{cls.__name__} does not handle specializations (got {tag!r})")

    @classmethod
    def resolve_backend(cls, spec):
        """
        Resolve `'name'` or `'name:tag'` to `(backend_class, options)`.

        :rtype: (type, dict)
        """
        tag = None
        if ':' in spec:
            spec, tag = spec.split(':', 1)
        backends = cls.get_backends()
        key = spec.lower()
        if key not in backends:
            raise ValueError(f"unknown electronic-structure backend {spec!r} (have {sorted(backends)})")
        backend = backends[key]
        opts = {} if tag is None else backend.handle_specialization(tag)
        return backend, opts

    @classmethod
    def from_spec(cls, spec, atoms, **opts):
        """
        Build a backend from a name/`'name:tag'` (or return `spec` if it already is one).

        :rtype: ElectronicStructureBackend
        """
        if isinstance(spec, ElectronicStructureBackend):
            return spec
        backend, base_opts = cls.resolve_backend(spec)
        return backend(atoms, **dict(base_opts, **opts))

    @classmethod
    def discover_dependencies(cls, **opts):
        """
        The packages/programs this backend needs, as a `McUtils.ExternalPrograms.JobDependencies`.
        """
        from McUtils.ExternalPrograms import JobDependencies
        return JobDependencies()

    @classmethod
    def get_init_options(cls):
        """
        The names this backend's constructor accepts and whether it also takes arbitrary
        keywords (which are then passed on to the underlying model).

        :rtype: (set[str], bool)
        """
        names = set()
        var_kw = False
        for klass in cls.__mro__:
            init = klass.__dict__.get('__init__')
            if init is None:
                continue
            for p in inspect.signature(init).parameters.values():
                if p.kind == p.VAR_KEYWORD:
                    var_kw = True
                elif p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY) and p.name != 'self':
                    names.add(p.name)
            if not var_kw:
                break  # only follow `**kwargs` chains up the MRO
        return names, var_kw

    # endregion

    # region properties

    @property
    def natoms(self):
        return len(self.atoms)

    @property
    def supported_properties(self):
        return tuple(self.analytic_derivative_orders.keys())

    def get_analytic_derivative_order(self, prop):
        """
        The highest coordinate-derivative order available analytically for `prop`.
        """
        try:
            return self.analytic_derivative_orders[prop]
        except KeyError:
            raise ValueError(f"{type(self).__name__} doesn't provide {prop!r} (has {list(self.supported_properties)})") from None

    def get_property_units(self, prop):
        """
        The native units of `prop` (a unit name, a tuple of unit names, or `None` for model-native/atomic).
        """
        return self.property_units.get(prop)

    def get_property_shape(self, prop):
        shape = self.property_shapes.get(prop, ())
        if callable(shape):
            shape = shape(self.natoms)
        return tuple(shape)

    def resolve_state(self, charge=None, multiplicity=None):
        return (
            self.charge if charge is None else charge,
            self.multiplicity if multiplicity is None else multiplicity
        )

    def check_electron_count(self, charge, multiplicity):
        """
        Raise if `charge` and `multiplicity` are inconsistent with the electron count.
        """
        mult = 1 if multiplicity is None else multiplicity
        n_electrons = int(self.numbers.sum()) - int(charge)
        if (n_electrons - (mult - 1)) % 2 != 0 or n_electrons < mult - 1:
            raise ValueError(
                f"charge {charge} and multiplicity {mult} are inconsistent with {n_electrons} electrons"
            )

    # endregion

    # region evaluation

    def prep_coords(self, coords):
        """
        Flatten `(..., natoms, 3)` or `(..., 3 natoms)` coordinates to `(nstruct, natoms, 3)`.

        :return: `(base_shape, flat_coords)`
        """
        coords = np.asanyarray(coords, dtype=float)
        n = self.natoms
        if coords.shape[-2:] == (n, 3):
            base_shape = coords.shape[:-2]
        elif coords.shape[-1] == 3 * n:
            base_shape = coords.shape[:-1]
        else:
            raise ValueError(f"coordinates of shape {coords.shape} don't match {n} atoms")
        return base_shape, coords.reshape((-1, n, 3))

    def canonicalize_request(self, order, properties):
        if properties is None:
            properties = self.supported_properties[:1]
        elif isinstance(properties, str):
            properties = (properties,)
        properties = tuple(properties)
        if nput.is_int(order):
            orders = {p: order for p in properties}
        else:
            orders = dict(order)
            properties = tuple(orders.keys())
        for p, o in orders.items():
            ad = self.get_analytic_derivative_order(p)
            if o > ad:
                raise ValueError(
                    f"{type(self).__name__} provides {p!r} derivatives through order {ad}, not {o}"
                )
        return orders

    @abc.abstractmethod
    def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
        """
        Evaluate properties on a flat batch of structures.

        :param coords: `(nstruct, natoms, 3)`, in `distance_units`
        :type coords: np.ndarray
        :param orders: the derivative order wanted for each property
        :type orders: dict[str, int]
        :param charge: the total charge
        :param multiplicity: the spin multiplicity
        :return: for each property, `[value, d1, ...]` with a leading `nstruct` axis and
            derivative axes of size `3 * natoms`
        :rtype: dict[str, list[np.ndarray]]
        """
        ...

    def evaluate(self, coords, order=0, properties=None, charge=None, multiplicity=None, **opts):
        """
        Evaluate properties and their analytic Cartesian derivatives.

        :param coords: `(..., natoms, 3)` or `(..., 3 natoms)`, in `distance_units`
        :type coords: np.ndarray
        :param order: the derivative order (for every property), or `{property: order}`
        :type order: int | dict[str, int]
        :param properties: property names (default: the backend's first supported property)
        :type properties: str | Iterable[str] | None
        :param charge: total charge override for this call
        :param multiplicity: spin multiplicity override for this call
        :param opts: backend-specific evaluation options
        :return: `{property: [value, d1, ...]}`, with the leading batch shape of `coords`
        :rtype: dict[str, list[np.ndarray]]
        """
        orders = self.canonicalize_request(order, properties)
        charge, multiplicity = self.resolve_state(charge, multiplicity)
        base_shape, flat = self.prep_coords(coords)
        res = self.evaluate_properties(flat, orders, charge, multiplicity, **opts)
        return {
            p: [np.asanyarray(t).reshape(base_shape + np.shape(t)[1:]) for t in res[p][:orders[p] + 1]]
            for p in orders
        }

    def evaluate_property(self, coords, prop, order=0, **opts):
        """
        `evaluate` for a single property.

        :rtype: list[np.ndarray]
        """
        return self.evaluate(coords, order=order, properties=(prop,), **opts)[prop]

    def evaluate_field_response(self, coords, coord_order=0, field_order=2, **opts):
        """
        The response to a uniform electric field, `[dipole_expansion, polarizability_expansion][:field_order]`,
        each a list over coordinate-derivative orders through `coord_order`. Backends with their own
        response models (charge equilibration) override this.

        :rtype: list[list[np.ndarray]]
        """
        props = ('dipole', 'polarizability')[:field_order]
        res = self.evaluate(coords, order=coord_order, properties=props, **opts)
        return [res[p] for p in props]

    # endregion

    def __repr__(self):
        return "{}(natoms={}, charge={}, multiplicity={})".format(
            type(self).__name__, self.natoms, self.charge, self.multiplicity
        )
