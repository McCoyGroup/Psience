from __future__ import annotations

import abc
import enum
import functools
import tempfile
import itertools
import math
import sys, os
import uuid
import time
import contextlib
import numpy as np
import warnings
import subprocess

from McUtils.Data import AtomData, UnitsData
from McUtils.Zachary import TensorDerivativeConverter, FiniteDifferenceDerivative, CoordinateFunction
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
import McUtils.Iterators as itut
from McUtils.Scaffolding import Logger, read_flat_tree
from McUtils.ExternalPrograms import SingularityLauncher

from ..Data import PotentialSurface, DipoleSurface
from .CoordinateSystems import MolecularEmbedding
from .Properties import NormalModesManager

__all__ = [
    "MolecularEvaluator",
    "EnergyEvaluator",
    "RDKitEnergyEvaluator",
    "AIMNet2EnergyEvaluator",
    "XTBEnergyEvaluator",
    "PySCFEnergyEvaluator",
    "PotentialFunctionEnergyEvaluator",
    "DipoleEvaluator",
    "DipoleFunctionDipoleEvaluator",
    "ReducedDimensionalPotentialHandler"
]

class MolecularEvaluator:
    def __init__(self, embedding:MolecularEmbedding, normal_modes:NormalModesManager):
        """
        **LLM Docstring**

        Store the coordinate embedding and normal-modes manager this evaluator will use to evaluate functions and build displacements.

        :param embedding: the molecular coordinate embedding to evaluate against
        :type embedding: MolecularEmbedding
        :param normal_modes: the normal-modes manager used for mode-based displacement/nearest-point queries
        :type normal_modes: NormalModesManager
        :return: None
        :rtype: None
        """
        self.embedding = embedding
        self.normal_modes = normal_modes
    @property
    def coords(self):
        """
        **LLM Docstring**

        The Cartesian coordinates, delegating to `self.embedding.coords`.

        :return: the Cartesian coordinates
        :rtype: CoordinateSet
        """
        return self.embedding.coords

    @property
    def masses(self):
        """
        **LLM Docstring**

        The atomic masses, delegating to `self.embedding.masses`.

        :return: the atomic masses
        :rtype: np.ndarray
        """
        return self.embedding.masses

    def evaluate(self,
                 func,
                 use_internals=None,
                 order=None,
                 strip_embedding=False
                 ):
        """
        **LLM Docstring**

        Evaluate a coordinate function (and, optionally, its derivatives) either directly in Cartesian coordinates or in internal coordinates with the result re-expanded back through the internals-by-Cartesians Jacobian, optionally stripping/restoring the fixed embedding coordinates around the internal-coordinate evaluation.

        :param func: the function to evaluate, expected to accept `order=` for derivative evaluation
        :type func: callable
        :param use_internals: whether to evaluate in internal coordinates; defaults to whether the embedding has internal coordinates defined
        :type use_internals: bool | None
        :param order: the highest derivative order to compute; if `None`, only the function value is returned
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates before evaluating in internal coordinates (and restore zeroed placeholders for them in the resulting derivative tensors)
        :type strip_embedding: bool
        :return: the function value, or `[value, deriv1, deriv2, ...]` if `order` is given
        :rtype: object | list
        """
        if use_internals is None:
            use_internals = self.embedding.internals is not None
        if use_internals:
            coords = self.embedding.internal_coordinates
            if strip_embedding:
                embedding_coords = self.embedding.embedding_coords
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3 * len(self.embedding.masses)), embedding_coords)
                    coords = coords.reshape(coords.shape[:-2] + (3 * len(self.embedding.masses),))
                    coords = coords[..., good_coords]
            if order is None:
                return func(coords).view(np.ndarray)

            terms = func(coords, order=order)
            # raise Exception([t.shape for t in terms], coords.shape)
            if strip_embedding:
                embedding_coords = self.embedding.embedding_coords
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3 * len(self.embedding.masses)), embedding_coords)

                    const = terms[0]
                    terms = terms[1:]
                    new = []
                    ncs = 3 * len(self.embedding.masses)
                    if self.embedding.coords.ndim > 2:
                        npts = len(coords)
                    else:
                        npts = 1
                    for n, ders in enumerate(terms):
                        dt = np.zeros((npts,) + (ncs,)*(n+1))
                        idx_pos = (...,) + np.ix_(*[good_coords]*(n+1))
                        dt[idx_pos] = ders.view(np.ndarray)
                        if self.embedding.coords.ndim == 2:
                            dt = dt[0]
                        new.append(dt)
                    terms = [const] + new

            const = terms[0]
            jacs = self.embedding.get_internals_by_cartesians(order)

            terms = TensorDerivativeConverter(
                jacs,
                terms[1:],
                jacobians_name='dXdR',
                values_name='f'
            ).convert()  # , check_arrays=True)

            return [const.view(np.ndarray)] + [t.view(np.ndarray) for t in terms]
        else:
            if order is None:
                return func(self.embedding.coords).view(np.ndarray)
            else:
                return [x.view(np.ndarray) for x in func(self.embedding.coords, order=order)]

    def evaluate_at(self,
                    func,
                    coords,
                    use_internals=None,
                    order=None,
                    strip_embedding=False
                    ):
        """
        **LLM Docstring**

        Evaluate a coordinate function at an alternate set of coordinates, by building a fresh `MolecularEvaluator` for those coordinates (keeping the same masses/internal-coordinate spec/normal modes) and calling `evaluate` on it.

        :param func: the function to evaluate
        :type func: callable
        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param use_internals: whether to evaluate in internal coordinates
        :type use_internals: bool | None
        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates around the internal-coordinate evaluation
        :type strip_embedding: bool
        :return: the function value, or the value/derivative expansion if `order` is given
        :rtype: object | list
        """
        return type(self)(
            MolecularEmbedding(self.embedding.masses, coords, self.embedding.internals),
            self.normal_modes
        ).evaluate(func, use_internals=use_internals, order=order, strip_embedding=strip_embedding)

    def get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None,
                                  use_internals:"bool|Literal['reembed', 'convert']"=False,
                                  coordinate_expansion=None,
                                  expansion_active_positions=None,
                                  strip_embedding=False,
                                  shift=True,
                                  coords=None):
        """
        **LLM Docstring**

        Build displaced coordinates from a reference geometry by applying `displacements` -- either directly to selected flat coordinate indices/atom-axis pairs, to the full coordinate array, or (if `coordinate_expansion` is given) as a Taylor-series expansion along an arbitrary coordinate-transformation direction (e.g. a normal mode) -- optionally stripping the fixed embedding coordinates and handling the resulting internal-coordinate re-embedding/conversion convention.

        :param displacements: the displacement value(s) to apply
        :type displacements: np.ndarray
        :param which: explicit flat coordinate indices (or `(atom, axis)` pairs) to displace; mutually exclusive with `sel`/`axes`
        :type which: Iterable | None
        :param sel: a selection of atoms to restrict the displacement to (used together with `axes`)
        :type sel: Iterable[int] | None
        :param axes: which Cartesian axes to displace along (used together with `sel`)
        :type axes: Iterable[int] | None
        :param use_internals: whether the displacement (and/or resulting coordinates) should be interpreted/handled in internal coordinates; `'reembed'`/`'convert'` select how the internal-coordinate result is converted back
        :type use_internals: bool | str
        :param coordinate_expansion: a Taylor-series coordinate-transformation expansion (e.g. a normal-mode Jacobian) to displace along instead of raw coordinates
        :type coordinate_expansion: list[np.ndarray] | None
        :param expansion_active_positions: restrict `coordinate_expansion`'s output to these positions before applying
        :type expansion_active_positions: Iterable[int] | None
        :param strip_embedding: whether the coordinate expansion/target coordinates have had the fixed embedding coordinates stripped
        :type strip_embedding: bool
        :param shift: whether `displacements` are relative shifts (`True`) added to the reference, or absolute target values (`False`)
        :type shift: bool
        :param coords: alternate reference coordinates to displace from, instead of `self.embedding.coords`/`internal_coordinates`
        :type coords: np.ndarray | None
        :return: the displaced coordinates
        :rtype: CoordinateSet
        :raises ValueError: if `use_internals` is requested but the embedding has no internal-coordinate spec, or if `displacements` doesn't match the reference coordinates' shape
        """
        displacements = np.asanyarray(displacements)

        if which is not None:
            which = tuple(
                np.ravel_multi_index(idx, (len(self.embedding.masses), 3))
                    if not nput.is_int(idx) else
                idx
                for idx in which
            )

        if use_internals and self.embedding.internals is None:
            raise ValueError("can't displace in internals without internal coordinate spec")
        if coords is None:
            base_shape = self.embedding.coords.shape[:-2]
            base_coords = self.embedding.coords if not use_internals else self.embedding.internal_coordinates
        else:
            coords = np.asanyarray(coords)
            base_shape = coords.shape[:-2]
            if not use_internals:
                base_coords = coords
            else:
                base_coords = self.embedding.get_internals(coords=coords, strip_embedding=False)
        base_dim = len(base_shape)



        if coordinate_expansion is not None:
            if which is not None:
                base_disps = np.zeros(displacements.shape[:-1] + coordinate_expansion[0].shape[:1])
                base_disps[..., which] = displacements
                displacements = base_disps
            _ = []
            new_disps = 0
            shared = displacements.ndim-1
            for n,disp in enumerate(coordinate_expansion):
                disp = disp / math.factorial(n+1)
                if nput.is_numeric(disp) and disp == 0: continue
                for i in range(n+1):
                    if i == 0:
                        disp = np.tensordot(displacements, disp, axes=[-1, 0])
                    else:
                        disp = nput.vec_tensordot(displacements, disp, axes=[-1, shared], shared=shared)
                new_disps = new_disps + disp

            if expansion_active_positions is not None:
                new_disps = new_disps[..., expansion_active_positions]
                if strip_embedding:
                    ecs = self.embedding.embedding_coords
                    ntot = np.prod(base_coords.shape[base_dim:], dtype=int)
                    if len(ecs) > 0 and new_disps.shape[-1] < ntot:
                        rem_coords = np.setdiff1d(np.arange(ntot), ecs)
                        which = rem_coords[expansion_active_positions,]
                    else:
                        which = expansion_active_positions
                else:
                    which = expansion_active_positions
            elif strip_embedding:
                ecs = self.embedding.embedding_coords
                ntot = np.prod(base_coords.shape[base_dim:], dtype=int)
                if len(ecs) > 0 and new_disps.shape[-1] < ntot:
                    rem_coords = np.setdiff1d(np.arange(ntot), ecs)
                    which = rem_coords
                else:
                    which = None
            else:
                which = None

            if which is None:
                displacements = new_disps.reshape(
                    new_disps.shape[:-1] + base_coords.shape[base_dim:]
                )
            else:
                displacements = new_disps
        elif strip_embedding:
            ecs = self.embedding.embedding_coords
            all_coords = np.arange(len(self.embedding.masses) * 3)
            which = np.setdiff1d(all_coords, ecs)[which,]

        if which is not None:
            if displacements.shape[-1] != len(which):  # displacements provided in atom coordinates
                displacements = displacements.reshape(
                    displacements.shape[:-2] +
                    (np.prod(displacements.shape[-2:], dtype=int),)
                )
            # pad displacements so we can add them together
            if displacements.ndim > 1:
                for _ in range(displacements.ndim - 1):
                    base_coords = np.expand_dims(base_coords, 0)
                base_coords = np.broadcast_to(base_coords, displacements.shape[:-1] + base_coords.shape[-2:])
            base_coords = base_coords.copy()
            flat_coords = base_coords.reshape(displacements.shape[:-1] + (-1,))

            if shift:
                flat_coords[..., which] += displacements
            else:
                flat_coords[..., which] = displacements
            base_coords = flat_coords.reshape(base_coords.shape)

        elif sel is not None or axes is not None:
            if displacements.ndim > 2:
                for _ in range(displacements.ndim - 2):
                    base_coords = np.expand_dims(base_coords, 0)
                base_coords = np.broadcast_to(base_coords, displacements.shape[:-2] + base_coords.shape[-2:])
            base_coords = base_coords.copy()

            if sel is None:
                sel = np.arange(len(self.embedding.masses))
            if axes is None:
                axes = np.arange(3)

            sel = np.asanyarray(sel)[:, np.newaxis]
            axes = np.asanyarray(axes)[np.newaxis, :]

            if shift:
                base_coords[..., sel, axes] += displacements
            else:
                base_coords[..., sel, axes] = displacements

        else:
            if displacements.shape[-base_coords.ndim:] != base_coords.shape:
                raise ValueError("displacements with shape {} passed but coordinates have shape {}".format(
                    displacements.shape,
                    base_coords.shape
                ))
            if shift:
                bnd = base_coords.ndim
                base_coords = np.expand_dims(base_coords, list(range(displacements.ndim - bnd)))
                base_coords = base_coords + displacements
            else:
                base_coords = displacements

        if use_internals:
            # track the embedding info...
            base_coords = self.embedding.internal_coordinates.system(
                base_coords,
                **self.embedding.internal_coordinates.converter_options
            )
            if isinstance(use_internals, str):
                if use_internals == 'convert':
                    base_coords = base_coords.convert(self.embedding.coords.system)
                elif use_internals == 'reembed':
                    base_coords = self.embedding.embed_coords(
                        base_coords.convert(self.embedding.coords.system),
                        proper_rotation=True
                    )
        else:
            base_coords = self.embedding.coords.system(base_coords)
        return base_coords

    def get_scan_coordinates(self,
                             domains,
                             internals=False,
                             which=None, sel=None, axes=None,
                             coordinate_expansion=None,
                             strip_embedding=False,
                             shift=True,
                             return_displacements=False
                             ):
        """
        **LLM Docstring**

        Build a grid of displaced coordinates spanning the given `(start, stop, num)` domains for each scanned coordinate, via `get_displaced_coordinates` applied to a meshgrid of displacement values.

        :param domains: the `(start, stop, num)` ranges to scan over, one per scanned coordinate
        :type domains: Iterable[tuple]
        :param internals: whether the scan coordinates/displacements are internal coordinates
        :type internals: bool
        :param which: explicit coordinate indices to scan; forwarded to `get_displaced_coordinates`
        :type which: Iterable | None
        :param sel: atom selection to scan over; forwarded to `get_displaced_coordinates`
        :type sel: Iterable[int] | None
        :param axes: axes to scan over; forwarded to `get_displaced_coordinates`
        :type axes: Iterable[int] | None
        :param coordinate_expansion: a coordinate-transformation expansion to scan along instead of raw coordinates
        :type coordinate_expansion: list[np.ndarray] | None
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :param shift: whether the scan values are relative shifts or absolute targets
        :type shift: bool
        :param return_displacements: whether to also return the raw displacement mesh used
        :type return_displacements: bool
        :return: the scan coordinates, or `(displacement_mesh, scan_coordinates)` if `return_displacements` is set
        :rtype: np.ndarray | tuple
        """

        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        disps = self.get_displaced_coordinates(displacement_mesh, shift=shift,
                                              use_internals=internals, which=which, sel=sel, axes=axes,
                                              coordinate_expansion=coordinate_expansion,
                                              strip_embedding=strip_embedding
                                              )
        if displacement_mesh.shape[-1] == 1 and disps.shape[1] == 1:
            disps = disps.reshape(disps.shape[:1] + disps.shape[2:])
        if return_displacements:
            return displacement_mesh, disps
        else:
            return disps

    def get_nearest_displacement_atoms(self,
                                       points,
                                       sel=None, axes=None, weighting_function=None,
                                       return_distances=False
                                       ):
        """
        **LLM Docstring**

        For each of a set of query points, find the atom (restricted to `sel`, along `axes`) whose (optionally mass-weighted) distance to the point is smallest.

        :param points: the query points to find nearest atoms for
        :type points: np.ndarray
        :param sel: a selection of atoms to restrict the search to; defaults to all atoms
        :type sel: Iterable[int] | None
        :param axes: which Cartesian axes to compute distances over; defaults to all three
        :type axes: Iterable[int] | None
        :param weighting_function: a function of the atomic masses used to weight distances (e.g. to prefer displacing lighter atoms); defaults to `np.sqrt`
        :type weighting_function: callable | None
        :param return_distances: whether to also return the (weighted) distances to the nearest atoms
        :type return_distances: bool
        :return: the nearest atom index for each point, or `(atom_idx, distances)` if `return_distances` is set
        :rtype: np.ndarray | tuple
        """

        pts = np.asanyarray(points)
        smol = pts.ndim == 1
        if smol: pts = pts[np.newaxis]
        pts = pts.reshape(-1, pts.shape[-1])

        if axes is None:
            axes = [0, 1, 2]
        axes = np.asanyarray(axes)

        if sel is None:
            sel = np.arange(len(self.masses))
        sel = np.asanyarray(sel)

        ref = self.coords
        masses = self.masses
        dists = np.linalg.norm(
            pts[:, np.newaxis, :] - ref[np.newaxis, sel[:, np.newaxis], axes[np.newaxis, :]],
            axis=-1
        )
        if weighting_function is None:
            weighting_function = np.sqrt
        dists = dists * weighting_function(masses)[np.newaxis, sel]
        nearest = np.argmin(dists, axis=-1)
        atom_idx = sel[nearest]

        if return_distances:
            return atom_idx, dists[np.arange(len(nearest)), nearest]
        else:
            return atom_idx

    def get_nearest_displacement_coordinates(self,
                                             points,
                                             sel=None, axes=None, weighting_function=None,
                                             modes_nearest=False,
                                             return_distances=False
                                             ):
        """
        Displaces the _nearest_ atom (in a mass-weighted sense) to the given point
        This allows for the development of functions / potentials that are in the small-displacement limit
        where everything _except_ the nearest atom to the given Cartesian position is at its equilibrium value

        :param points:
        :param sel:
        :param axes:
        :param weighting_function:
        :return:
        """

        pts = np.asanyarray(points)
        smol = pts.ndim == 1
        if smol: pts = pts[np.newaxis]
        base_shape = pts.shape[:-1]
        pts = pts.reshape(-1, pts.shape[-1])

        if axes is None:
            axes = [0, 1, 2]
        axes = np.asanyarray(axes)

        if sel is None:
            sel = np.arange(len(self.masses))
        sel = np.asanyarray(sel)

        ref = self.coords
        atom_idx = self.get_nearest_displacement_atoms(
            points,
            sel=sel, axes=axes, weighting_function=weighting_function,
            return_distances=return_distances
        )
        if return_distances:
            atom_idx, dists = atom_idx

        if modes_nearest:
            npts = len(points)
            diag = np.arange(npts)
            mode_origin = self.normal_modes.modes.basis.origin
            mode_origin = mode_origin[:, axes]
            origin_coords = np.broadcast_to(mode_origin[np.newaxis], (npts,) + mode_origin.shape)

            origin_coords = origin_coords[diag, atom_idx]
            nearest_disps = pts - origin_coords

            modes = self.normal_modes.modes.basis.matrix
            mat = np.reshape(modes, (1, len(self.masses), -1, modes.shape[-1]))
            mat = mat[:, :, axes, :]
            mat = np.broadcast_to(mat, (npts,) + mat.shape[1:])
            mat_bits = mat[diag, atom_idx, :, :]  # N x 3 x N_modes
            pseudo_inverses = np.linalg.inv(mat_bits @ mat_bits.transpose(0, 2, 1))
            # print(mat_bits.shape, pseudo_inverses.shape, nearest_disps.shape)
            ls_coords = mat_bits.transpose(0, 2, 1) @ pseudo_inverses @ nearest_disps[:, :, np.newaxis]
            # ls_coords = np.reshape(ls_coords, ls_coords.shape[:2])
            # print(modes.shape, ls_coords.shape, pseudo_inverses.shape, modes.shape)
            test = modes[np.newaxis, :, :] @ ls_coords
            ls_coords = np.reshape(ls_coords, ls_coords.shape[:2])
            # print(test.shape, nearest_disps.shape)
            # print(test.reshape(test.shape[0], -1, 3))
            # print(nearest_disps)
            # print(pts)
            # print(ls_coords)
            norms = np.linalg.norm(ls_coords, axis=-1)
            sorting = np.argsort(norms)
            # print(norms[sorting[:5]])
            # print(ls_coords[sorting[:5]])
            # print(pts[sorting[:5]])


            # raise Exception(ls_coords)
            coords = self.normal_modes.modes.basis.to_new_modes().unembed_coords(
                np.reshape(ls_coords, ls_coords.shape[:2])
            )

            # print(coords)

        else:
            coords = np.broadcast_to(ref[np.newaxis], (len(pts),) + ref.shape).copy()
            coords[np.arange(len(atom_idx))[:, np.newaxis], atom_idx[:, np.newaxis], axes[np.newaxis, :]] = pts
            coords = coords.reshape(base_shape + coords.shape[-2:])
        if smol: coords = coords[0]

        if return_distances:
            return coords, dists
        else:
            return coords

    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        """
        **LLM Docstring**

        Build a grid of "nearest displacement" coordinates (see `get_nearest_displacement_coordinates`) spanning the given `(start, stop, num)` domains, via a meshgrid of query points.

        :param domains: the `(start, stop, num)` ranges to scan over
        :type domains: Iterable[tuple]
        :param sel: a selection of atoms to restrict the nearest-atom search to
        :type sel: Iterable[int] | None
        :param axes: which Cartesian axes the domains correspond to
        :type axes: Iterable[int] | None
        :return: the grid of nearest-displacement coordinates
        :rtype: np.ndarray
        """
        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        return self.get_nearest_displacement_coordinates(displacement_mesh, axes=axes, sel=sel)

class InternalHandlingMode(enum.Enum):
    ReembedCartesians = "reembed"
    Convert = "convert"
    Cartesians = "cartesians"

    @classmethod
    def resolve(cls, mode):
        """
        **LLM Docstring**

        Normalize a user-supplied internal-coordinate handling flag into an `InternalHandlingMode` enum value: `True` maps to `Convert`, `False` maps to `Cartesians`, and any other value is passed through to the enum constructor (e.g. an already-valid mode string).

        :param mode: the raw handling-mode flag to resolve
        :type mode: bool | str | InternalHandlingMode
        :return: the resolved enum member
        :rtype: InternalHandlingMode
        """
        if mode is True:
            return cls.Convert
        elif mode is False:
            return cls.Cartesians
        else:
            return cls(mode)

class PropertyEvaluator(metaclass=abc.ABCMeta):

    def __init__(self,
                 embedding=None,
                 permutation=None,
                 use_internals=True,
                 strip_embedding=True,
                 flatten_internals=True,
                 reembed_cartesians=False,
                 supports_internals=None,
                 **defaults):
        """
        **LLM Docstring**

        Set up the shared configuration for evaluating a molecular property: which coordinate embedding to use, whether/how to work in internal coordinates (falling back to `'reembed'` if the concrete evaluator doesn't support internals directly), whether to strip fixed embedding coordinates, whether to flatten internal-coordinate arrays, whether to Eckart-reembed Cartesians, an atom permutation to apply, and any evaluator-specific default keyword options.

        :param embedding: the coordinate embedding to evaluate the property against
        :type embedding: MolecularEmbedding | None
        :param permutation: an atom (or flattened-coordinate) permutation to apply when embedding/unembedding coordinates and derivatives
        :type permutation: np.ndarray | None
        :param use_internals: whether to evaluate in internal coordinates; if the concrete evaluator subclass doesn't support internals directly, this is coerced to `'reembed'`
        :type use_internals: bool | str
        :param strip_embedding: whether to strip the fixed embedding coordinates when converting to/from internal coordinates
        :type strip_embedding: bool
        :param flatten_internals: whether to flatten internal-coordinate arrays to a single trailing axis
        :type flatten_internals: bool
        :param reembed_cartesians: whether to Eckart-reembed Cartesian coordinates before evaluating
        :type reembed_cartesians: bool
        :param supports_internals: whether this evaluator can evaluate directly in internal coordinates; defaults to the class attribute `self.supports_internals`
        :type supports_internals: bool | None
        :param defaults: extra default keyword options passed through to `evaluate_expansion`/`evaluate_term`
        :type defaults: dict
        :return: None
        :rtype: None
        """
        self.defaults = defaults
        self.embedding = embedding
        if supports_internals is None:
            supports_internals = self.supports_internals
        self.supports_internals = supports_internals
        if use_internals and not self.supports_internals:
            use_internals = 'reembed'
        self.use_internals = use_internals
        self.strip_embedding = strip_embedding
        self.flatten_internals = flatten_internals
        self.reembed_cartesians = reembed_cartesians
        self.permutation = permutation

    def use_internal_coordinate_handlers(self):
        """
        **LLM Docstring**

        Whether this evaluator should route coordinate handling through the internal-coordinate machinery: true when `use_internals` is enabled, an embedding is set, and that embedding actually has internal coordinates defined.

        :return: whether internal-coordinate handling applies
        :rtype: bool
        """
        return (self.use_internals
                and self.embedding is not None
                and self.embedding.internals is not None)

    def embed_coords(self, coords, embed_reembedded=True):
        """
        **LLM Docstring**

        Convert raw Cartesian coordinates into whatever representation this evaluator actually evaluates in: internal coordinates (optionally stripped of fixed embedding coordinates and flattened) if `use_internals` is enabled and the embedding has internals, or Eckart-reembedded Cartesians if `reembed_cartesians` is set; then applies any configured atom/coordinate permutation.

        :param coords: the raw Cartesian coordinates to embed
        :type coords: np.ndarray
        :param embed_reembedded: whether to still perform the internal-coordinate conversion even when `use_internals` is specifically `'reembed'` (which otherwise defers embedding to elsewhere)
        :type embed_reembedded: bool
        :return: the embedded (and possibly permuted) coordinates, in whatever representation the evaluator uses
        :rtype: np.ndarray
        :raises NotImplementedError: if permuting non-flattened internal coordinates is requested
        """
        coords = np.asanyarray(coords)
        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                base_coords = coords
                coords = self.embedding.get_internals(coords=coords, strip_embedding=self.strip_embedding)
                if not self.strip_embedding and self.flatten_internals:
                    if base_coords.ndim == coords.ndim:
                        coords = coords.reshape(base_coords.shape[:-2] + (-1,))
            elif self.reembed_cartesians:
                coords = self.embedding.embed_coords(coords, proper_rotation=True)
        if self.permutation is not None:
            if self.use_internals and self.embedding.internals is not None:
                if not self.strip_embedding and not self.flatten_internals:
                    raise NotImplementedError("don't know how to permute not-flat internals...")
                coords = coords[..., self.permutation]
            else:
                coords = coords[..., self.permutation, :]
        return coords

    def unembed_coords(self, coords, embed_reembedded=True):
        """
        **LLM Docstring**

        Undo `embed_coords`: convert the evaluator's working coordinate representation back to plain Cartesians (via `self.embedding.get_cartesians` if internal coordinates were used), and undo any applied atom/coordinate permutation.

        :param coords: the evaluator's working-representation coordinates to unembed
        :type coords: np.ndarray
        :param embed_reembedded: whether to perform the internal-to-Cartesian conversion even when `use_internals` is specifically `'reembed'`
        :type embed_reembedded: bool
        :return: the Cartesian coordinates
        :rtype: np.ndarray
        :raises NotImplementedError: if un-permuting non-flattened internal coordinates is requested
        """
        # coords = np.asanyarray(coords)
        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                coords = self.embedding.get_cartesians(coords=coords, strip_embedding=self.strip_embedding)

        if self.permutation is not None:
            inv = np.argsort(self.permutation)
            if self.use_internals and self.embedding.internals is not None:
                if not self.strip_embedding and not self.flatten_internals:
                    raise NotImplementedError("don't know how to permute not-flat internals...")
                coords = coords[..., inv]
            else:
                coords = coords[..., inv, :]
        return coords

    def embed_derivs(self, coords, derivs, embed_reembedded=True):
        """
        **LLM Docstring**

        Re-express a set of derivative tensors (computed with respect to the evaluator's working coordinates) in terms of Cartesian coordinates, by re-expanding through the Cartesians-by-internals Jacobian if internal coordinates are in use.

        :param coords: the coordinates the derivatives were computed at (in the evaluator's working representation)
        :type coords: np.ndarray
        :param derivs: the derivative tensors to re-express
        :type derivs: list[np.ndarray]
        :param embed_reembedded: whether to perform the re-expansion even when `use_internals` is specifically `'reembed'`
        :type embed_reembedded: bool
        :return: the derivative tensors, re-expressed in Cartesian coordinates if applicable, or unchanged
        :rtype: list[np.ndarray]
        """
        if len(derivs) == 0: return derivs

        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                expansion = self.embedding.get_cartesians_by_internals(
                    order=len(derivs),
                    coords=coords,
                    strip_embedding=self.strip_embedding
                )
                derivs = nput.tensor_reexpand(expansion, derivs, len(derivs))

        return derivs

    def unembed_derivs(self, base_coords, coords, derivs, embed_reembedded=True):
        """
        **LLM Docstring**

        Undo `embed_derivs` (and any atom/coordinate permutation) on a set of Cartesian-coordinate derivative tensors, converting them back into the evaluator's working coordinate representation (internal coordinates, or Eckart-unrotated Cartesians) as appropriate.

        :param base_coords: the reference coordinates (in the evaluator's working representation) the returned derivatives should be expressed relative to
        :type base_coords: np.ndarray
        :param coords: the Cartesian coordinates the input derivatives were computed at
        :type coords: np.ndarray
        :param derivs: the Cartesian-coordinate derivative tensors to convert
        :type derivs: list[np.ndarray]
        :param embed_reembedded: whether to perform the re-expansion even when `use_internals` is specifically `'reembed'`
        :type embed_reembedded: bool
        :return: the derivative tensors, converted to the evaluator's working representation
        :rtype: list[np.ndarray]
        """
        if len(derivs) == 0: return derivs

        if self.permutation is not None:
            inv = np.argsort(self.permutation)
            new_derivs = []
            base_shape = base_coords.shape[:-2]
            ndim = len(base_shape)
            f_shape = derivs[0].shape[ndim + 1:]
            fdim = len(f_shape)
            ncoord = derivs[0].shape[ndim]
            nat = ncoord // 3  # not always used, but never an error
            for i, d in enumerate(derivs):
                # TODO: fix assumption that we always request all orders
                if self.use_internals and self.embedding.internals is not None:
                    for ax in range(-(i + 1), 0):
                        d = np.take(d, inv, axis=(ax - fdim))
                else:
                    for ax in range(i + 1):
                        coord_shape = tuple(ncoord for _ in range(ax)) + (nat, 3) + tuple(ncoord for _ in range(i - ax))
                        d_shape = d.shape
                        d = d.reshape(base_shape + coord_shape + f_shape)
                        d = np.take(d, inv, axis=ndim + ax).reshape(d_shape)
                new_derivs.append(d)
            derivs = new_derivs

        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                expansion = self.embedding.get_internals_by_cartesians(
                    order=len(derivs),
                    coords=base_coords,
                    strip_embedding=self.strip_embedding
                )
                derivs = nput.tensor_reexpand(expansion, derivs, len(derivs))
            elif self.reembed_cartesians:
                derivs = self.embedding.unembed_derivs(coords, derivs)
        return derivs

    @abc.abstractmethod
    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Abstract hook for evaluating the property's value (or a specific analytic derivative order) at the given coordinates. Concrete evaluator subclasses must implement this.

        :param coords: the coordinates to evaluate at, in the evaluator's working representation
        :type coords: np.ndarray
        :param order: the analytic derivative order to evaluate
        :type order: int
        :param opts: extra evaluator-specific options
        :type opts: dict
        :return: the evaluated term(s)
        :rtype: object
        """
        ...

    @classmethod
    def prep_mol_opts(cls, mol, embedding=None, charge=None, multiplicity=None, **opts):
        """
        **LLM Docstring**

        Build the keyword-argument dict used to construct an evaluator from a molecule, pulling the embedding, charge, and spin/multiplicity off the molecule wherever not explicitly overridden.

        :param mol: the molecule to build evaluator options from
        :type mol: AbstractMolecule
        :param embedding: an explicit embedding to use instead of `mol.embedding`
        :type embedding: MolecularEmbedding | None
        :param charge: an explicit charge to use instead of `mol.charge`
        :type charge: int | None
        :param multiplicity: an explicit spin multiplicity to use instead of `mol.spin`
        :type multiplicity: object | None
        :param opts: additional options passed through unchanged
        :type opts: dict
        :return: the resolved options dict
        :rtype: dict
        """
        if embedding is None: embedding = mol.embedding
        if embedding is not None: opts['embedding'] = embedding
        if charge is None: charge = mol.charge
        if charge is not None: opts['charge'] = charge
        if multiplicity is None: multiplicity = mol.spin
        if multiplicity is not None: opts['multiplicity'] = multiplicity
        return opts
    @classmethod
    @abc.abstractmethod
    def from_mol(cls, mol, **opts):
        """
        **LLM Docstring**

        Abstract constructor for building an evaluator instance from a molecule. Concrete evaluator subclasses must implement this.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param opts: evaluator-specific construction options
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PropertyEvaluator
        """
        ...

    fd_defaults=dict(
        stencil=None,
        mesh_spacing=.005,
        displacement_function=None,
        prep=None,
        lazy=False,
        cache_evaluations=True,
        parallelizer=None,
    )
    def get_fd_opts(self, **opts):
        """
        **LLM Docstring**

        Merge the class-level `fd_defaults` with any explicitly passed finite-difference options.

        :param opts: explicit per-call finite-difference option overrides
        :type opts: dict
        :return: the merged finite-difference options
        :rtype: dict
        """
        return dict(self.fd_defaults, **opts)

    def finite_difference_derivs(self, coords, order,
                                 batched_orders=None,
                                 displacement_generator=None,
                                 coordinate_prep=None,
                                 index_filter=None,
                                 analytic_derivative_order=None,
                                 **opts):
        """
        **LLM Docstring**

        Compute derivative tensors of the property numerically via finite differences in the evaluator's working (typically Cartesian) coordinates, wrapping `evaluate_term` as the base function fed to a `FiniteDifferenceDerivative`.

        :param coords: the reference coordinates to differentiate around
        :type coords: np.ndarray
        :param order: the highest total derivative order needed
        :type order: int
        :param batched_orders: whether `evaluate_term` returns all derivative orders in one batched call rather than being called once per order
        :type batched_orders: bool | None
        :param displacement_generator: a callable to transform the raw finite-difference displacement points before evaluation (e.g. to apply a coordinate expansion)
        :type displacement_generator: callable | None
        :param coordinate_prep: a callable to preprocess the flattened center coordinates before finite-differencing
        :type coordinate_prep: callable | None
        :param index_filter: a callable/filter restricting which derivative-tensor index combinations are computed
        :type index_filter: callable | None
        :param analytic_derivative_order: the derivative order already available analytically (via `evaluate_term`), used as the finite-difference base point rather than order 0
        :type analytic_derivative_order: int | None
        :param opts: extra options, split between `FiniteDifferenceDerivative` construction options and options forwarded to `evaluate_term`
        :type opts: dict
        :return: the finite-difference-computed derivative tensors, from `analytic_derivative_order + 1` up to `order`
        :rtype: list[np.ndarray]
        :raises NotImplementedError: if `multi_expansion_order > 0` (multi-term expansion unsplitting isn't supported yet)
        """
        opts = dev.OptionsSet(opts)
        fd_opts = opts.filter(FiniteDifferenceDerivative)
        opts = opts.exclude(FiniteDifferenceDerivative)

        if batched_orders is None:
            batched_orders = self.batched_orders

        base_shape = coords.shape[:-2]
        coord_shape = coords.shape[-2:]
        flat_coords = coords.reshape(-1, np.prod(coord_shape, dtype=int))

        if analytic_derivative_order is None:
            analytic_derivative_order = self.analytic_derivative_order

        expansion_shape = [None]
        def derivs(structs, center=flat_coords, expansion_shape=expansion_shape):
            """
            **LLM Docstring**

            Evaluate `evaluate_term` at a batch of displaced structures (after applying `displacement_generator`, if given) for use as the base function inside the enclosing finite-difference derivative computation.

            :param structs: the (possibly displacement-generator-transformed) flattened coordinate points to evaluate at
            :type structs: np.ndarray
            :param center: the reference (un-displaced) flattened coordinates, from the enclosing scope
            :type center: np.ndarray
            :param expansion_shape: a mutable single-element list used to record the shape of a multi-term expansion, from the enclosing scope
            :type expansion_shape: list
            :return: the evaluated term(s) at the given structures
            :rtype: np.ndarray
            """
            if displacement_generator is not None:
                structs = displacement_generator(structs, evaluator=self, center=center)
            reconst = structs.reshape(structs.shape[:-1] + coord_shape)
            ders = self.evaluate_term(reconst, analytic_derivative_order, **opts)

            if self.multi_expansion_order > 0:
                nrec = reconst.shape[:-2]
                if expansion_shape[0] is not None:
                    expansion_shape[0] = [d[-1].shape for d in ders]
                ders = [
                    np.concatenate([d[-1].reshape(d[-1].shape[:nrec] + (-1,))], axis=-1)
                    for d in ders
                ]

            if batched_orders:
                return ders[-1]
            else:
                return ders

        der = FiniteDifferenceDerivative(derivs,
                                         function_shape=((0,), (0,) * analytic_derivative_order),
                                         **self.get_fd_opts(**fd_opts)
                                         )
        if coordinate_prep is not None:
            flat_coords = coordinate_prep(flat_coords, evaluator=self)
        tensors = der.derivatives(flat_coords).derivative_tensor(
            list(range(1, (order - analytic_derivative_order) + 1)),
            pos_filter=index_filter
        )
        if flat_coords.shape[0] > 1:
            tensors = [np.moveaxis(d, n + 1, 0) for n, d in enumerate(tensors)]

        res = [
            (
                t.reshape(base_shape + t.shape[1:])
                    if flat_coords.shape[0] != 1 else
                t.reshape(base_shape + t.shape)
            )
                if t.ndim > 0 else
            t
            for t in tensors
        ]

        if self.multi_expansion_order > 0:
            #TODO: un-split expansions
            ...
            raise NotImplementedError("unsplitting not yet supported")

        return res

    def internal_finite_difference_derivs(self,
                                          coords,
                                          order,
                                          batched_orders=None,
                                          displacement_generator=None,
                                          coordinate_prep=None,
                                          index_filter=None,
                                          analytic_derivative_order=None,
                                          **opts):
        """
        **LLM Docstring**

        Compute derivative tensors of the property numerically via finite differences directly in internal coordinates, converting each finite-difference-displaced internal-coordinate point back to Cartesians (and re-embedding any analytic derivatives) before calling `evaluate_term`, when `use_internals` is set to `'reembed'`.

        :param coords: the reference coordinates (internal or Cartesian, converted to internal coordinates if needed) to differentiate around
        :type coords: np.ndarray
        :param order: the highest total derivative order needed
        :type order: int
        :param batched_orders: whether `evaluate_term` returns all derivative orders in one batched call
        :type batched_orders: bool | None
        :param displacement_generator: a callable to transform the raw finite-difference displacement points before evaluation
        :type displacement_generator: callable | None
        :param coordinate_prep: a callable to preprocess the flattened center coordinates before finite-differencing
        :type coordinate_prep: callable | None
        :param index_filter: a callable/filter restricting which derivative-tensor index combinations are computed
        :type index_filter: callable | None
        :param analytic_derivative_order: the derivative order already available analytically
        :type analytic_derivative_order: int | None
        :param opts: extra options, split between `FiniteDifferenceDerivative` construction options and options forwarded to `evaluate_term`
        :type opts: dict
        :return: the finite-difference-computed internal-coordinate derivative tensors
        :rtype: list[np.ndarray]
        """
        opts = dev.OptionsSet(opts)
        fd_opts = opts.filter(FiniteDifferenceDerivative)
        opts = opts.exclude(FiniteDifferenceDerivative)

        if batched_orders is None:
            batched_orders = self.batched_orders
        if analytic_derivative_order is None:
            analytic_derivative_order = self.analytic_derivative_order

        reembed = dev.str_is(self.use_internals, 'reembed')
        if (
                reembed
                and coords.ndim > 1
                and coords.shape[-2:] == (len(self.embedding.masses), 3)
        ):
            internals = self.embedding.get_internals(coords=coords, strip_embedding=self.strip_embedding)
        else:
            internals = coords

        base_shape = internals.shape[:-1]
        coord_shape = internals.shape[-1:]
        flat_coords = internals.reshape(-1, np.prod(coord_shape, dtype=int))

        def derivs(structs, center=flat_coords):
            """
            **LLM Docstring**

            Evaluate `evaluate_term` at a batch of displaced internal-coordinate structures, converting each displaced point back to Cartesians (and re-embedding any lower-order analytic derivatives via `embed_derivs`) when `use_internals` is `'reembed'`, for use as the base function inside the enclosing finite-difference derivative computation.

            :param structs: the (possibly displacement-generator-transformed) flattened internal-coordinate points to evaluate at
            :type structs: np.ndarray
            :param center: the reference (un-displaced) flattened internal coordinates, from the enclosing scope
            :type center: np.ndarray
            :return: the evaluated (and, if reembedding, Cartesian-re-expressed) term(s)
            :rtype: np.ndarray
            """
            if displacement_generator is not None:
                structs = displacement_generator(structs, evaluator=self, center=center)
            reconst = structs.reshape(structs.shape[:-1] + coord_shape)
            if reembed:
                reconst = self.unembed_coords(reconst)
                if not batched_orders:
                    res = [
                        self.evaluate_term(reconst, o, **opts)
                        for o in range(analytic_derivative_order + 1)
                    ]
                else:
                    res = self.evaluate_term(reconst, analytic_derivative_order, **opts)
                res = res[:1] + self.embed_derivs(reconst, res[1:])
                return res[-1]
            else:
                ders = self.evaluate_term(reconst, analytic_derivative_order, **opts)
                if batched_orders:
                    return ders[-1]
                else:
                    return ders

        der = FiniteDifferenceDerivative(derivs,
                                         function_shape=((0,), (0,) * analytic_derivative_order),
                                         **self.get_fd_opts(**fd_opts)
                                         )

        if coordinate_prep is not None:
            flat_coords  = coordinate_prep(flat_coords, evaluator=self)

        tensors = der.derivatives(flat_coords).derivative_tensor(
            list(range(1, (order - analytic_derivative_order) + 1)),
            pos_filter=index_filter
        )
        if flat_coords.shape[0] > 1:
            tensors = [np.moveaxis(d, n+1, 0) for n,d in enumerate(tensors)]
        res = [
            (
                t.reshape(base_shape + t.shape[1:])
                    if flat_coords.shape[0] != 1 else
                t.reshape(base_shape + t.shape)
            )
                if t.ndim > 0 else
            t
            for t in tensors
        ]

        return res

    supports_internals = False
    property_units = None
    target_property_units = None
    distance_units = 'Angstroms'
    batched_orders = False
    analytic_derivative_order = 1
    multi_expansion_order = 0
    def evaluate_expansion(self,
                           coords,
                           order=0,
                           analytic_derivative_order=None,
                           batched_orders=None,
                           logger=None,
                           fd_displacement_generator=None,
                           fd_coordinate_prep=None,
                           fd_index_filter=None,
                           fd_handler=None,
                           **opts):
        """
        **LLM Docstring**

        Build the property's Taylor expansion up to the requested derivative order(s), combining analytic derivatives (via `evaluate_term`, for orders up to `analytic_derivative_order`) with finite-difference derivatives (via `fd_handler`, for any higher orders requested), then applying the evaluator's energy/distance unit conversions.

        :param coords: the coordinates to evaluate the expansion at
        :type coords: np.ndarray
        :param order: the derivative order(s) to compute, either a single integer (meaning `0..order`) or an explicit sorted list of orders
        :type order: int | list[int]
        :param analytic_derivative_order: the highest order to compute analytically via `evaluate_term`; defaults to `self.analytic_derivative_order`
        :type analytic_derivative_order: int | None
        :param batched_orders: whether `evaluate_term` computes all analytic orders in one batched call
        :type batched_orders: bool | None
        :param logger: accepted for interface consistency but not used in this method's body
        :type logger: Logger | None
        :param fd_displacement_generator: forwarded to `fd_handler` as `displacement_generator`
        :type fd_displacement_generator: callable | None
        :param fd_coordinate_prep: forwarded to `fd_handler` as `coordinate_prep`
        :type fd_coordinate_prep: callable | None
        :param fd_index_filter: forwarded to `fd_handler` as `index_filter`
        :type fd_index_filter: callable | None
        :param fd_handler: the finite-difference routine to use for orders beyond `analytic_derivative_order`; defaults to `self.finite_difference_derivs`
        :type fd_handler: callable | None
        :param opts: extra options, merged with `self.defaults` and split between finite-difference construction options and options forwarded to `evaluate_term`
        :type opts: dict
        :return: the expansion terms, one per requested order (in units converted per `property_units`/`distance_units`)
        :rtype: list
        """

        opts = dev.OptionsSet(dict(self.defaults, **opts))
        fd_opts = opts.filter(FiniteDifferenceDerivative)
        opts = opts.exclude(FiniteDifferenceDerivative)

        if nput.is_numeric(order):
            order = list(range(order+1))
        # TODO: we assume order sorted, handle when it's not
        expansion = []
        ad = opts.pop('analytic_derivative_order', None)
        if analytic_derivative_order is None:
            analytic_derivative_order = ad
        if analytic_derivative_order is None:
            analytic_derivative_order = self.analytic_derivative_order
        if batched_orders is None:
            batched_orders = self.batched_orders
        if batched_orders:
            terms = self.evaluate_term(coords, min(analytic_derivative_order, order[-1]), **opts)
            if self.multi_expansion_order > 0:
                expansion.extend(
                    [t[o] for o in order if o <= analytic_derivative_order]
                    for t in terms
                )
            else:
                expansion.extend([terms[o] for o in order if o <= analytic_derivative_order])
        else:
            for i in order:
                if i > analytic_derivative_order: break
                expansion.append(self.evaluate_term(coords, i, **opts))
        if order[-1] > analytic_derivative_order:
            if fd_handler is None: fd_handler = self.finite_difference_derivs

            # self.fd_derivs(**opts)
            #TODO: handle disjoin list of orders...
            expansion.extend(
                fd_handler(coords, order[-1], batched_orders=batched_orders,
                           displacement_generator=fd_displacement_generator,
                           coordinate_prep=fd_coordinate_prep,
                           index_filter=fd_index_filter,
                           analytic_derivative_order=analytic_derivative_order,
                           **fd_opts)
            )

        if self.property_units is None:
            eng_conv = 1
        else:
            eng_conv = UnitsData.convert(self.property_units, self.target_property_units)

        if self.distance_units is None:
            bohr_conv = 1
        else:
            bohr_conv = UnitsData.convert(self.distance_units, "BohrRadius")

        if self.multi_expansion_order > 0:
            expansion = [
                [
                    e * eng_conv / (bohr_conv ** o)
                    for o, e in zip(order, subexpansion)
                ] for subexpansion in expansion
            ]
        else:
            expansion = [
                e * eng_conv / (bohr_conv ** o)
                for o, e in zip(order, expansion)
            ]
        return expansion

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 fd_handler=None,
                 analytic_derivative_order=None,
                 unembed_derivatives=True,
                 **opts):
        """
        **LLM Docstring**

        Top-level property-evaluation entry point: embeds the input coordinates into the evaluator's working representation, computes the Taylor expansion (temporarily disabling internal double-handling/permutation/re-embedding so `evaluate_expansion` operates on already-embedded coordinates), re-expands analytic-order derivatives back through the Cartesians-by-internals Jacobian if working in internal coordinates, and finally un-embeds the resulting derivatives back to plain Cartesian coordinates (unless `unembed_derivatives` is `False`).

        :param coords: the (Cartesian) coordinates to evaluate the property at
        :type coords: np.ndarray
        :param order: the derivative order(s) to compute
        :type order: int | list[int]
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param fd_handler: an explicit finite-difference handler to use; defaults to `internal_finite_difference_derivs` if working in internal coordinates
        :type fd_handler: callable | None
        :param analytic_derivative_order: the highest order to compute analytically
        :type analytic_derivative_order: int | None
        :param unembed_derivatives: whether to convert the resulting derivatives back to plain Cartesian coordinates
        :type unembed_derivatives: bool
        :param opts: extra options forwarded to `evaluate_expansion`
        :type opts: dict
        :return: the property expansion (value plus derivatives), in Cartesian coordinates unless `unembed_derivatives` is `False`
        :rtype: list
        :raises NotImplementedError: if `multi_expansion_order > 0` and internal-coordinate re-expansion is needed
        """
        use_internals = self.use_internals
        reembed = self.reembed_cartesians
        perm = self.permutation
        base_coords = coords
        coords = self.embed_coords(coords, embed_reembedded=False)
        in_internals = self.use_internal_coordinate_handlers()

        if fd_handler is None and self.use_internal_coordinate_handlers():
            fd_handler = self.internal_finite_difference_derivs

        if analytic_derivative_order is None:
            analytic_derivative_order = self.defaults.get('analytic_derivative_order')
        if analytic_derivative_order is None:
            analytic_derivative_order = self.analytic_derivative_order

        try:
            if not dev.str_is(self.use_internals, 'reembed'):
                self.use_internals = False
            self.reembed_cartesians = False
            self.permutation = None

            # TODO: make this less hacky...
            expansion = self.evaluate_expansion(
                coords,
                order=order,
                logger=logger,
                fd_handler=fd_handler,
                analytic_derivative_order=analytic_derivative_order,
                **opts
            )
            if in_internals and analytic_derivative_order > 0:
                tf = self.embedding.get_cartesians_by_internals(
                    coords=base_coords,
                    order=analytic_derivative_order,
                    strip_embedding=self.strip_embedding
                )

                if self.multi_expansion_order > 0:
                    raise NotImplementedError("transformations of multi-expansions not implemented yet")
                else:
                    expansion = expansion[:1] + nput.tensor_reexpand(
                        tf,
                        expansion[1:analytic_derivative_order+1]
                    ) + expansion[analytic_derivative_order+1:]

        finally:
            self.use_internals = use_internals
            self.reembed_cartesians = reembed
            self.permutation = perm

        if unembed_derivatives:
            if self.multi_expansion_order > 0:
                expansion = [
                    [e[0]] + self.unembed_derivs(base_coords, coords, e[1:], embed_reembedded=True)
                    for e in expansion
                ]
            else:
                expansion = [expansion[0]] + self.unembed_derivs(base_coords, coords, expansion[1:], embed_reembedded=True)
        return expansion

    @staticmethod
    def _diag_inds(pos):
        """
        **LLM Docstring**

        Filter a list of index tuples down to just the "diagonal" ones, where every index in the tuple is the same.

        :param pos: the index tuples to filter
        :type pos: list[tuple]
        :return: the diagonal index tuples
        :rtype: list[tuple]
        """
        return [p for p in pos if len(np.unique(p)) == 1]
    def partial_force_field(self,
                            coords:np.ndarray,
                            modes,
                            mode_distance_units='BohrRadius',
                            analytic_derivative_order=None,
                            order=4,
                            index_filter=None,
                            fd_handler=None,
                            mesh_spacing=0.5, # nms require much larger steps
                            **opts
                            ):
        """
        **LLM Docstring**

        Compute a partial force-field-style expansion of the property restricted to displacements along a given set of normal modes (rather than full Cartesian/internal derivatives), by generating finite-difference displacement points along the mode-coordinate directions and only requesting the "diagonal" derivative terms by default.

        :param coords: the reference coordinates to expand around
        :type coords: np.ndarray
        :param modes: the normal modes to expand the property in
        :type modes: object
        :param mode_distance_units: the units the mode displacement coordinates are given in, converted to `self.distance_units`
        :type mode_distance_units: str
        :param analytic_derivative_order: the highest order to compute analytically
        :type analytic_derivative_order: int | None
        :param order: the highest derivative order to compute
        :type order: int
        :param index_filter: which derivative-tensor index combinations to compute; defaults to `_diag_inds` (only same-mode-repeated terms)
        :type index_filter: callable | None
        :param fd_handler: an explicit finite-difference handler to use
        :type fd_handler: callable | None
        :param mesh_spacing: the finite-difference step size (normal modes typically need larger steps than Cartesian/internal coordinates)
        :type mesh_spacing: float
        :param opts: extra options forwarded to `evaluate`
        :type opts: dict
        :return: the mode-basis partial force-field expansion, with derivative orders beyond `analytic_derivative_order` rescaled back from mode units
        :rtype: list
        """
        conv = UnitsData.convert(mode_distance_units, self.distance_units)
        if index_filter is None:
            index_filter = self._diag_inds

        if analytic_derivative_order is None:
            analytic_derivative_order = self.analytic_derivative_order

        def get_displacements(displacements, *, evaluator, center):
            """
            **LLM Docstring**

            Convert normal-mode-basis displacement values into Cartesian (or working-representation) coordinate displacements, added onto the enclosing `center` reference coordinates, using the mode-to-coordinate transformation matrix from the enclosing scope.

            :param displacements: the mode-basis displacement values
            :type displacements: np.ndarray
            :param evaluator: the evaluator invoking this displacement generator (unused directly in the body)
            :type evaluator: PropertyEvaluator
            :param center: the reference coordinates to displace from
            :type center: np.ndarray
            :return: the displaced coordinates
            :rtype: np.ndarray
            """
            disps = conv * np.tensordot(displacements, modes.coords_by_modes, axes=[-1, 0])
            dd = disps.ndim - center.ndim
            if dd > 0:
                center = np.expand_dims(center, list(range(dd)))
            return center + disps

        def prep_coordinates(flat_coords, *, evaluator):
            """
            **LLM Docstring**

            Build a zeroed placeholder array shaped to hold mode-basis coordinates (matching the batch shape of `flat_coords` but with the mode dimension from the enclosing `modes`), used as the finite-difference center point for a mode-basis expansion.

            :param flat_coords: the flattened Cartesian coordinates whose batch shape should be matched
            :type flat_coords: np.ndarray
            :param evaluator: the evaluator invoking this preparation function (unused directly in the body)
            :type evaluator: PropertyEvaluator
            :return: the zeroed mode-basis coordinate array
            :rtype: np.ndarray
            """
            return np.zeros(
                flat_coords.shape[:-1] +(modes.coords_by_modes.shape[0],),
                dtype=float
            )

        base_expansion = self.evaluate(coords,
                                       order=order,
                                       fd_displacement_generator=get_displacements,
                                       fd_coordinate_prep=prep_coordinates,
                                       fd_index_filter=index_filter,
                                       fd_handler=fd_handler,
                                       analytic_derivative_order=analytic_derivative_order,
                                       mesh_spacing=mesh_spacing,
                                       **opts
                                       )

        expansion = base_expansion[:(analytic_derivative_order+1)] + [
            e / (conv**(n+1))
            for n,e in enumerate(base_expansion[analytic_derivative_order+1:])
        ]

        return expansion


    evaluator_registry: dict
    @classmethod
    def register(cls, name, method=None):
        """
        **LLM Docstring**

        Register a named evaluator-construction method (or function) in the class's `evaluator_registry`, usable directly or as a decorator.

        :param name: the registry key to register under, or (if `method` is `None` and `name` has a `.name` attribute) the method/callable itself
        :type name: str | callable
        :param method: the evaluator-construction method/callable to register; if `None`, `register` returns a decorator instead
        :type method: callable | None
        :return: the registered method (if called directly), or a decorator function (if used as `@register(name)`)
        :rtype: callable
        """
        if method is None and hasattr(name, 'name'):
            method = name
            name = method.name
        if method is not None:
            cls.evaluator_registry[name] = method
            return method
        else:
            def register(method, name=name):
                """
                **LLM Docstring**

                Decorator form of `register`, closing over `name` from the enclosing call, used when `register` is invoked without an explicit `method`.

                :param method: the evaluator-construction method/callable to register
                :type method: callable
                :param name: the registry key to register under, from the enclosing scope
                :type name: str
                :return: the result of `cls.register(name, method)`
                :rtype: callable
                """
                return cls.register(name, method)
            return register

    # evaluator_types = {}
    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The set of built-in named evaluator constructors for this evaluator type. Returns an empty dict on the base class; concrete subclasses override this to expose their built-in evaluator options.

        :return: mapping from evaluator name to constructor
        :rtype: dict
        """
        return {}

    @classmethod
    def get_evaluator_map(cls):
        """
        **LLM Docstring**

        The full set of named evaluator constructors available for this evaluator type: the built-ins from `get_evaluators`, overlaid with any user-registered ones from `evaluator_registry`.

        :return: the combined evaluator-name-to-constructor mapping
        :rtype: dict
        """
        return cls.evaluator_registry | cls.get_evaluators()

    @classmethod
    def get_evaluators_by_attributes(cls):
        """
        **LLM Docstring**

        Optional mapping used to resolve an evaluator by a set of attribute-based hints rather than by exact name, for use by the dispatch machinery. Returns `None` on the base class (no attribute-based dispatch); concrete subclasses may override this.

        :return: the attribute-based dispatch mapping, or `None`
        :rtype: dict | None
        """
        return None

    _profile_dispatch = dev.uninitialized
    default_evaluator_type = 'expansion'
    @classmethod
    def profile_generator_dispatch(cls):
        """
        **LLM Docstring**

        Lazily build (and cache) the `OptionsMethodDispatch` object used to resolve a named/attribute-described evaluator specification into a concrete evaluator constructor and options.

        :return: the cached dispatch object
        :rtype: dev.OptionsMethodDispatch
        """
        cls._profile_dispatch = dev.handle_uninitialized(
            cls._profile_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_evaluator_map,),
            kwargs=dict(
                # default_method=cls.default_evaluator_type,
                attributes_map=cls.get_evaluators_by_attributes(),
            )
        )
        return cls._profile_dispatch

    @classmethod
    @abc.abstractmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        Abstract hook returning the `PropertyFunctionEvaluator`-family class to use when a bare callable is supplied as an evaluator specification. Concrete evaluator subclasses must implement this.

        :return: the default function-evaluator class
        :rtype: type
        """
        ...

    @classmethod
    def _resolve_torch_device(cls, device=None):
        """
        **LLM Docstring**

        Resolve which PyTorch device to run a model on, auto-detecting CUDA availability if not explicitly given.

        :param device: an explicit device string; if `None`, resolves to `'cuda'` if available, otherwise `'cpu'`
        :type device: str | None
        :return: the resolved device string
        :rtype: str
        """
        if device is None:
            import torch
            device = (
                # 'mps'
                #     if torch.backends.mps.is_available() else
                'cuda'
                    if torch.cuda.is_available() else
                'cpu'
            )
        return device
    @classmethod
    def handle_specialization(cls, tag):
        """
        **LLM Docstring**

        Hook for resolving an evaluator "specialization" tag (the part after a `:` in an evaluator name string, e.g. `'method:tag'`) into extra options. Not supported on the base class; concrete evaluator subclasses may override this to support tags.

        :param tag: the specialization tag to handle
        :type tag: str
        :return: never returns on the base class
        :rtype: dict
        :raises NotImplementedError: always, on classes that don't override this
        """
        raise NotImplementedError(f"{cls.__name__} does not handle specializations (got {tag})")

    @classmethod
    def resolve_evaluator(cls, name):
        """
        **LLM Docstring**

        Resolve an evaluator specification -- a bare callable, an object already exposing `evaluate`, or a (possibly `'name:tag'`-qualified) name/dict/tuple understood by the evaluator-dispatch machinery -- into a concrete evaluator class/instance plus its construction options.

        :param name: the evaluator specification to resolve
        :type name: object | str | dict | tuple | list | None
        :return: `(evaluator, opts)` -- the resolved evaluator class/instance and its construction options
        :rtype: tuple
        :raises ValueError: if `name` is a non-callable, non-string/dict/tuple/list object that also doesn't expose `evaluate`
        """
        if (
                name is not None
                and not isinstance(name, (str, dict, tuple, list))
        ):
            if hasattr(name, 'evaluate'):
                return name, {}
            elif callable(name):
                eval = cls.get_default_function_evaluator_type()
                eval.bind_default(name)
                return eval, {}
            else:
                raise ValueError(f"can't get construct evalautor of type {cls} from {name}")
        elif isinstance(name, str) and ':' in name:
            name, tag = name.split(':', 1)
        else:
            tag = None
        if dev.is_dict_like(name) and hasattr(name, 'method'):
            if ':' in name['method']:
                name = name.copy()
                n, tag = name['method'].split(':', 1)
                name['method'] = n

        base, opts = cls.profile_generator_dispatch().resolve(name)
        if tag is not None:
            opts = opts | base.handle_specialization(tag)
        return base, opts

    class quiet_mode:
        def __init__(self, quiet=True):
            """
            **LLM Docstring**

            Set up a context manager that optionally silences stdout for the duration of its `with` block.

            :param quiet: whether to actually silence stdout; if `False`, the context manager is a no-op
            :type quiet: bool
            :return: None
            :rtype: None
            """
            self.quiet = quiet
            self._stdout = None
            self._devnull = None
        def __enter__(self):
            """
            **LLM Docstring**

            Redirect `sys.stdout` to `os.devnull` if `self.quiet` is set.

            :return: None
            :rtype: None
            """
            if self.quiet:
                self._stdout = sys.stdout
                self._devnull = open(os.devnull, 'w+').__enter__()
                sys.stdout = self._devnull
        def __exit__(self, exc_type, exc_val, exc_tb):
            """
            **LLM Docstring**

            Restore `sys.stdout` (and close the devnull stream) if `self.quiet` is set.

            :param exc_type: the exception type, if an exception occurred in the `with` block
            :type exc_type: type | None
            :param exc_val: the exception instance, if any
            :type exc_val: BaseException | None
            :param exc_tb: the exception traceback, if any
            :type exc_tb: object | None
            :return: None
            :rtype: None
            """
            if self.quiet:
                self._devnull.__exit__(exc_type, exc_val, exc_tb)
                sys.stdout = self._stdout

class PropertyFunctionEvaluator(PropertyEvaluator):
    # Slightly strange use case...we mostly use this as a component to
    # implement individual methods

    def __init__(self,
                 potential_function,
                 property_units=None,
                 distance_units='Angstroms',
                 batched_orders=False,
                 analytic_derivative_order=np.inf,
                 **defaults
                 ):
        """
        **LLM Docstring**

        Wrap a bare Python function as a `PropertyEvaluator`, storing its unit conventions and derivative-computation characteristics.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)`
        :type potential_function: callable
        :param property_units: the energy/property units the function's output is given in
        :type property_units: str | None
        :param distance_units: the distance units the function expects its input coordinates in
        :type distance_units: str
        :param batched_orders: whether the function returns all requested derivative orders in a single batched call
        :type batched_orders: bool
        :param analytic_derivative_order: the highest order the function can compute analytically; defaults to `inf` (the function is assumed to support any order directly)
        :type analytic_derivative_order: float | int
        :param defaults: extra default options forwarded to `PropertyEvaluator.__init__` (note: dispatched via `super(type(self).__bases__[0], self)`, skipping the immediate `PropertyEvaluator.__init__` override chain)
        :type defaults: dict
        :return: None
        :rtype: None
        """
        super(type(self).__bases__[0], self).__init__(**defaults)
        self.property_function = potential_function
        self.batched_orders = batched_orders
        self.analytic_derivative_order = analytic_derivative_order
        self.property_units = property_units
        self.distance_units = distance_units

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped function directly at the given coordinates/order.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        return self.property_function(coords, order=order, **opts)

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default property function to be used the next time `initialize_from_mol` is called without an explicit `property_function`.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        cls.default_property_function = potential
    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        """
        **LLM Docstring**

        Resolve the property function to wrap. On this base class, simply returns `prop_func` unchanged; subclasses may override this to do additional resolution (e.g. binding molecule-specific parameters).

        :param prop_func: the candidate property function
        :type prop_func: callable
        :param mol: the molecule the evaluator is being built for
        :type mol: AbstractMolecule
        :param opts: extra options, unused on the base class
        :type opts: dict
        :return: the resolved property function
        :rtype: callable
        """
        return prop_func
    @classmethod
    def initializer(cls, mol, **opts):
        """
        **LLM Docstring**

        Build a bare (function-less) instance of this evaluator type; used as a lightweight constructor hook by subclasses/dispatch code that will bind a function afterward.

        :param mol: the molecule the evaluator is being built for (unused directly in this method's body)
        :type mol: AbstractMolecule
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed evaluator, with `potential_function=None`
        :rtype: PropertyFunctionEvaluator
        """
        return cls(None, **opts)
    @classmethod
    def initialize_from_mol(root,
                            cls,
                            mol,
                            property_function=None,
                            *,
                            embedding=None,
                            use_internals=None,
                            reembed_cartesians=None,
                            property_units=None,
                            distance_units=None,
                            batched_orders=None,
                            analytic_derivative_order=None,
                            **opts):
        """
        **LLM Docstring**

        Build a `PropertyFunctionEvaluator` for a molecule, resolving the property function to wrap (falling back to a previously `bind_default`-registered default function if none is given) and filling in construction options (embedding, unit conventions, derivative-order handling) from the molecule wherever not explicitly overridden.

        :param root: the class this method is unbound from (Python descriptor mechanics artifact; effectively unused directly)
        :type root: type
        :param cls: the evaluator class to instantiate
        :type cls: type
        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the function to wrap; falls back to `cls.default_property_function` if not given
        :type property_function: callable | None
        :param embedding: an explicit embedding to use instead of `mol.embedding`
        :type embedding: MolecularEmbedding | None
        :param use_internals: whether to evaluate in internal coordinates
        :type use_internals: bool | None
        :param reembed_cartesians: whether to Eckart-reembed Cartesian coordinates
        :type reembed_cartesians: bool | None
        :param property_units: the units the function's output is given in
        :type property_units: str | None
        :param distance_units: the units the function expects its input in
        :type distance_units: str | None
        :param batched_orders: whether the function returns all derivative orders in one batched call
        :type batched_orders: bool | None
        :param analytic_derivative_order: the highest order the function computes analytically
        :type analytic_derivative_order: int | None
        :param opts: extra options forwarded to `get_property_function` and the constructor
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PropertyFunctionEvaluator
        :raises ValueError: if no property function could be resolved
        """
        if property_function is None:
            property_function = cls.default_property_function
            cls.default_property_function = None
        property_function, opts = cls.get_property_function(property_function, mol, **opts)
        if property_function is None:
            raise ValueError(f"can't evaluate with property function {property_function}")
        init_opts = {
            k: v for k, v in dict(
                embedding=mol.embedding if embedding is None else embedding,
                use_internals=use_internals,
                reembed_cartesians=reembed_cartesians,
                property_units=property_units,
                distance_units=distance_units,
                batched_orders=batched_orders,
                analytic_derivative_order=analytic_derivative_order,
            ).items()
            if v is not None
        }
        return cls(
            property_function,
            **opts,
            **init_opts
        )

class EnergyEvaluator(PropertyEvaluator):

    target_property_units = 'Hartrees'
    evaluator_registry = {}

    # @classmethod
    # def get_evaluators(cls):
    #     return {
    #         'rdkit': RDKitEnergyEvaluator,
    #         'aimnet2': AIMNet2EnergyEvaluator,
    #         'mace': MACEEnergyEvaluator,
    #         'uma': UMAEnergyEvaluator,
    #         'hipnn':HIPNNEnergyEvaluator,
    #         'ase': ASECalcEnergyEvaluator,
    #         'xtb': XTBEnergyEvaluator,
    #         'pyscf': PySCFEnergyEvaluator,
    #         'expansion': PotentialExpansionEnergyEvaluator
    #     }
    @classmethod
    def get_evaluators_by_attributes(cls):
        """
        **LLM Docstring**

        Attribute-based dispatch mapping for energy evaluators: recognizes a `potential_function` keyword as identifying a `PotentialFunctionEnergyEvaluator`.

        :return: the attribute-tuple-to-evaluator-class mapping
        :rtype: dict
        """
        return {
            ('potential_function',):PotentialFunctionEnergyEvaluator
        }
    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as an energy-evaluator specification.

        :return: `PotentialFunctionEnergyEvaluator`
        :rtype: type
        """
        return PotentialFunctionEnergyEvaluator

    def _modify_ase_calc(self, #TODO: add in orthogonal projections as well
                         calc,
                         gradient_modification_function=None,
                         gradient_modification_mode=None,
                         orthogonal_projection_generator=None,
                         **opts
                         ):
        """
        **LLM Docstring**

        Wrap an ASE calculator's `calculate` method so that, if a gradient-modification function is supplied, the computed forces are intercepted, converted to a flat gradient, passed through `_modify_gradient`, and written back into the calculator's results in ASE's expected shape.

        :param calc: the ASE calculator to modify in place
        :type calc: object
        :param gradient_modification_function: the function used to modify the raw gradient/forces; if `None`, `calc` is returned unmodified
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode forwarded to `_modify_gradient` (e.g. `'shift'`)
        :type gradient_modification_mode: str | None
        :param orthogonal_projection_generator: an optional per-geometry projector generator forwarded to `_modify_gradient`
        :type orthogonal_projection_generator: callable | None
        :param opts: extra options forwarded to `_modify_gradient`
        :type opts: dict
        :return: the (possibly modified) calculator
        :rtype: object
        """
        if gradient_modification_function is not None:
            calc._old_calculate = calc.calculate
            def prep_coords(atoms=None, *etc, **kwetc):
                """
                **LLM Docstring**

                Extract the Cartesian atom positions from an ASE `Atoms` object (or `calc.atoms` if none given), for use as the coordinate argument passed to the gradient-modification function.

                :param atoms: the ASE atoms object to read positions from; defaults to `calc.atoms`
                :type atoms: object | None
                :param etc: extra positional arguments, unused
                :type etc: tuple
                :param kwetc: extra keyword arguments, unused
                :type kwetc: dict
                :return: the atom positions
                :rtype: np.ndarray
                """
                if atoms is None: atoms = calc.atoms
                return atoms.positions
            def prep_grad(_):
                """
                **LLM Docstring**

                Extract (and cache) the computed forces from an ASE calculator's results dict, flattening them for the gradient-modification function; returns `None` (short-circuiting the modification) if forces aren't present or have already been modified once.

                :param _: the calculator's results dict (positional, unnamed in the source)
                :type _: dict
                :return: `(flat_forces, results)`, or `(None, results)` if there's nothing new to modify
                :rtype: tuple
                """
                if 'forces' not in calc.results or 'old_forces' in calc.results:
                    return None, calc.results
                else:
                    f = calc.results.pop('forces')
                    calc.results['_cached_force'] = f
                    return f.reshape(f.shape[:-2] + (-1,)), calc.results
            def post_grad(grad, res):
                """
                **LLM Docstring**

                Write a modified gradient back into an ASE calculator's results dict in the expected `(natoms, 3)`-shaped `'forces'` entry, preserving the original (pre-modification) forces under `'old_forces'`.

                :param grad: the modified flat gradient
                :type grad: np.ndarray
                :param res: the results dict to update, from the enclosing `prep_grad` call
                :type res: dict
                :return: the updated results dict
                :rtype: dict
                """
                # from .Molecule import Molecule
                # from McUtils.ExternalPrograms import ASEMolecule
                #
                # Molecule.from_ase(ASEMolecule.from_atoms(calc.atoms)).animate_coordinate(
                #     0,
                #     coordinate_expansion=[nput.vec_normalize(grad.reshape((1, -1)))],
                #     backend='x3d'
                # ).show()

                f0 = res.get('_cached_force')
                if f0 is not None:
                    res['old_forces'] = f0
                res['forces'] = grad.reshape(grad.shape[:-1] + (-1, 3))
                return res
            calc.calculate = self._modify_gradient(
                calc.calculate,
                gradient_modification_function,
                modification_mode=gradient_modification_mode,
                coord_prep=prep_coords,
                grad_prep=prep_grad,
                grad_post=post_grad,
                use_forces=True,
                orthogonal_projection_generator=orthogonal_projection_generator,
                **opts
            )
        return calc
    def to_ase(self,
               gradient_modification_function=None,
               gradient_modification_mode=None,
               convert_modification_distances=None,
               convert_modification_energies=None,
               orthogonal_projection_generator=None,
               **etc
               ):
        """
        **LLM Docstring**

        Wrap this evaluator as an ASE `Calculator` object (via `ASECalculator`), optionally intercepting/modifying its gradient via `_modify_ase_calc`.

        :param gradient_modification_function: an optional function to modify the computed gradient/forces
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode to use
        :type gradient_modification_mode: str | None
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool | None
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool | None
        :param orthogonal_projection_generator: an optional per-geometry projector generator applied to the modified gradient
        :type orthogonal_projection_generator: callable | None
        :param etc: extra options forwarded to `ASECalculator`
        :type etc: dict
        :return: the constructed (and possibly gradient-modified) ASE calculator
        :rtype: object
        """
        from McUtils.ExternalPrograms import ASECalculator

        cacl = ASECalculator(self.evaluate_term, **etc)
        return self._modify_ase_calc(cacl,
                                     gradient_modification_function=gradient_modification_function,
                                     modification_mode=gradient_modification_mode,
                                     convert_modification_distances=convert_modification_distances,
                                     convert_modification_energies=convert_modification_energies,
                                     orthogonal_projection_generator=orthogonal_projection_generator)

    def _modify_pysis_calc(self,
                           calc,
                           gradient_modification_function=None,
                           gradient_modification_mode=None,
                           convert_modification_distances=None,
                           convert_modification_energies=None,
                           **opts
                           ):
        """
        **LLM Docstring**

        Wrap a pysisyphus calculator's `get_forces` (and `get_hessian`, if present) so that, if a gradient-modification function is supplied, the computed forces/Hessian are intercepted, passed through `_modify_gradient`, and written back.

        :param calc: the pysisyphus calculator to modify in place
        :type calc: object
        :param gradient_modification_function: the function used to modify the raw forces; if `None`, `calc` is returned unmodified
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode forwarded to `_modify_gradient`
        :type gradient_modification_mode: str | None
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool | None
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool | None
        :param opts: extra options forwarded to `_modify_gradient`
        :type opts: dict
        :return: the (possibly modified) calculator
        :rtype: object
        """
        if gradient_modification_function is not None:
            calc._old_get_forces = calc.get_forces
            if hasattr(calc, 'get_hessian'):
                calc._old_get_hessian = calc.get_hessian

            f_cache = [None]
            def prep_coords(atoms, coords):
                """
                **LLM Docstring**

                Pass through the coordinates unchanged, matching the `(atoms, coords) -> coords` signature expected by the pysisyphus gradient wrapper.

                :param atoms: unused
                :type atoms: object
                :param coords: the coordinates to pass through
                :type coords: np.ndarray
                :return: `coords`, unchanged
                :rtype: np.ndarray
                """
                return coords
            def prep_grad(res):
                """
                **LLM Docstring**

                Extract the `'forces'` entry from a pysisyphus results dict for use by the gradient-modification function.

                :param res: the results dict
                :type res: dict
                :return: `(forces, res)`
                :rtype: tuple
                """
                f = res.pop('forces')
                return f, res
            def post_grad(grad, res):
                """
                **LLM Docstring**

                Write a modified gradient back into a pysisyphus results dict's `'forces'` entry.

                :param grad: the modified gradient
                :type grad: np.ndarray
                :param res: the results dict to update
                :type res: dict
                :return: the updated results dict
                :rtype: dict
                """
                res['forces'] = grad
                return res

            calc.get_forces = self._modify_gradient(
                calc.get_forces,
                gradient_modification_function,
                coord_prep=prep_coords,
                use_forces=True,
                modification_mode=gradient_modification_mode,
                convert_modification_distances=convert_modification_distances,
                convert_modification_energies=convert_modification_energies,
                grad_prep=prep_grad,
                grad_post=post_grad,
                **opts
            )
            if hasattr(calc, 'get_hessian'):
                calc.get_hessian = self._modify_gradient(
                    calc.get_hessian,
                    gradient_modification_function,
                    coord_prep=prep_coords,
                    use_forces=True,
                    modification_mode=gradient_modification_mode,
                    convert_modification_distances=convert_modification_distances,
                    convert_modification_energies=convert_modification_energies,
                    grad_prep=prep_grad,
                    grad_post=post_grad,
                    **opts
                )
        return calc
    def to_pysis(self,
                 gradient_modification_function=None,
                 gradient_modification_mode=None,
                 convert_modification_distances=None,
                 convert_modification_energies=None,
                 orthogonal_projection_generator=None,
                 **etc):
        """
        **LLM Docstring**

        Wrap this evaluator as a pysisyphus calculator (via `PysisCalculator`), optionally intercepting/modifying its gradient via `_modify_pysis_calc`.

        :param gradient_modification_function: an optional function to modify the computed forces
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode to use
        :type gradient_modification_mode: str | None
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool | None
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool | None
        :param orthogonal_projection_generator: forwarded to `_modify_pysis_calc`
        :type orthogonal_projection_generator: callable | None
        :param etc: extra options forwarded to `PysisCalculator`
        :type etc: dict
        :return: the constructed (and possibly gradient-modified) pysisyphus calculator
        :rtype: object
        """
        from McUtils.ExternalPrograms import PysisCalculator

        calc = PysisCalculator(self.evaluate_term,
                               distance_units=self.distance_units,
                               energy_units=self.target_property_units,
                               batched_orders=self.batched_orders,
                               **etc)
        calc = self._modify_pysis_calc(calc,
                                       use_forces=True,
                                       gradient_modification_function=gradient_modification_function,
                                       gradient_modification_mode=gradient_modification_mode,
                                       convert_modification_distances=convert_modification_distances,
                                       convert_modification_energies=convert_modification_energies,
                                       orthogonal_projection_generator=orthogonal_projection_generator)

        return calc

    def minimizer_function_by_order(self, order, allow_fd=False, modifier=None, **opts):
        """
        **LLM Docstring**

        Build a callable of the `scipy.optimize`-style `(coords, *extra) -> value` form for a specific derivative order (energy, gradient, or Hessian), using analytic evaluation when available and (optionally) finite differences otherwise, converting the result to Hartrees.

        :param order: the derivative order the function should evaluate (`0` energy, `1` gradient, `2` Hessian)
        :type order: int
        :param allow_fd: whether to fall back to `finite_difference_derivs` when the analytic derivative order isn't high enough
        :type allow_fd: bool
        :param modifier: an optional function to post-process the returned value
        :type modifier: callable | None
        :param opts: extra options forwarded to `evaluate_term`/`finite_difference_derivs`
        :type opts: dict
        :return: the constructed minimizer-compatible function, or `None` if the requested order isn't available and `allow_fd` is `False`
        :rtype: callable | None
        """
        conv = UnitsData.convert(self.property_units, "Hartrees")
        if self.analytic_derivative_order >= order:
            def func(crd, _):
                """
                **LLM Docstring**

                Evaluate the energy/gradient/Hessian (per the enclosing `order`) at a set of Cartesian coordinates using either analytic evaluation (`evaluate_term`) or finite differences (`finite_difference_derivs`), converting the result to Hartrees; for use as a `scipy.optimize`-style minimizer callback. Three variants of this closure are defined in the enclosing method depending on whether analytic evaluation, finite-difference evaluation, or a `modifier`-wrapped call applies.

                :param crd: the flattened coordinates to evaluate at
                :type crd: np.ndarray
                :param _: an extra positional argument accepted for `scipy.optimize` signature compatibility, unused
                :type _: object
                :return: the evaluated (and unit-converted) term
                :rtype: np.ndarray
                """
                res = self.evaluate_term(crd.reshape(1, -1, 3), order, **opts)
                if self.batched_orders:
                    base = res[-1] * conv
                else:
                    base = res * conv
                return base
        elif allow_fd:
            def func(crd, _):
                """
                **LLM Docstring**

                Evaluate the energy/gradient/Hessian (per the enclosing `order`) at a set of Cartesian coordinates using either analytic evaluation (`evaluate_term`) or finite differences (`finite_difference_derivs`), converting the result to Hartrees; for use as a `scipy.optimize`-style minimizer callback. Three variants of this closure are defined in the enclosing method depending on whether analytic evaluation, finite-difference evaluation, or a `modifier`-wrapped call applies.

                :param crd: the flattened coordinates to evaluate at
                :type crd: np.ndarray
                :param _: an extra positional argument accepted for `scipy.optimize` signature compatibility, unused
                :type _: object
                :return: the evaluated (and unit-converted) term
                :rtype: np.ndarray
                """
                res = self.finite_difference_derivs(crd.reshape(1, -1, 3), order, **opts)[-1]
                return res[np.newaxis] * conv # it obliterates the 1
        else:
            func = None
        if modifier is not None:
            def func(crd, _, _caller=func):
                """
                **LLM Docstring**

                Evaluate the energy/gradient/Hessian (per the enclosing `order`) at a set of Cartesian coordinates using either analytic evaluation (`evaluate_term`) or finite differences (`finite_difference_derivs`), converting the result to Hartrees; for use as a `scipy.optimize`-style minimizer callback. Three variants of this closure are defined in the enclosing method depending on whether analytic evaluation, finite-difference evaluation, or a `modifier`-wrapped call applies.

                :param crd: the flattened coordinates to evaluate at
                :type crd: np.ndarray
                :param _: an extra positional argument accepted for `scipy.optimize` signature compatibility, unused
                :type _: object
                :return: the evaluated (and unit-converted) term
                :rtype: np.ndarray
                """
                return modifier(_caller(crd, _))
        return func
    def minimizer_func(self, **opts):
        """
        **LLM Docstring**

        Build a `scipy.optimize`-compatible energy-evaluation function, via `minimizer_function_by_order(0, ...)`.

        :param opts: extra options forwarded to `minimizer_function_by_order`
        :type opts: dict
        :return: the constructed energy function
        :rtype: callable
        """
        return self.minimizer_function_by_order(0, **opts)
    def minimizer_jacobian(self, **opts):
        """
        **LLM Docstring**

        Build a `scipy.optimize`-compatible gradient-evaluation function (falling back to finite differences if needed), via `minimizer_function_by_order(1, allow_fd=True, ...)`.

        :param opts: extra options forwarded to `minimizer_function_by_order`
        :type opts: dict
        :return: the constructed gradient function
        :rtype: callable
        """
        return self.minimizer_function_by_order(1, allow_fd=True, **opts)
    def minimizer_hessian(self, **opts):
        """
        **LLM Docstring**

        Build a `scipy.optimize`-compatible Hessian-evaluation function, via `minimizer_function_by_order(2, ...)`.

        :param opts: extra options forwarded to `minimizer_function_by_order`
        :type opts: dict
        :return: the constructed Hessian function
        :rtype: callable
        """
        return self.minimizer_function_by_order(2, **opts)

    def minimizer_internal_function_by_order(self, order, allow_fd=False, modifier=None, **opts):
        """
        **LLM Docstring**

        Like `minimizer_function_by_order`, but for evaluating in internal coordinates: when re-embedding, converts the working-representation coordinates back to Cartesians before calling `evaluate_term`, re-embeds any lower-order derivatives (with the appropriate Bohr-to-`distance_units` scaling), and converts the final result to Hartrees.

        :param order: the derivative order the function should evaluate
        :type order: int
        :param allow_fd: whether to fall back to `internal_finite_difference_derivs` when the analytic order isn't high enough
        :type allow_fd: bool
        :param modifier: an optional function to post-process the returned value
        :type modifier: callable | None
        :param opts: extra options forwarded to `evaluate_term`/`internal_finite_difference_derivs`
        :type opts: dict
        :return: the constructed minimizer-compatible function, or `None` if unavailable and `allow_fd` is `False`
        :rtype: callable | None
        """
        conv = UnitsData.convert(self.property_units, "Hartrees")
        if self.analytic_derivative_order >= order:
            opts = dev.OptionsSet(opts)
            # fd_opts = opts.filter(FiniteDifferenceDerivative)
            opts = opts.exclude(FiniteDifferenceDerivative)

            reembed = dev.str_is(self.use_internals, 'reembed')
            def func(crd, _):
                """
                **LLM Docstring**

                Evaluate the energy/gradient/Hessian (per the enclosing `order`) at a set of internal (or re-embeddable Cartesian) coordinates, handling the reembedding/unembedding of coordinates and derivatives and unit conversion to Hartrees, for use as a `scipy.optimize`-style minimizer callback.

                :param crd: the coordinates to evaluate at
                :type crd: np.ndarray
                :param _: an extra positional argument accepted for `scipy.optimize` signature compatibility, unused
                :type _: object
                :return: the evaluated (and unit-converted) term
                :rtype: np.ndarray
                """
                # res = self.evaluate_term(crd, order, **opts)
                if reembed:
                    # og = crd
                    crd = self.unembed_coords(crd)
                    # print("="*100)
                    # print(np.round(og[0], 3))
                    # print(np.round(self.embed_coords(crd)[0], 3))

                    if not self.batched_orders:
                        res = [
                            self.evaluate_term(crd, o, **opts)
                            for o in range(order + 1)
                        ]
                    else:
                        res = self.evaluate_term(crd, order, **opts)
                    dist_conv = UnitsData.convert("BohrRadius", self.distance_units)
                    res = [res[0]] + self.embed_derivs(crd, [r*dist_conv**(n+1) for n,r in enumerate(res[1:])])
                    # if order > 0:
                    #     print(np.round(res[1][0], 3))
                    if not self.batched_orders:
                        res = res[-1]
                else:
                    res = self.evaluate_term(crd, order, **opts)
                if self.batched_orders:
                    return res[-1] * conv
                else:
                    return res * conv
        elif allow_fd:
            def func(crd, _):
                """
                **LLM Docstring**

                Evaluate the energy/gradient/Hessian (per the enclosing `order`) at a set of internal (or re-embeddable Cartesian) coordinates, handling the reembedding/unembedding of coordinates and derivatives and unit conversion to Hartrees, for use as a `scipy.optimize`-style minimizer callback.

                :param crd: the coordinates to evaluate at
                :type crd: np.ndarray
                :param _: an extra positional argument accepted for `scipy.optimize` signature compatibility, unused
                :type _: object
                :return: the evaluated (and unit-converted) term
                :rtype: np.ndarray
                """
                res = self.internal_finite_difference_derivs(crd, order, **opts)[-1]
                return res[np.newaxis] * conv  # it obliterates the 1
        else:
            func = None

        if modifier is not None:
            def func(crd, _, _caller=func):
                """
                **LLM Docstring**

                Evaluate the energy/gradient/Hessian (per the enclosing `order`) at a set of internal (or re-embeddable Cartesian) coordinates, handling the reembedding/unembedding of coordinates and derivatives and unit conversion to Hartrees, for use as a `scipy.optimize`-style minimizer callback.

                :param crd: the coordinates to evaluate at
                :type crd: np.ndarray
                :param _: an extra positional argument accepted for `scipy.optimize` signature compatibility, unused
                :type _: object
                :return: the evaluated (and unit-converted) term
                :rtype: np.ndarray
                """
                return modifier(_caller(crd, _))
        return func

    def get_internal_coordinate_indices(self, coord_spec):
        """
        **LLM Docstring**

        Resolve which internal-coordinate indices (and the total internal-coordinate dimension) correspond to a given coordinate specification, supporting both generic-internal (`'specs'`) and Z-matrix-based internal-coordinate systems.

        :param coord_spec: the coordinate(s) to look up, in whatever spec form the underlying internal-coordinate system expects
        :type coord_spec: object
        :return: `(ndim, inds)` -- the total number of internal coordinates and the resolved indices for `coord_spec`
        :rtype: tuple[int, list[int]]
        :raises ValueError: if the embedding's internal-coordinate spec isn't recognized
        """
        base_internals = self.embedding.internals
        if base_internals.get('specs') is not None:
            specs = base_internals.get('specs')
            ndim = len(specs)
            inds = [coordops.find_internal(specs, c) for c in coord_spec]
        elif base_internals.get('zmatrix') is not None:
            zmat = base_internals.get('zmatrix')
            ndim = coordops.num_zmatrix_coords(zmat)
            inds = coordops.zmatrix_indices(
                zmat,
                coord_spec
            )
        else:
            raise ValueError(f"can't get {coord_spec} for {base_internals}")

        return ndim, inds


    def get_internal_coordinate_projector(self, coord_spec, mask=None):
        """
        **LLM Docstring**

        Build a constraint function returning the orthogonal-projection matrix that removes the directions corresponding to a given set of internal coordinates (accounting for any redundant-coordinate transformation), for use in constrained optimization.

        :param coord_spec: the internal coordinate(s) to project out
        :type coord_spec: object
        :param mask: accepted for interface consistency but not used in this method's body
        :type mask: callable | None
        :return: a function `(coords) -> projector` broadcasting the projection matrix to the batch shape of `coords`
        :rtype: callable
        :raises ValueError: if the embedding's internal-coordinate spec isn't recognized
        """

        base_internals = self.embedding.internals
        if base_internals.get('specs') is not None:
            specs = base_internals.get('specs')
            # sidx = [tuple(s) for s in specs]
            # cpos = [tuple(c) for c in coord_spec]
            inds = [coordops.find_internal(specs, c) for c in coord_spec]
            num_coords = len(inds)
            redundant_transformation = base_internals.get('redundant_transformation')
            base_tensor = np.zeros((len(specs), num_coords))
            for n, i in enumerate(inds):
                base_tensor[..., i, n] = 1

            if redundant_transformation is not None:
                base_tensor = nput.maximum_similarity_transformation(
                    redundant_transformation,
                    base_tensor,
                    apply_transformation=True
                )

            proj = nput.orthogonal_projection_matrix(base_tensor)
            def constraint(coords):
                """
                **LLM Docstring**

                Broadcast the enclosing scope's precomputed orthogonal-projection matrix to match the batch shape of the given coordinates.

                :param coords: the coordinates whose batch shape the projector should be broadcast to
                :type coords: np.ndarray
                :return: the broadcast projection matrix
                :rtype: np.ndarray
                """
                nstruct = len(coords.shape[:-1])
                return np.broadcast_to(
                    np.expand_dims(proj, list(range(nstruct))),
                    coords.shape[:-1] + proj.shape
                )

        elif base_internals.get('zmatrix') is not None:
            zmat = base_internals.get('zmatrix')
            num_specs = coordops.num_zmatrix_coords(
                zmat,
                coord_spec
            )
            inds = coordops.zmatrix_indices(
                zmat,
                coord_spec
            )
            num_coords = len(inds)

            base_tensor = np.zeros((num_specs, num_coords))
            for n, i in enumerate(inds):
                base_tensor[..., i, n] = 1
            proj = nput.orthogonal_projection_matrix(base_tensor)

            def constraint(coords):
                """
                **LLM Docstring**

                Broadcast the enclosing scope's precomputed orthogonal-projection matrix to match the batch shape of the given coordinates.

                :param coords: the coordinates whose batch shape the projector should be broadcast to
                :type coords: np.ndarray
                :return: the broadcast projection matrix
                :rtype: np.ndarray
                """
                nstruct = len(coords.shape[:-1])
                return np.broadcast_to(
                    np.expand_dims(proj, list(range(nstruct))),
                    coords.shape[:-1] + proj.shape
                )
        else:
            raise ValueError(f"can't get {coord_spec} for {base_internals}")

        return constraint

    def get_internal_coordinate_constraints(self, constraints):
        """
        **LLM Docstring**

        Resolve a constraints specification for internal coordinates into a projector-generating constraint function: a dict of per-coordinate mask functions, a bare list of coordinates, a raw projection-matrix-like array, or an already-resolved constraint object passed through unchanged.

        :param constraints: the constraints specification to resolve
        :type constraints: dict | list | np.ndarray | object
        :return: the resolved constraint function/matrix
        :rtype: callable | np.ndarray | object
        """
        if dev.is_dict_like(constraints):
            coord_list = list(constraints.keys())
            funs = list(constraints.values())

            def mask_func(internal_vals):
                """
                **LLM Docstring**

                Apply each per-coordinate mask function to the corresponding internal-coordinate value, from the enclosing `constraints` dict.

                :param internal_vals: the internal-coordinate values to test, aligned with the enclosing `funs` list
                :type internal_vals: Iterable
                :return: a generator of the per-coordinate mask results
                :rtype: Iterator
                """
                return (f(i) for i, f in zip(internal_vals, funs))

            return self.get_internal_coordinate_projector(coord_list, mask=mask_func)
        elif coordops.is_coordinate_list_like(constraints):
            return self.get_internal_coordinate_projector(constraints)
        elif nput.is_numeric_array_like(constraints):
            return nput.orthogonal_projection_matrix(constraints)
        else:
            return constraints
    internal_optimizer_defaults = dict(
        max_displacement=.01
    )

    optimizer_defaults = dict(
        method='quasi-newton',
        unitary=False,
        # generate_rotation=False,
        # dtype='float64',
        orthogonal_directions=None,
        convergence_metric=None,
        tol=1e-8,
        max_iterations=25,
        damping_parameter=None,
        damping_exponent=None,
        restart_interval=None,
        max_displacement=.2,
        line_search=None,
        optimizer_settings=None,
        mode=None,
        func=None,
        jacobian=None,
        hessian=None,
        logger=None,

        generate_rotation=False,
        dtype='float64',
        orthogonal_projection_generator=None,
        region_constraints=None,
        max_displacement_norm=None,
        oscillation_damping_factor=None,
        termination_function=None,
        prevent_oscillations=None,
        use_max_for_error=True,
        track_best=True
    )
    @classmethod
    def get_optimizer_options(self):
        """
        **LLM Docstring**

        The full set of recognized keyword-argument names for `optimize`/`optimize_iterative`/`relaxed_scan`, combining the optimizer defaults, a few extra named options, and the `FiniteDifferenceDerivative` option names -- used to split a caller's `**opts` between optimizer-specific and evaluator-specific options.

        :return: the tuple of recognized option names
        :rtype: tuple[str]
        """
        # import scipy.optimize._optimize
        return (
                tuple(self.optimizer_defaults.keys())
                + ('coordinate_constraints', 'gradient_modification_function', "initialization_function", "return_trajectory")
                + FiniteDifferenceDerivative.__props__
        )

    def get_coordinate_projector(self, coord_spec, region_constraints=None, mask=None):
        """
        **LLM Docstring**

        Build a constraint function returning the orthogonal-projection matrix that removes the directions spanned by a general (Cartesian-space) coordinate specification, optionally restricting the projection to a specific region (via `region_constraints`) or an arbitrary mask function.

        :param coord_spec: the coordinate specification (e.g. bond/angle/dihedral definitions) whose directions should be projected out
        :type coord_spec: object
        :param region_constraints: per-coordinate `(min, max)` bounds restricting where the projection is applied
        :type region_constraints: Iterable[tuple] | None
        :param mask: an explicit mask function selecting which basis vectors to zero out, in addition to any region mask
        :type mask: callable | None
        :return: a function `(coords) -> projector` computing the (region-restricted) orthogonal projection matrix at the given coordinates
        :rtype: callable
        """
        tf_fun = nput.internal_conversion_function(coord_spec, order=0)
        if region_constraints is not None:
            def region_mask(coords):
                """
                **LLM Docstring**

                Check whether each coordinate value in `coords` falls within its corresponding `region_constraints` bounds.

                :param coords: the coordinate values to check
                :type coords: np.ndarray
                :return: a boolean mask of the same shape, `True` where the constraint is satisfied
                :rtype: np.ndarray
                """
                mask_success = np.full(coords.shape, True)
                for i,(m,M) in enumerate(region_constraints):
                    mask_success[..., i] = m <= coords[..., i] & coords[..., i] <= M
                return mask_success
            if mask is None:
                mask = region_mask
            else:
                def mask2(coords, mask=mask, region_mask=region_mask):
                    """
                    **LLM Docstring**

                    Combine an existing mask function with `region_mask`, requiring both to pass.

                    :param coords: the coordinate values to check
                    :type coords: np.ndarray
                    :return: the combined boolean mask
                    :rtype: np.ndarray
                    """
                    return mask(coords) & region_mask(coords)
                mask = mask2

        order = None if (mask is None and region_constraints is None) else 1
        def constraint(coords):
            """
            **LLM Docstring**

            Compute the (possibly region/mask-restricted) orthogonal-projection matrix for the enclosing coordinate specification at a batch of Cartesian geometries.

            :param coords: the Cartesian geometries to compute the projector at
            :type coords: np.ndarray
            :return: the projection matrix/matrices
            :rtype: np.ndarray
            """
            coords = coords.reshape((-1, len(self.embedding.masses), 3))
            bases, _, _ = nput.internal_basis(coords, coord_spec)#, order=order)
            basis_vectors = np.concatenate(bases, axis=-1)
            if mask is not None:
                base_tensors = tf_fun(coords)
                checks = ~mask(base_tensors[0])
                basis_vectors[..., checks] = 0 # no orthogonal projection when within the region

            # print(bases[0])
            proj = nput.orthogonal_projection_matrix(basis_vectors, orthonormal=False)
            # print(proj)
            # raise Exception(...)

            return proj

        return constraint

    def get_coordinate_constraints(self, constraints, region_constraints=None):
        """
        **LLM Docstring**

        Resolve a constraints specification for general (Cartesian-space) coordinates into a projector-generating constraint function: a dict of per-coordinate mask functions, a bare coordinate-list specification, a raw projection-matrix-like array, or an already-resolved constraint object passed through unchanged.

        :param constraints: the constraints specification to resolve
        :type constraints: dict | list | np.ndarray | object
        :param region_constraints: per-coordinate bounds forwarded to `get_coordinate_projector` when `constraints` is a coordinate-list spec
        :type region_constraints: Iterable[tuple] | None
        :return: the resolved constraint function/matrix
        :rtype: callable | np.ndarray | object
        """
        if dev.is_dict_like(constraints):
            coord_list = list(constraints.keys())
            funs = list(constraints.values())
            def mask_func(internal_vals):
                """
                **LLM Docstring**

                Apply each per-coordinate mask function to the corresponding coordinate value, from the enclosing `constraints` dict.

                :param internal_vals: the coordinate values to test, aligned with the enclosing `funs` list
                :type internal_vals: Iterable
                :return: a generator of the per-coordinate mask results
                :rtype: Iterator
                """
                return (f(i) for i,f in zip(internal_vals, funs))

            return self.get_coordinate_projector(coord_list, mask=mask_func)
        elif coordops.is_coordinate_list_like(constraints):
            return self.get_coordinate_projector(constraints, region_constraints=region_constraints)
        elif nput.is_numeric_array_like(constraints):
            return nput.orthogonal_projection_matrix(constraints)
        else:
            return constraints

    def _modify_gradient(self,
                         gradient_function,
                         modification_function,
                         modification_mode='shift',
                         convert_modification_distances=True,
                         convert_modification_energies=True,
                         orthogonal_projector=None,
                         orthogonal_projection_generator=None,
                         use_forces=False,
                         coord_prep=None,
                         grad_prep=None,
                         grad_post=None
                         ):
        """
        **LLM Docstring**

        Build a wrapped version of a raw gradient-computing function that applies a user-supplied modification function to the gradient/forces (in `'shift'` mode, adding the modification's output to the raw gradient with appropriate unit conversions; otherwise, replacing the gradient with the modification function's output directly), then optionally projects the result through an orthogonal projector/projector-generator, converting units as needed on the way in and out.

        :param gradient_function: the underlying function that computes the raw gradient/forces
        :type gradient_function: callable
        :param modification_function: the function used to modify the gradient, called as `(coords_in_bohr, gradient)` (or with a sign flip if `use_forces`)
        :type modification_function: callable
        :param modification_mode: `'shift'` to add the modification's output to the raw gradient, or any other value to replace the gradient with it directly
        :type modification_mode: str
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function (`'shift'` mode only)
        :type convert_modification_energies: bool
        :param orthogonal_projector: a fixed projection matrix to apply to the modified gradient
        :type orthogonal_projector: np.ndarray | None
        :param orthogonal_projection_generator: a per-geometry projector-generating function, combined with `orthogonal_projector` if both are given
        :type orthogonal_projection_generator: callable | None
        :param use_forces: whether the wrapped function operates on forces (sign-flipped gradients) rather than raw gradients
        :type use_forces: bool
        :param coord_prep: a function to extract/prepare the coordinates from the wrapped function's raw arguments
        :type coord_prep: callable | None
        :param grad_prep: a function to extract the raw gradient (and any extra state) from the wrapped function's result, or signal a short-circuit by returning `(None, state)`
        :type grad_prep: callable | None
        :param grad_post: a function to reinsert the modified gradient back into the original result structure
        :type grad_post: callable | None
        :return: the wrapped gradient-computing function
        :rtype: callable
        :raises ValueError: if `modification_mode` is `None`
        """
        if modification_mode is None:
            raise ValueError("no gradient modification mode supplied")
        if convert_modification_distances:
            d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
        else:
            d_conv = 1
        if orthogonal_projector is not None:
            orthogonal_projector = np.asarray(orthogonal_projector)
        if modification_mode == 'shift':
            if convert_modification_energies:
                e_conv = UnitsData.convert(self.property_units, "Hartrees")
            else:
                e_conv = 1
            def grad(crds, *rem):
                """
                **LLM Docstring**

                Compute the (possibly modification-function-adjusted and orthogonally-projected) gradient/forces by calling the enclosing `gradient_function`, applying `modification_function` per `modification_mode`, and applying any orthogonal projector, converting units as configured in the enclosing scope.

                :param crds: the coordinates to evaluate at
                :type crds: np.ndarray
                :param rem: any additional positional arguments forwarded to `gradient_function`
                :type rem: tuple
                :return: the modified gradient/forces (or the short-circuited state, if `grad_prep` signaled nothing new to modify)
                :rtype: np.ndarray | object
                """
                res = gradient_function(crds, *rem)
                if coord_prep is not None:
                    crds = coord_prep(crds, *rem)
                if grad_prep is not None:
                    res, state = grad_prep(res)
                    if res is None: return state
                else:
                    state = None
                re_res = res * e_conv / d_conv
                if use_forces:
                    supp = -modification_function(crds * d_conv, -re_res)
                else:
                    supp = modification_function(crds * d_conv, re_res)
                res = re_res + supp
                if orthogonal_projector is not None:
                    projector = orthogonal_projector[np.newaxis]
                    if orthogonal_projection_generator is not None:
                        projector = projector @ orthogonal_projection_generator(crds)
                elif orthogonal_projection_generator is not None:
                    projector = orthogonal_projection_generator(crds)
                else:
                    projector = None

                if projector is not None:
                    res_shape = res.shape
                    res = res.reshape(projector.shape[:-2] + (-1, 1))
                    res = projector @ res
                    res = res.reshape(res_shape)
                res = res * d_conv / e_conv
                if grad_post is not None:
                    res = grad_post(res, state)
                return res
        else:
            def grad(crds, *rem):
                """
                **LLM Docstring**

                Compute the (possibly modification-function-adjusted and orthogonally-projected) gradient/forces by calling the enclosing `gradient_function`, applying `modification_function` per `modification_mode`, and applying any orthogonal projector, converting units as configured in the enclosing scope.

                :param crds: the coordinates to evaluate at
                :type crds: np.ndarray
                :param rem: any additional positional arguments forwarded to `gradient_function`
                :type rem: tuple
                :return: the modified gradient/forces (or the short-circuited state, if `grad_prep` signaled nothing new to modify)
                :rtype: np.ndarray | object
                """
                res = gradient_function(crds, *rem)
                if coord_prep is not None:
                    crds = coord_prep(crds, *rem)
                if grad_prep is not None:
                    res, state = grad_prep(res)
                    if res is None: return state
                else:
                    state = None
                if use_forces:
                    res = -modification_function(crds * d_conv, -res)
                else:
                    res = modification_function(crds * d_conv, res)
                if orthogonal_projector is not None:
                    projector = orthogonal_projector[np.newaxis]
                    if orthogonal_projection_generator is not None:
                        projector = projector @ orthogonal_projection_generator(crds)
                elif orthogonal_projection_generator is not None:
                    projector = orthogonal_projection_generator(crds)
                else:
                    projector = None
                if projector is not None:
                    res = projector @ res
                res = res * d_conv
                if grad_post is not None:
                    res = grad_post(res, state)
                return res
        return grad

    scipy_no_hessian_methods = {'cg', 'bfgs'}
    scipy_no_grad_methods = {'nelder-mead'}
    use_scipy_linesearch = False
    def optimize_iterative(self,
                           coords,
                           coordinate_constraints=None,
                           gradient_modification_function=None,
                           gradient_modification_mode='shift',
                           convert_modification_distances=True,
                           convert_modification_energies=True,
                           return_trajectory=False,
                           initialization_function=None,
                           line_search_step=None,
                           **opts
                           ):
        """
        **LLM Docstring**

        Core geometry-optimization driver supporting several backends (a custom iterative step-minimization scheme by default, or delegating to `scipy.optimize.minimize`, ASE's optimizers, or pysisyphus, depending on `mode`), building the energy/gradient/Hessian callables (falling back to the evaluator's own `minimizer_func`/`minimizer_jacobian`/`minimizer_hessian` if not supplied), applying any coordinate constraints as an orthogonal projector, and optionally recording/returning the optimization trajectory.

        :param coords: the initial coordinates to optimize from
        :type coords: np.ndarray
        :param coordinate_constraints: coordinates to hold fixed/project out during optimization
        :type coordinate_constraints: object | None
        :param gradient_modification_function: an optional function to modify the gradient/forces during optimization (e.g. to bias the search)
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode forwarded to `_modify_gradient`
        :type gradient_modification_mode: str
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool
        :param return_trajectory: whether to record and return the optimization trajectory
        :type return_trajectory: bool
        :param initialization_function: an optional function applied to the initial coordinates (in Bohr) before optimization begins
        :type initialization_function: callable | None
        :param line_search_step: a fixed step size to force during (non-`'scipy'`-line-search) line searches, bypassing the Wolfe-condition search
        :type line_search_step: float | None
        :param opts: the full set of optimizer options (method, tolerances, damping, region/coordinate constraints, backend-specific settings, etc.), merged with `self.optimizer_defaults`
        :type opts: dict
        :return: `(converged, coords_or_(coords,trajectory), info)` -- whether the optimization converged, the optimized (and optionally trajectory) coordinates, and backend-specific optimization metadata
        :rtype: tuple
        """
        from McUtils.Numputils import iterative_step_minimize
        opts = dict(self.optimizer_defaults, **{o:k for o,k in opts.items() if k is not None})

        (
            method,
            unitary,
            # generate_rotation=False,
            # dtype='float64',
            orthogonal_directions,
            orthogonal_projection_generator,
            convergence_metric,
            tol,
            max_iterations,
            damping_parameter,
            damping_exponent,
            restart_interval,
            region_constraints,
            max_displacement,
            line_search,
            optimizer_settings,
            mode,
            func,
            jacobian,
            hessian,
            logger
        ) = (
            opts.pop(k) for k in
            [
                "method",
                "unitary",
                # generate_rotation=False,
                # dtype='float64',
                "orthogonal_directions",
                "orthogonal_projection_generator",
                "convergence_metric",
                "tol",
                "max_iterations",
                "damping_parameter",
                "damping_exponent",
                "restart_interval",
                "region_constraints",
                "max_displacement",
                "line_search",
                "optimizer_settings",
                "mode",
                "func",
                "jacobian",
                "hessian",
                "logger"
            ]
        )

        fopts = dev.OptionsSet(opts).exclude(None, props=self.optimizer_defaults.keys())
        opt_opts = dev.OptionsSet(opts).filter(None, props=self.optimizer_defaults.keys())

        if optimizer_settings is None:
            optimizer_settings = {}

        #TODO: let evaluations cache results when batched
        if func is None:
            func = self.minimizer_func(**fopts)
        if jacobian is None:
            jacobian = self.minimizer_jacobian(**fopts)
        if hessian is None:
            hessian = self.minimizer_hessian(**fopts)

        logger = Logger.lookup(logger, construct=True)

        if orthogonal_projection_generator is None:
            orthogonal_projection_generator = self.get_coordinate_constraints(
                coordinate_constraints,
                region_constraints
            )

        if initialization_function is not None:
            if convert_modification_distances:
                d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
            else:
                d_conv = 1
            coords = initialization_function(d_conv*coords)
            coords = coords / d_conv
        if mode == 'scipy' or method == 'scipy':
            from scipy.optimize import minimize, _optimize, _minimize

            if line_search is None:
                line_search = self.use_scipy_linesearch

            if not line_search:
                optimizer_settings = {'c1': 0.00001, 'c2': 0.999} | optimizer_settings

            def sfunc(crd):
                """
                **LLM Docstring**

                Evaluate the enclosing energy function at a flattened coordinate vector and unwrap a length-1 batch result, for use as a `scipy.optimize.minimize` objective function.

                :param crd: the flattened coordinates to evaluate at
                :type crd: np.ndarray
                :return: the scalar energy value
                :rtype: float
                """
                wat = func(crd, None)
                if isinstance(wat, np.ndarray):
                    return wat[0]
                else:
                    return wat
            def sjacobian(crd):
                """
                **LLM Docstring**

                Evaluate the enclosing gradient function at a flattened coordinate vector, apply any orthogonal projection, and unwrap a length-1 batch result, for use as a `scipy.optimize.minimize` Jacobian.

                :param crd: the flattened coordinates to evaluate at
                :type crd: np.ndarray
                :return: the (projected) gradient vector
                :rtype: np.ndarray
                """
                huh = jacobian(crd, None)
                if orthogonal_projection_generator is not None:
                    gen = orthogonal_projection_generator(crd)
                    huh = np.reshape(gen @ huh[..., np.newaxis], huh.shape)
                # print(huh)
                if huh.ndim == 2:
                    res = huh[0]
                else:
                    res = huh
                return res

            if gradient_modification_function is not None:
                sjacobian = self._modify_gradient(sjacobian,
                                                  gradient_modification_function,
                                                  modification_mode=gradient_modification_mode,
                                                  convert_modification_distances=convert_modification_distances,
                                                  convert_modification_energies=convert_modification_energies,
                                                  )

            if self.analytic_derivative_order > 1:
                def shessian(crd):
                    """
                    **LLM Docstring**

                    Evaluate the enclosing Hessian function at a flattened coordinate vector and unwrap a length-1 batch result, for use as a `scipy.optimize.minimize` Hessian.

                    :param crd: the flattened coordinates to evaluate at
                    :type crd: np.ndarray
                    :return: the Hessian matrix
                    :rtype: np.ndarray
                    """
                    return hessian(crd, None)[0]
            else:
                shessian = None

            callback = optimizer_settings.pop('callback', None)
            if return_trajectory:
                traj = []
                if callback is None:
                    def callback(x):
                        """
                        **LLM Docstring**

                        Record each optimizer iterate into the enclosing trajectory list, for use as a `scipy.optimize.minimize` callback when `return_trajectory` is set and no user callback was supplied.

                        :param x: the current iterate's coordinates
                        :type x: np.ndarray
                        :return: None
                        :rtype: None
                        """
                        traj.append(x)
                else:
                    def append_callback(intermediate_result, cb):
                        """
                        **LLM Docstring**

                        Record each optimizer iterate into the enclosing trajectory list, then forward the call on to a user-supplied callback, combining trajectory-recording with the user's own callback behavior.

                        :param intermediate_result: the current optimizer iterate (as passed by `scipy.optimize.minimize`)
                        :type intermediate_result: object
                        :param cb: the user-supplied callback to forward to
                        :type cb: callable
                        :return: the result of calling `cb(intermediate_result)`
                        :rtype: object
                        """
                        traj.append(intermediate_result)
                        return cb(intermediate_result)
                    callback = functools.partial(append_callback, cb=callback)
            if logger.active:
                prev_re=[coords.flatten().view(np.ndarray)]
                if callback is None:
                    def log_callback(intermediate_result, prev_re):
                        """
                        **LLM Docstring**

                        Log the current optimizer iterate and its step relative to the previous one via the enclosing `logger`, updating the running history list; two variants of this closure are defined depending on whether a user callback also needs to be chained.

                        :param intermediate_result: the current optimizer iterate
                        :type intermediate_result: object
                        :param cb: (one variant only) a user-supplied callback to additionally invoke
                        :type cb: callable
                        :param prev_re: the running list of previous iterates, from the enclosing scope
                        :type prev_re: list
                        :return: `None`, or the result of calling `cb(intermediate_result)` in the variant that chains a user callback
                        :rtype: None | object
                        """
                        logger.log_print(
                            [
                                "Struct: {intermediate_result}",
                                "Step: {intermediate_step}"
                            ],
                            intermediate_result=intermediate_result,
                            intermediate_step=intermediate_result - prev_re[-1]
                        ),
                        prev_re.append(intermediate_result)
                    callback = functools.partial(log_callback, prev_re=prev_re)
                else:
                    def log_callback(intermediate_result, cb, prev_re):
                        """
                        **LLM Docstring**

                        Log the current optimizer iterate and its step relative to the previous one via the enclosing `logger`, updating the running history list; two variants of this closure are defined depending on whether a user callback also needs to be chained.

                        :param intermediate_result: the current optimizer iterate
                        :type intermediate_result: object
                        :param cb: (one variant only) a user-supplied callback to additionally invoke
                        :type cb: callable
                        :param prev_re: the running list of previous iterates, from the enclosing scope
                        :type prev_re: list
                        :return: `None`, or the result of calling `cb(intermediate_result)` in the variant that chains a user callback
                        :rtype: None | object
                        """
                        logger.log_print(
                            [
                                "Struct: {intermediate_result}",
                                "Step: {intermediate_step}"
                            ],
                            intermediate_result=intermediate_result,
                            intermediate_step=intermediate_result - prev_re[-1]
                        ),
                        prev_re.append(intermediate_result)
                        return cb(intermediate_result)
                    callback = functools.partial(log_callback, cb=callback, prev_re=prev_re)

            method = (
                method
                    if dev.str_is(mode, 'scipy') else
                opts.pop('scipy_method', self.optimizer_defaults.get('method', 'bfgs'))
            )
            min_ops = {'maxiter':max_iterations}
            method = 'bfgs' if method=='quasi-newton' else method
            if method in self.scipy_no_hessian_methods or method in self.scipy_no_grad_methods:
                shessian = None
            if method in self.scipy_no_grad_methods:
                sjacobian = None
                if coordinate_constraints is not None:
                    base_internals = self.embedding.internals
                    if base_internals.get('specs') is not None:
                        specs = base_internals.get('specs')
                        inds = [coordops.find_internal(specs, c) for c in coordinate_constraints]

                    elif base_internals.get('zmatrix') is not None:
                        zmat = base_internals.get('zmatrix')
                        # num_specs = coordops.num_zmatrix_coords(
                        #     zmat,
                        #     coordinate_constraints
                        # )
                        inds = coordops.zmatrix_indices(
                            zmat,
                            coordinate_constraints
                        )

                    else:
                        raise ValueError(f"can't get {coordinate_constraints} for {base_internals}")

                    min_ops['bounds'] = [
                        (c, c)
                            if i in inds else
                        (None, None)
                        for i,c in enumerate(coords.flatten())
                    ]

                if region_constraints is not None:
                    base_internals = self.embedding.internals
                    if base_internals.get('specs') is not None:
                        specs = base_internals.get('specs')
                        cons = {
                            coordops.find_internal(specs, c):v
                            for c,v in region_constraints.items()
                        }

                    elif base_internals.get('zmatrix') is not None:
                        zmat = base_internals.get('zmatrix')
                        # num_specs = coordops.num_zmatrix_coords(
                        #     zmat,
                        #     coordinate_constraints
                        # )
                        inds = coordops.zmatrix_indices(
                            zmat,
                            list(region_constraints.keys())
                        )
                        cons = {
                            i: v
                            for i, v in zip(inds, region_constraints.values())
                        }
                    else:
                        raise ValueError(f"can't get {coordinate_constraints} for {base_internals}")

                    #TODO: decide if I want to apply the region bounds as an offset or not...
                    min_ops['bounds'] = [
                        (
                            (c + cons[i][0] if b1 is None else b1)
                            (c + cons[i][1] if b1 is None else b2)
                        )
                            if i in cons else
                        (b1, b2)
                        for i,(c, (b1, b2)) in enumerate(zip(
                            coords.flatten(),
                            min_ops.get('bonds', [None] * len(coords))
                        ))
                    ]

            scipy_meth = 'bfgs' if method=='quasi-newton' else method
            if 'bfgs' in scipy_meth or scipy_meth in {'cg'}:
                min_ops['gtol'] = tol
                # min_ops['ftol'] = 0
                # min_ops['xtol'] = 0
            min_ops = dict(
                min_ops,
                **optimizer_settings
            )
            bounds = min_ops.pop('bounds', None)
            constraints = min_ops.pop('constraints', ())
            try:
                if not line_search:
                    old_wolfe = _optimize.line_search_wolfe1
                    def _find_max_displacement_step(
                            f, fprime, xk, pk, gfk=None,
                            old_fval=None, old_old_fval=None,
                            args=(), c1=1e-4, c2=0.9, amax=50, amin=1e-8,
                            xtol=1e-14
                    ):
                        """
                        **LLM Docstring**

                        Replacement for `scipy.optimize`'s internal Wolfe-condition line search that instead returns either a fixed `line_search_step` or a step size capping the maximum per-coordinate displacement at `max_displacement`, used to prevent scipy's default line search from taking overly large steps.

                        :param f: the objective function (unused, accepted for interface compatibility with `scipy.optimize._optimize.line_search_wolfe1`)
                        :type f: callable
                        :param fprime: the gradient function (unused)
                        :type fprime: callable
                        :param xk: the current point (unused)
                        :type xk: np.ndarray
                        :param pk: the search direction, used to compute the maximum per-coordinate step
                        :type pk: np.ndarray
                        :param gfk: the current gradient (unused)
                        :type gfk: np.ndarray | None
                        :param old_fval: the previous function value, passed through unchanged
                        :type old_fval: float | None
                        :param old_old_fval: the function value before that, passed through unchanged
                        :type old_old_fval: float | None
                        :param args: extra arguments (unused)
                        :type args: tuple
                        :param c1: Wolfe condition parameter (unused)
                        :type c1: float
                        :param c2: Wolfe condition parameter (unused)
                        :type c2: float
                        :param amax: maximum step size (unused)
                        :type amax: float
                        :param amin: minimum step size (unused)
                        :type amin: float
                        :param xtol: tolerance (unused)
                        :type xtol: float
                        :return: `(step_size, None, None, old_fval, old_old_fval, None)`, matching the tuple shape `scipy.optimize`'s line search would return
                        :rtype: tuple
                        """
                        if line_search_step is not None:
                            return line_search_step
                        else:
                            max_pk = np.max(np.abs(pk))
                            if max_pk < 1e-6:
                                return max_displacement, None, None, old_fval, old_old_fval, None
                            else:
                                return max_displacement/max_pk, None, None, old_fval, old_old_fval, None
                    _optimize.line_search_wolfe1 = _find_max_displacement_step
                min = minimize(sfunc,
                               coords.flatten(),
                               method=scipy_meth,
                               tol=tol,
                               jac=sjacobian,
                               hess=shessian,
                               callback=callback,
                               bounds=bounds,
                               constraints=constraints,
                               options=min_ops
                               )
            finally:
                if not line_search:
                    _optimize.line_search_wolfe1 = old_wolfe
            if return_trajectory:
                return min.success, (
                    min.x.reshape(coords.shape),
                    [t.reshape(coords.shape) for t in traj],
                ), min
            else:
                return min.success, min.x.reshape(coords.shape), min
        elif mode == 'ase':
            from McUtils.ExternalPrograms import ASEMolecule
            # calc = self.to_ase()

            if orthogonal_projection_generator is None:
                orthogonal_projection_generator = self.get_coordinate_constraints(
                    coordinate_constraints,
                    region_constraints
                )

            calc = self.to_ase(
                gradient_modification_function=gradient_modification_function,
                gradient_modification_mode=gradient_modification_mode,
                convert_modification_distances=convert_modification_distances,
                convert_modification_energies=convert_modification_energies,
                orthogonal_projection_generator=orthogonal_projection_generator
            )

            if initialization_function is not None:
                d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
                coords = initialization_function(d_conv * coords)
                coords = coords / d_conv

            if return_trajectory:
                traj = tempfile.NamedTemporaryFile().name
            else:
                traj = None

            mol = ASEMolecule.from_coords(
                self.atoms,
                coords.reshape((-1,) + coords.shape[-2:])[0],
                calculator=calc
            )
            new_opt = mol.optimize_structure(coords,
                                             fmax=tol, steps=max_iterations,
                                             maxstep=max_displacement,
                                             trajectory=traj,
                                             logger=logger,
                                             method=method)
            if return_trajectory:
                from ase.io.trajectory import Trajectory
                try:
                    with Trajectory(traj, mode='r') as reader:
                        traj_data = [
                            a.positions for a in reader
                        ]
                except FileNotFoundError:
                    ...
                else:
                    os.remove(traj)
                cond, coords, d = new_opt
                new_opt = (cond, (coords, traj_data), d)
            return new_opt
        elif mode == 'pysis':
            from McUtils.ExternalPrograms import run_pysisyphus, prep_pysis_images, patch_pysis_logging
            # calc = self.to_ase()
            patch_pysis_logging()

            if orthogonal_projection_generator is None:
                orthogonal_projection_generator = self.get_coordinate_constraints(
                    coordinate_constraints,
                    region_constraints
                )

            calc = self.to_pysis(
                gradient_modification_function=gradient_modification_function,
                gradient_modification_mode=gradient_modification_mode,
                convert_modification_distances=False,
                convert_modification_energies=False,
                orthogonal_projection_generator=orthogonal_projection_generator
            )

            if initialization_function is not None:
                d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
                coords = initialization_function(d_conv * coords)
                coords = coords / d_conv

            if return_trajectory:
                traj = tempfile.NamedTemporaryFile().name
            else:
                traj = None

            geom = prep_pysis_images(
                self.atoms,
                coords * UnitsData.convert(self.distance_units, "BohrRadius")
            )
            geom.set_calculator(calc)
            use_max_for_error = opts.pop('use_max_for_error', None)
            if method == 'ts':
                method, optimizer = method, None
                opts['images'] = [geom]
            elif method.startswith('ts-'):
                method, optimizer = method.split("-", 1)
                opts['images'] = [geom]
            else:
                method, optimizer = 'optimize', method
                opts['geom'] = geom
            geom, optimizer, logs = run_pysisyphus(
                None,
                method,
                optimizer=optimizer,
                max_cycles=max_iterations,
                max_step=max_displacement,
                tol=tol,
                use_max_for_error=use_max_for_error,
                return_logs=return_trajectory,
                logger=logger,
                **opts
            )
            # new_opt = mol.optimize_structure(coords,
            #                                  fmax=tol, steps=max_iterations,
            #                                  maxstep=max_displacement,
            #                                  trajectory=traj,
            #                                  method=method)
            # if return_trajectory:
            #     new_opt,
            #     from ase.io.trajectory import Trajectory
            #     try:
            #         with Trajectory(traj, mode='r') as reader:
            #             traj_data = [
            #                 a.positions for a in reader
            #             ]
            #     except FileNotFoundError:
            #         ...
            #     else:
            #         os.remove(traj)
            #     cond, coords, d = new_opt
            #     new_opt = (cond, (coords, traj_data), d)
            coords = geom.cart_coords.reshape(coords.shape)
            if return_trajectory:
                raise NotImplementedError("kinda annoying")
            else:
                return optimizer.is_converged, coords * UnitsData.convert("BohrRadius", self.distance_units), {}
        else:
            if gradient_modification_function is not None:
                jacobian = self._modify_gradient(jacobian,
                                                 gradient_modification_function,
                                                 modification_mode=gradient_modification_mode,
                                                 convert_modification_distances=False,#convert_modification_distances,
                                                 convert_modification_energies=False,#convert_modification_energies,
                                                 )

            if method is None or isinstance(method, str):
                method = {
                        'method': method,
                        'func': func,
                        'jacobian': jacobian,
                        'hessian': hessian,
                        'damping_parameter': damping_parameter,
                        'damping_exponent': damping_exponent,
                        'restart_interval': restart_interval,
                        'line_search': line_search,
                    }
                method = {k:v for k,v in method.items() if v is not None}
                method = dict(method, **optimizer_settings)
            opt_data = iterative_step_minimize(
                coords.flatten(),
                method,
                unitary=unitary,
                orthogonal_directions=orthogonal_directions,
                convergence_metric=convergence_metric,
                tol=tol,
                function=func,
                max_iterations=max_iterations,
                max_displacement=max_displacement,
                logger=logger,
                orthogonal_projection_generator=orthogonal_projection_generator,
                return_trajectory=return_trajectory,
                region_constraints=region_constraints,
                **opt_opts
            )
            if return_trajectory:
                (opt_coords, converged, (errs, its)), traj = opt_data
                return (
                    converged,
                    (opt_coords.reshape(coords.shape), [t.reshape(coords.shape) for _,t in traj]),
                    {'errors': errs, 'iterations': its}
                )
            else:
                opt_coords, converged, (errs, its) = opt_data
                return converged, opt_coords.reshape(coords.shape), {'errors':errs, 'iterations':its}


    def optimize(self,
                 coords,
                 use_internals=None,
                 func=None,
                 jacobian=None,
                 hessian=None,
                 logger=None,
                 coordinate_constraints=None,
                 orthogonal_projection_generator=None,
                 gradient_modification_function=None,
                 initialization_function=None,
                 return_trajectory=False,
                 mode=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Geometry-optimize a structure, automatically routing through internal coordinates (building the appropriate internal-coordinate minimizer functions and constraints, then converting the result back to Cartesians) when `use_internals` is enabled and internal coordinates are available and the backend `mode` supports it, otherwise delegating directly to `optimize_iterative`.

        :param coords: the initial coordinates to optimize from
        :type coords: np.ndarray
        :param use_internals: whether to optimize in internal coordinates; defaults to `self.use_internals`
        :type use_internals: bool | None
        :param func: an explicit energy function to use instead of the default minimizer function
        :type func: callable | None
        :param jacobian: an explicit gradient function to use
        :type jacobian: callable | None
        :param hessian: an explicit Hessian function to use
        :type hessian: callable | None
        :param logger: logger to report optimization progress to
        :type logger: Logger | str | None
        :param coordinate_constraints: coordinates to hold fixed/project out during optimization
        :type coordinate_constraints: object | None
        :param orthogonal_projection_generator: an explicit per-geometry constraint projector, bypassing automatic construction from `coordinate_constraints`
        :type orthogonal_projection_generator: callable | None
        :param gradient_modification_function: an optional function to modify the gradient during optimization
        :type gradient_modification_function: callable | None
        :param initialization_function: an optional function applied to the initial coordinates before optimization
        :type initialization_function: callable | None
        :param return_trajectory: whether to record and return the optimization trajectory
        :type return_trajectory: bool
        :param mode: the optimization backend to use (`'scipy'`, `'ase'`, `'pysis'`, or `None` for the default iterative scheme); internal-coordinate optimization is skipped for `'ase'`/`'pysis'`
        :type mode: str | None
        :param opts: extra options forwarded to `optimize_iterative`
        :type opts: dict
        :return: `(converged, coords_or_(coords,trajectory), info)`
        :rtype: tuple
        """

        if use_internals is None: use_internals = self.use_internals
        if (
                mode not in {'pysis', 'ase'}
                and use_internals
                and self.embedding.internals is not None
        ):
            opts = dict(self.internal_optimizer_defaults, **opts)
            fopts = dev.OptionsSet(opts).exclude(None, props=self.optimizer_defaults.keys())
            if func is None:
                func = self.minimizer_internal_function_by_order(0, **fopts)
            if jacobian is None:
                jacobian = self.minimizer_internal_function_by_order(1, allow_fd=True,
                                                                     **fopts)
            if hessian is None:
                hessian = self.minimizer_internal_function_by_order(2, **fopts)
            try:
                self.use_internals = use_internals
                coords = self.embed_coords(coords)
                if orthogonal_projection_generator is None:
                    orthogonal_projection_generator = self.get_internal_coordinate_constraints(coordinate_constraints)
                convergence, opt_coords, opt_settings = self.optimize_iterative(coords,
                                                                                func=func,
                                                                                jacobian=jacobian,
                                                                                hessian=hessian,
                                                                                logger=logger,
                                                                                coordinate_constraints=coordinate_constraints,
                                                                                orthogonal_projection_generator=orthogonal_projection_generator,
                                                                                gradient_modification_function=gradient_modification_function,
                                                                                convert_modification_distances=False,
                                                                                initialization_function=initialization_function,
                                                                                return_trajectory=return_trajectory,
                                                                                **opts
                                                                                )
                coords = self.unembed_coords(opt_coords)
                return convergence, coords, opt_settings
            finally:
                self.use_internals = use_internals

        else:
            return self.optimize_iterative(coords,
                                           logger=logger,
                                           coordinate_constraints=coordinate_constraints,
                                           orthogonal_projection_generator=orthogonal_projection_generator,
                                           gradient_modification_function=gradient_modification_function,
                                           initialization_function=initialization_function,
                                           return_trajectory=return_trajectory,
                                           mode=mode,
                                           **opts
                                           )

    @classmethod
    def get_relaxed_scan_options(cls):
        """
        **LLM Docstring**

        The full set of recognized keyword-argument names for `relaxed_scan`, extending `get_optimizer_options` with the scan-specific option names.

        :return: the tuple of recognized option names
        :rtype: tuple[str]
        """
        # import scipy.optimize._optimize
        return cls.get_optimizer_options() + (
            'absolute_mesh',
            'coordinate_expansion',
            'expansion_active_positions',
            'include_displacement_constraints',
            'shift',
            'adjust_displacements',
            'split_scan_mesh'
        )
    def relaxed_scan(self,
                     coords,
                     displacement_specs,
                     displacement_coords,
                     coordinate_constraints=None,
                     absolute_mesh=False,
                     orthogonal_projection_generator=None,
                     coordinate_expansion=None,
                     expansion_active_positions=None,
                     include_displacement_constraints=True,
                     shift=True,
                     adjust_displacements=True,
                     split_scan_mesh=True,
                     **optimization_settings
                     ):
        """
        **LLM Docstring**

        Perform a constrained ("relaxed") potential-energy scan: for each point on a mesh of displacement values along the specified coordinate(s), build a displaced starting structure (via a coordinate-index-based, coordinate-expansion-based, or user-supplied displacement-generator scheme) and re-optimize it subject to constraints holding the scanned coordinates fixed, splitting the scan mesh at the point nearest zero displacement so both directions can be scanned outward (zig-zag) from the reference geometry if `split_scan_mesh` is set.

        :param coords: the reference (starting) coordinates for the scan
        :type coords: np.ndarray
        :param displacement_specs: the `(start, stop, num)` range(s) (or explicit displacement values, if `absolute_mesh`) to scan each displaced coordinate over
        :type displacement_specs: Iterable
        :param displacement_coords: which coordinate(s) to displace along -- an internal-coordinate index/spec, a dict mapping coordinates to base displacement weights (defining a `coordinate_expansion` direction), or a callable/`(callable, callable)` pair generating custom displacement directions
        :type displacement_coords: object | dict | callable | tuple
        :param coordinate_constraints: additional coordinates to hold fixed during each relaxation, beyond the scanned ones
        :type coordinate_constraints: object | None
        :param absolute_mesh: whether `displacement_specs` gives explicit displacement values rather than `(start, stop, num)` ranges
        :type absolute_mesh: bool
        :param orthogonal_projection_generator: an explicit constraint projector generator, bypassing automatic construction
        :type orthogonal_projection_generator: callable | None
        :param coordinate_expansion: an explicit coordinate-transformation expansion direction to displace along
        :type coordinate_expansion: list[np.ndarray] | None
        :param expansion_active_positions: restrict `coordinate_expansion`'s output to these positions
        :type expansion_active_positions: Iterable[int] | None
        :param include_displacement_constraints: whether the scanned coordinates should also be added to `coordinate_constraints` (held fixed during relaxation)
        :type include_displacement_constraints: bool
        :param shift: whether displacement values are relative shifts or absolute targets
        :type shift: bool
        :param adjust_displacements: whether to convert the mesh of absolute displacement values into incremental steps between successive optimizations
        :type adjust_displacements: bool
        :param split_scan_mesh: whether to split the mesh at the point nearest zero and scan outward in both directions from the reference geometry
        :type split_scan_mesh: bool
        :param optimization_settings: extra options forwarded to `optimize` for each relaxation step
        :type optimization_settings: dict
        :return: `(did_optimize, geometries, metadata)` -- per-point convergence flags, the resulting geometries (ordered to match `displacement_specs`), and per-point optimizer metadata
        :rtype: tuple
        :raises NotImplementedError: if a Cartesian-only displacement generator is used without internal-coordinate handling, or if `coordinate_expansion` is combined with a custom displacement-generator pair
        """
        if isinstance(displacement_coords, dict):
            displacement_values = list(displacement_coords.values())
            displacement_coords = list(displacement_coords.keys())
            ndim, idx = self.get_internal_coordinate_indices(displacement_coords)
            expansion_vector = np.zeros(ndim)
            for i,v in zip(idx, displacement_values):
                expansion_vector[i] = v
            coordinate_expansion = [expansion_vector[np.newaxis]]
            expansion_active_positions = idx
        if not callable(displacement_coords) and callable(displacement_coords[0]):
            displacement_function, projected_coordinate_generator = displacement_coords
        else:
            if not callable(displacement_coords):
                if nput.is_int(displacement_coords[0]):
                    displacement_coords = [displacement_coords]
                if coordinate_constraints is None:
                    coordinate_constraints = displacement_coords
                elif include_displacement_constraints:
                    coordinate_constraints = list(coordinate_constraints) + list(displacement_coords)

                if self.use_internal_coordinate_handlers():
                    d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
                    evaluator = MolecularEvaluator(self.embedding, None)
                    if coordinate_expansion is not None:
                        idx = [0]
                    else:
                        ndim, idx = self.get_internal_coordinate_indices(displacement_coords)
                    def displacement_function(coords, disp, idx=idx):
                        """
                        **LLM Docstring**

                        Build a displaced starting structure for one scan point, given the current structure and a displacement value/tuple, using whichever displacement-generation scheme (index-based, coordinate-expansion, or user-supplied generator) was configured in the enclosing `relaxed_scan` call. Several closures with this same name are defined for the different displacement-generation branches.

                        :param coords: the current (pre-displacement) coordinates, in `self.distance_units`
                        :type coords: np.ndarray
                        :param disp: the displacement value(s) for this scan point
                        :type disp: object
                        :return: the displaced coordinates, in `self.distance_units`
                        :rtype: np.ndarray
                        """
                        return evaluator.get_displaced_coordinates(
                            [disp],
                            which=idx,
                            use_internals='reembed',
                            shift=shift,
                            coords=coords * d_conv,
                            strip_embedding=True,
                            coordinate_expansion=coordinate_expansion,
                            expansion_active_positions=expansion_active_positions
                        )[0] / d_conv
                else:
                    raise NotImplementedError("`displacement_coords` must be a displacement generator in Cartesians")
            else: # displacement generator
                if coordinate_expansion is not None:
                    raise NotImplementedError("using `displaced_coords` as a displacement generator with a `coordinate_expansion` not supported")

                d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
                evaluator = MolecularEvaluator(self.embedding, None)

                cur_disp = [None]

                def subprojector(coords, _):
                    """
                    **LLM Docstring**

                    Build the orthogonal-projection matrix removing the direction of the most recently applied custom displacement, used as (part of) the constraint projector when a user-supplied displacement generator is in use.

                    :param coords: the coordinates at the current scan point (unused directly, accepted for interface compatibility)
                    :type coords: np.ndarray
                    :param _: an unused mask argument, accepted for interface compatibility
                    :type _: object
                    :return: the projection matrix for the most recent displacement direction
                    :rtype: np.ndarray
                    """
                    cd: np.ndarray = cur_disp[0]
                    return nput.orthogonal_projection_matrix(cd)

                if orthogonal_projection_generator is None:
                    orthogonal_projection_generator = subprojector
                else:
                    def projector(coords, mask, base=orthogonal_projection_generator):
                        """
                        **LLM Docstring**

                        Combine a user-supplied orthogonal-projection generator with the enclosing `subprojector` (removing the most recent custom displacement direction), so both constraints apply together.

                        :param coords: the coordinates to compute the projector at
                        :type coords: np.ndarray
                        :param mask: the mask argument forwarded to the base projector generator
                        :type mask: object
                        :param base: the user-supplied projector generator, captured from the enclosing scope
                        :type base: callable
                        :return: the combined projection matrix
                        :rtype: np.ndarray
                        """
                        m = orthogonal_projection_generator(coords, mask)
                        m = m @ subprojector(coords, mask)
                        return m
                    orthogonal_projection_generator = projector
                if self.use_internal_coordinate_handlers():
                    def displacement_function(coords, disp):
                        """
                        **LLM Docstring**

                        Build a displaced starting structure for one scan point, given the current structure and a displacement value/tuple, using whichever displacement-generation scheme (index-based, coordinate-expansion, or user-supplied generator) was configured in the enclosing `relaxed_scan` call. Several closures with this same name are defined for the different displacement-generation branches.

                        :param coords: the current (pre-displacement) coordinates, in `self.distance_units`
                        :type coords: np.ndarray
                        :param disp: the displacement value(s) for this scan point
                        :type disp: object
                        :return: the displaced coordinates, in `self.distance_units`
                        :rtype: np.ndarray
                        """
                        old_coords = coords * d_conv
                        expansion, which, active_sites = displacement_coords(old_coords, disp)
                        new_coords = evaluator.get_displaced_coordinates(
                            [disp],
                            which=which,
                            use_internals='reembed',
                            shift=shift,
                            coords=coords * d_conv,
                            strip_embedding=True,
                            coordinate_expansion=expansion[np.newaxis],
                            expansion_active_positions=expansion_active_positions
                        )[0] / d_conv
                        cur_disp[0] = expansion
                        return new_coords
                else:
                    def displacement_function(coords, disp):
                        """
                        **LLM Docstring**

                        Build a displaced starting structure for one scan point, given the current structure and a displacement value/tuple, using whichever displacement-generation scheme (index-based, coordinate-expansion, or user-supplied generator) was configured in the enclosing `relaxed_scan` call. Several closures with this same name are defined for the different displacement-generation branches.

                        :param coords: the current (pre-displacement) coordinates, in `self.distance_units`
                        :type coords: np.ndarray
                        :param disp: the displacement value(s) for this scan point
                        :type disp: object
                        :return: the displaced coordinates, in `self.distance_units`
                        :rtype: np.ndarray
                        """
                        old_coords = coords * d_conv
                        expansion, which = displacement_coords(old_coords, disp)
                        new_coords = evaluator.get_displaced_coordinates(
                            [disp],
                            which=which,
                            use_internals=False,
                            shift=shift,
                            coords=old_coords,
                            coordinate_expansion=expansion,
                            expansion_active_positions=expansion_active_positions
                        )[0] / d_conv
                        cur_disp[0] = (new_coords - old_coords).reshape((-1, new_coords.shape[-2] * 3, 1))
                        return new_coords

        if nput.is_numeric(displacement_specs[0]):
            displacement_specs = [displacement_specs]

        if absolute_mesh:
            mesh = displacement_specs
        else:
            mesh = [
                np.linspace(*m)
                for m in displacement_specs
            ]

        if split_scan_mesh:
            # we find the nearest mesh point to (0,...) along the first axis and if need be scan both directions from there
            z_idx = np.argmin(np.abs(mesh[0]))
            if z_idx > 0:
                submeshes = [
                    ([np.flip(mesh[0][:z_idx])] + mesh[1:], True),
                    ([mesh[0][z_idx:]] + mesh[1:], False)
                ]
            else:
                submeshes = [(mesh, False)]
        else:
            submeshes = [(mesh, False)]

        full_did_opt = []
        full_geoms = []
        full_meta = []
        for mesh, flip in submeshes:
            if shift and adjust_displacements:
                _ = []
                for m in mesh:
                    _.append(np.concatenate([m[:1], np.diff(m)]))
                mesh = _

            opt_res = self._optimize_along_displacements(
                coords,
                mesh,
                displacement_function,
                coordinate_constraints=coordinate_constraints,
                orthogonal_projection_generator=orthogonal_projection_generator,
                **optimization_settings
            )
            did_opt = [is_opt for is_opt, geom, meta in opt_res]
            meta = [meta for is_opt, geom, meta in opt_res]
            geoms = [geom for is_opt, geom, meta in opt_res]
            if flip:
                full_did_opt.extend(reversed(did_opt))
                full_geoms.extend(reversed(geoms))
                full_meta.extend(reversed(meta))
            else:
                full_did_opt.extend(did_opt)
                full_geoms.extend(geoms)
                full_meta.extend(meta)
        return full_did_opt, np.asanyarray(full_geoms), full_meta

    def _optimize_along_displacements(self,
                                      coords,
                                      displacement_meshes,
                                      displacement_function,
                                      zigzag=True,
                                      mesh_iterator=None,
                                      **optimization_settings):
        """
        **LLM Docstring**

        Walk a mesh of displacement values (in zig-zag or plain product order), building each successive displaced structure from the previous optimized geometry (via `displacement_function`) and re-optimizing it, accumulating the sequence of optimization results.

        :param coords: the starting coordinates for the first displacement
        :type coords: np.ndarray
        :param displacement_meshes: the per-coordinate displacement value sequences to iterate over
        :type displacement_meshes: list[np.ndarray]
        :param displacement_function: builds a displaced structure from the current coordinates and a displacement tuple
        :type displacement_function: callable
        :param zigzag: whether to traverse the mesh in a zig-zag order (rather than a plain nested-loop product) for smoother successive displacements
        :type zigzag: bool
        :param mesh_iterator: an explicit iterator over displacement tuples, bypassing the zig-zag/product construction
        :type mesh_iterator: Iterator | None
        :param optimization_settings: extra options forwarded to `optimize` at each mesh point
        :type optimization_settings: dict
        :return: the list of `(converged, coords, info)` results, one per mesh point
        :rtype: list[tuple]
        """
        if mesh_iterator is None:
            if not zigzag:
                mesh_iterator = itertools.product(*displacement_meshes)
            else:
                mesh_iterator = itut.zigzag_product(*displacement_meshes)

        results = []
        for disp in mesh_iterator:
            new_struct = displacement_function(coords, disp)
            res = self.optimize(
                new_struct,
                **optimization_settings
            )
            results.append(res)
            coords = res[1]

        return results

@EnergyEvaluator.register('rdkit')
class RDKitEnergyEvaluator(EnergyEvaluator):
    def __init__(self, rdmol, force_field='mmff', charge=None, multiplicity=None, **defaults):
        super().__init__(**defaults)
        # if hasattr(mol, 'rdmol'):
        #     mol = mol.rdmol
        self.rdmol = rdmol
        self.force_field = force_field

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.rdmol, **cls.prep_mol_opts(mol, **opts))

    property_units = 'Kilocalories/Mole'
    analytic_derivative_order = 1
    def evaluate_term(self, coords, order, force_field_type=None, **opts):
        if force_field_type is None:
            force_field_type = self.force_field
        if order == 0:
            return self.rdmol.calculate_energy(coords, force_field_type=force_field_type, **(self.defaults | opts))
        elif order == 1:
            return self.rdmol.calculate_gradient(coords, force_field_type=force_field_type, **(self.defaults | opts))
        else:
            raise ValueError(f"order {order} not supported")

    def optimize(self,
                 coords,
                 method=None,
                 force_field_type=None,
                 unitary=False,
                 # generate_rotation=False,
                 # dtype='float64',
                 # orthogonal_directions=None,
                 # convergence_metric=None,
                 tol=1e-8,
                 max_iterations=100,
                 # damping_parameter=None,
                 # damping_exponent=None,
                 # reset_interval=None,
                 **opts
                 ):
        if force_field_type is None:
            force_field_type = self.force_field

        if method == 'rdkit':
            maxIters = opts.pop('maxIters', max_iterations)
            optimizer = opts.pop('optimizer', None)
            return self.rdmol.optimize_structure(
                geoms=coords,
                force_field_type=force_field_type,
                optimizer=optimizer,
                maxIters=maxIters,
                **opts
            )
        else:
            return super().optimize(
                coords,
                method=method,
                tol=tol,
                max_iterations=max_iterations,
                **opts
            )

@EnergyEvaluator.register('aimnet2')
class AIMNet2EnergyEvaluator(EnergyEvaluator):
    """
    Borrows structure from AIMNet2ASE to call appropriately
    """
    def __init__(self, atoms, model='aimnet2', charge=0, device=None, model_dir=None, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.eval = self.setup_aimnet(model, device=device, model_dir=model_dir)
        self.model = model
        self.atoms = atoms
        self.numbers = [AtomData[atom, "Number"] for atom in atoms]
        self.charge = charge
        self.multiplicity = multiplicity
        self.quiet = quiet

    @classmethod
    def handle_specialization(cls, tag):
        return {'model':tag}

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    def to_ase(self,
               gradient_modification_function=None,
               gradient_modification_mode='shift',
               convert_modification_distances=True,
               convert_modification_energies=True,
               orthogonal_projection_generator=None,
               **etc
               ):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with dev.OutputRedirect():
                try:
                    from aimnet.calculators import AIMNet2ASE
                except ImportError:
                    from aimnet2calc import AIMNet2ASE

        mult=self.multiplicity
        if mult is None:
            mult = 1
        calc = AIMNet2ASE(base_calc=self.eval, charge=self.charge, mult=mult)
        calc = self._modify_ase_calc(calc,
                                     gradient_modification_function=gradient_modification_function,
                                     gradient_modification_mode=gradient_modification_mode,
                                     convert_modification_distances=convert_modification_distances,
                                     convert_modification_energies=convert_modification_energies,
                                     orthogonal_projection_generator=orthogonal_projection_generator)
        return calc

    def to_pysis(self,
                 gradient_modification_function=None,
                 gradient_modification_mode='shift',
                 convert_modification_distances=False,
                 convert_modification_energies=False,
                 orthogonal_projection_generator=None,
                 **etc):
        from pysisyphus.config import OUT_DIR_DEFAULT

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with dev.OutputRedirect():
                try:
                    from aimnet.calculators import AIMNet2Pysis
                except ImportError:
                    from aimnet2calc import AIMNet2Pysis

        mult=self.multiplicity
        if mult is None:
            mult = 1
        base_calc = AIMNet2Pysis(model=self.eval, charge=self.charge, mult=mult,
                                 out_dir=OUT_DIR_DEFAULT,
                                 **etc)
        base_calc.copy = functools.partial(self._copy_pysis, base_calc)
        base_calc = self._modify_pysis_calc(base_calc,
                                            gradient_modification_function=gradient_modification_function,
                                            gradient_modification_mode=gradient_modification_mode,
                                            convert_modification_distances=convert_modification_distances,
                                            convert_modification_energies=convert_modification_energies,
                                            orthogonal_projection_generator=orthogonal_projection_generator)
        return base_calc

    @classmethod
    def _copy_pysis(cls, base_calc, **opts):
        return type(base_calc)(
            **( opts | dict(model=base_calc.model, charge=base_calc.charge, mult=base_calc.mult, out_dir=base_calc.out_dir))
        )

    def _maybe_download_asset(self, file:str, url:str):
        ...


    @classmethod
    @contextlib.contextmanager
    def _overload_aimnet_modeldir(cls, root_dir):
        if root_dir is None:
            root_dir = os.environ.get("MODEL_CACHE_DIR")

        if root_dir is not None:
            import aimnet.calculators.model_registry as reg
            cur = reg._maybe_download_asset
            try:
                def _maybe_download_asset(file: str, url: str) -> str:
                    import requests
                    filename = os.path.join(root_dir, 'aimnet', file)
                    if not os.path.exists(filename):
                        print(f"Downloading {url} -> {filename}")
                        with open(filename, "wb") as f:
                            response = requests.get(url, timeout=60)
                            f.write(response.content)
                    return filename
                reg._maybe_download_asset = _maybe_download_asset
                yield _maybe_download_asset
            finally:
                reg._maybe_download_asset = cur
        else:
            yield None


    @classmethod
    def setup_aimnet(cls, model, device=None, model_dir=None):
        if isinstance(model, str):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                with dev.OutputRedirect():
                    try:
                        with cls._overload_aimnet_modeldir(model_dir):
                            from aimnet.calculators import AIMNet2Calculator
                    except ImportError:
                        from aimnet2calc import AIMNet2Calculator

            with cls.quiet_mode():
                device = cls._resolve_torch_device(device)

                try:
                    calc = AIMNet2Calculator(model, device=device)
                except TypeError:
                    calc = AIMNet2Calculator(model) # can't set device in older aimnet

        else:
            calc = model

        return calc

    @staticmethod
    def autodiff(expr, coord, pad_dim, create_graph=False, retain_graph=False):
        import torch
        # here forces have shape (N, 3) and coord has shape (N+1, 3)
        # return hessian with shape (N, 3, N, 3)
        coord_ndim = len(coord.shape)
        shape = coord.shape + expr.shape
        grads = []
        # N = coord.shape[0] - 1
        # for p in itertools.combinations(coord.shape, expr.ndim):
        npad = len(pad_dim)
        rav_shape = tuple(
            np.prod([expr.shape[coord_ndim*i+k] for k in range(coord_ndim)], dtype=int)
            for i in range((expr.ndim - npad) // coord_ndim)
        ) + pad_dim
        if len(rav_shape) > 0:
            cache = {}
            inds = np.array(np.unravel_index(np.arange(np.prod(expr.shape, dtype=int)), rav_shape)).T
            if npad == 0:
                shinds = np.sort(inds, axis=1)
            else:
                shinds = np.concatenate([np.sort(inds[:, :-(npad)], axis=1), inds[:, -npad:]], axis=-1)
            for i,_f in enumerate(expr.flatten()):
                ravioli = tuple(inds[i])
                key = tuple(shinds[i])
                if ravioli == key:
                    g = torch.autograd.grad(_f,
                                            coord,
                                            # allow_unused=True,
                                            retain_graph=True,
                                            create_graph=create_graph)
                    g = g[0]
                    grads.append(g)
                    cache[key] = g
                else:
                    grads.append(cache[key])
        else:
            grads = []
            if len(expr.shape) == 0:
                grads.append(
                    torch.autograd.grad(expr,
                                        coord,
                                        # allow_unused=True,
                                        retain_graph=True,
                                        create_graph=create_graph)[0]
                )
            else:
                for n, _f in enumerate(expr):
                    g = torch.autograd.grad(_f,
                                            coord,
                                            # allow_unused=True,
                                            retain_graph=True,
                                            create_graph=create_graph)[0]
                    grads.append(g)

        deriv = torch.stack(grads).view(shape)
        shape_tuple = ((slice(None, -1),) + (slice(None),) * (coord_ndim-1)) * ((len(shape) - npad) // coord_ndim)
        return deriv, shape_tuple

    # def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
    #     super().calculate(atoms, properties, system_changes)
    #     self.update_tensors()
    #
    #     if self.atoms.cell is not None and self.atoms.pbc.any():
    #         # assert self.base_calc.cutoff_lr < float('inf'), 'Long-range cutoff must be finite for PBC'
    #         cell = self.atoms.cell.array
    #     else:
    #         cell = None
    #
    #     results = self.base_calc({
    #         'coord': torch.tensor(self.atoms.positions, dtype=torch.float32, device=self.base_calc.device),
    #         'numbers': self._t_numbers,
    #         'cell': cell,
    #         'mol_idx': self._t_mol_idx,
    #         'charge': self._t_charge,
    #         'mult': self._t_mult,
    #     }, forces='forces' in properties, stress='stress' in properties)

    property_units = 'ElectronVolts'
    batched_orders = True
    analytic_derivative_order = 2
    def prep_eval(self, coords, order, keep_graph=None, forces=None, hessian=None, **opts):
        import torch

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])
        # numbers = np.broadcast_to(
        #     np.array(self.numbers)[np.newaxis],
        #     (len(coords), len(self.numbers))
        # )

        if forces is None:
            forces = order > 0
        if hessian is None:
            hessian = order > 1

        arg_dict = {
            'coord': torch.tensor(coords.reshape(-1, 3), dtype=torch.float32, device=self.eval.device),
            'numbers': torch.tensor(
                np.repeat(np.array([self.numbers]), coords.shape[0], axis=0).flatten(),
                dtype=torch.int64, device=self.eval.device
            ),
            'charge': torch.tensor(self.charge, dtype=torch.float32, device=self.eval.device),
            'cell': None,
            'mol_idx': torch.tensor(
                np.repeat(np.arange(coords.shape[0]), coords.shape[1]).flatten(),
                dtype=torch.int64,
                device=self.eval.device
            )
        }

        if self.multiplicity is not None:
            arg_dict['mult'] = torch.tensor(self.multiplicity, dtype=torch.float32, device=self.eval.device)

        data = self.eval.prepare_input(arg_dict)
        if hessian and data['mol_idx'][-1] > 0:
            raise NotImplementedError('higher-derivative calculation is not supported for multiple molecules')
        data = self.eval.set_grad_tensors(data,
                                          forces=keep_graph or forces,
                                          stress=False,
                                          hessian=keep_graph or hessian
                                          )
        if keep_graph and not (forces or hessian):
            data['coord'].requires_grad_(True)
            self.eval._saved_for_grad['coord'] = data['coord']
        with torch.jit.optimized_execution(False):
            data = self.eval.model(data)
        data = self.eval.get_derivatives(
            data,
            forces=forces,
            hessian=hessian,
            **opts
        )

        return base_shape, coords, data

    def _apply_autodiff(self, data, root_key, deriv_key, order, property_shape,
                        diff_var=None,
                        key_tag='deriv',
                        pad_coord=1,
                        create_graph=None,
                        retain_graph=None):
        if pad_coord != 1:
            raise NotImplementedError("all AIMNet coordinates are padded by 1, need to relax this condition")
        keys = []
        if order > 0:
            if diff_var is None:
                diff_var = self.eval._saved_for_grad['coord']
            if deriv_key is not None:
                keys.append(f'{key_tag}_1')
                data[f'{key_tag}_1'] = (-data[deriv_key], (slice(None), slice(None)))
                for o in range(2, order + 1):
                    keys.append(f'{key_tag}_{o}')
                    data[f'{key_tag}_{o}'] = self.autodiff(data[f'{key_tag}_{o - 1}'][0],
                                                       diff_var,
                                                       property_shape,
                                                       create_graph=create_graph if create_graph is not None else ((o==1) or (o < order)),
                                                       retain_graph=retain_graph if retain_graph is not None else (o < order)
                                                       )
            else:
                if 'forces' in data: data['forces'].unbind()
                data[f'{key_tag}_0'] = (data[root_key], (slice(None),))
                for o in range(1, order + 1):
                    keys.append(f'{key_tag}_{o}')
                    data[f'{key_tag}_{o}'] = self.autodiff(data[f'{key_tag}_{o - 1}'][0],
                                                       diff_var,
                                                       property_shape,
                                                       create_graph=create_graph if create_graph is not None else ((o==1) or (o < order)),
                                                       retain_graph=retain_graph if retain_graph is not None else (o < order)
                                                       )
        return keys
    def _replace_nans(self, e, replacement):
        if replacement is True: replacement = 1e8
        if replacement is False or not nput.is_numeric(replacement):
            return e
        else:
            e = np.asanyarray(e)
            if e.ndim == 0:
                if np.isnan(e):
                    return replacement
                else:
                    return e
            else:
                e[np.isnan(e)] = replacement
            return e
    def process_aimnet_derivs(self, base_shape, coords, data, root_key, order,
                              extra_deriv_coord=None,
                              deriv_key=None,
                              property_shape=(),
                              pad_coord=1,
                              replace_nans=True,
                              output_shape=None):
        if extra_deriv_coord is not None:
            deriv_keys = []
            root_keys = []
            if deriv_key is not None: raise ValueError("mixed derivatives with a deriv_key not supported")
            var, suborder = extra_deriv_coord
            self._apply_autodiff(data, root_key, deriv_key, suborder, property_shape,
                                 diff_var=var,
                                 pad_coord=pad_coord,
                                 retain_graph=True,
                                 create_graph=True
                                 )
            for n in range(suborder+1):
                deriv, subshape = data[f"deriv_{n}"]
                if len(subshape) > 0 and len(deriv.shape) > 0:
                    deriv = deriv.__getitem__(subshape)
                data[f"subderiv_{n}"] = deriv
                dkeys = self._apply_autodiff(data, f"subderiv_{n}", None, suborder,
                                             data[f"subderiv_{n}"].shape,
                                             key_tag=f"deriv_{n}",
                                             pad_coord=pad_coord
                                             )
                deriv_keys += dkeys
                root_keys.append(f"subderiv_{n}")
        else:
            deriv_keys = self._apply_autodiff(data, root_key, deriv_key, order, property_shape,
                                              diff_var=None, pad_coord=pad_coord)
            root_keys = [root_key]

        if output_shape is None:
            output_shape = property_shape
        else:
            output_shape = tuple(output_shape)
        slice_shape = tuple(slice(None, o) for o in output_shape)
        for key in deriv_keys:
            deriv,subshape = data[key]
            if len(subshape + slice_shape) > 0:
                deriv = deriv.__getitem__(subshape + slice_shape)
            data[key] = deriv

        ko = self.eval.keys_out
        afk = self.eval.atom_feature_keys
        try:
            self.eval.keys_out = self.eval.keys_out + root_keys + deriv_keys
            if root_key not in self.eval.keys_out:
                self.eval.keys_out.append(root_key)
            self.eval.atom_feature_keys = self.eval.atom_feature_keys + deriv_keys
            expansion = self.eval.process_output(data)
        finally:
            self.eval.keys_out = ko
            self.eval.atom_feature_keys = afk

        ndim = np.prod(coords.shape[-2:], dtype=int)

        if extra_deriv_coord is None:
            expansion = [
                self._replace_nans(
                    expansion[k].detach().cpu().numpy().reshape(
                        base_shape + (ndim,) * o + output_shape
                    ),
                    replace_nans
                )
                for o,k in enumerate(root_keys + deriv_keys)
            ]
        else:
            var, nsub = extra_deriv_coord
            targ_shape = list(var.shape[len(base_shape):])
            targ_shape[0] -= pad_coord
            vdim = np.prod(targ_shape, dtype=int)
            _ = []
            for n in range(extra_deriv_coord[1]+1):
                subexpansion = []
                for i in range(order + 1):
                    esub = expansion[f"subderiv_{n}" if i == 0 else f"deriv_{n}_{i}"].detach().cpu().numpy()
                    e_new = esub.reshape(
                        base_shape + (ndim,) * i + (vdim,) * n + output_shape
                    )
                    e_new = self._replace_nans(e_new, replace_nans)
                    subexpansion.append(e_new)
                _.append(subexpansion)
            expansion = _

        return expansion
    def evaluate_term(self, coords, order, **opts):
        if coords.ndim == 2 or order < 2:
            base_shape, coords, data = self.prep_eval(coords, order, **opts)

            return self.process_aimnet_derivs(
                base_shape, coords, data,
                'energy',
                order,
                deriv_key='forces'
            )
        else:
            base_shape = coords.shape[:-2]
            coords = coords.reshape((-1,) + coords.shape[-2:])
            expansions = [
                self.evaluate_term(c, order, **opts)
                for c in coords
            ]
            return [
                np.array(e).reshape(base_shape + e[0].shape)
                for e in zip(*expansions)
            ]

    # def optimize(self,
    #              coords,
    #              method=None,
    #              mode=None,
    #              initialization_function=None,
    #              gradient_modification_function=None,
    #              gradient_modification_mode='shift',
    #              convert_modification_distances=True,
    #              convert_modification_energies=True,
    #              return_trajectory=False,
    #              **opts
    #              ):
    #     tol, max_iterations, max_displacement = [
    #         opts.pop(k, self.optimizer_defaults[k])
    #         for k in ["tol", "max_iterations", "max_displacement"]
    #     ]
    #     if mode == 'ase':
    #         from McUtils.ExternalPrograms import ASEMolecule
    #         # calc = self.to_ase()
    #
    #         calc = self.to_ase(
    #             gradient_modification_function=gradient_modification_function,
    #             gradient_modification_mode=gradient_modification_mode,
    #             convert_modification_distances = convert_modification_distances,
    #             convert_modification_energies = convert_modification_energies,
    #         )
    #
    #         if initialization_function is not None:
    #             d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
    #             coords = initialization_function(d_conv * coords)
    #             coords = coords / d_conv
    #
    #         if return_trajectory:
    #             traj = tempfile.NamedTemporaryFile().name
    #         else:
    #             traj = None
    #
    #         mol = ASEMolecule.from_coords(
    #             self.atoms,
    #             coords.reshape((-1,) + coords.shape[-2:])[0],
    #             calculator=calc
    #         )
    #         # if gradient_modification_function is not None:
    #         #     old_get_forces = calc.get_forces
    #         #     def prep_coords(atoms=None):
    #         #         if atoms is None: atoms = mol.mol
    #         #         return atoms.positions
    #         #     calc.get_forces = self._modify_gradient(
    #         #         calc.get_forces,
    #         #         gradient_modification_function,
    #         #         modification_mode=gradient_modification_mode,
    #         #         convert_modification_distances=convert_modification_distances,
    #         #         convert_modification_energies=convert_modification_energies,
    #         #         coord_prep=prep_coords
    #         #     )
    #         #
    #         # try:
    #         new_opt = mol.optimize_structure(coords,
    #                                          fmax=tol, steps=max_iterations,
    #                                          maxstep=max_displacement,
    #                                          trajectory=traj,
    #                                          method=method)
    #         # finally:
    #         #     if gradient_modification_function is not None:
    #         #         calc.get_forces = old_get_forces
    #         if return_trajectory:
    #             from ase.io.trajectory import Trajectory
    #             try:
    #                 with Trajectory(traj, mode='r') as reader:
    #                     traj_data = [
    #                         a.positions for a in reader
    #                     ]
    #             except FileNotFoundError:
    #                 ...
    #             else:
    #                 os.remove(traj)
    #             cond, coords, d = new_opt
    #             new_opt = (cond, (coords, traj_data), d)
    #         return new_opt
    #     else:
    #         return super().optimize(
    #             coords,
    #             method=method,
    #             tol=tol,
    #             mode=mode,
    #             max_iterations=max_iterations,
    #             max_displacement=max_displacement,
    #             return_trajectory=return_trajectory,
    #             initialization_function=initialization_function,
    #             gradient_modification_function=gradient_modification_function,
    #             gradient_modification_mode=gradient_modification_mode,
    #             convert_modification_distances=convert_modification_distances,
    #             convert_modification_energies=convert_modification_energies,
    #             **opts
    #         )

@EnergyEvaluator.register('ase')
class ASECalcEnergyEvaluator(EnergyEvaluator):
    def __init__(self, atoms, charge=0, multiplicity=1, quiet=True, embedding=None, **defaults):
        d = defaults.copy()
        d.pop('calc', None)
        super().__init__(**d)
        self.eval = self.setup_calc(**defaults)
        self.atoms = atoms
        self.numbers = [AtomData[atom, "Number"] for atom in atoms]
        self.charge = 0 if charge is None else charge
        self.multiplicity = 1 if multiplicity is None else multiplicity
        self.quiet = quiet
        self.embedding = embedding

    def to_ase(self,
               gradient_modification_function=None,
               gradient_modification_mode='shift',
               convert_modification_distances=True,
               convert_modification_energies=True,
               orthogonal_projection_generator=None,
               **etc
               ):
        calc = self._modify_ase_calc(self.eval,
                                     gradient_modification_function=gradient_modification_function,
                                     gradient_modification_mode=gradient_modification_mode,
                                     convert_modification_distances=convert_modification_distances,
                                     convert_modification_energies=convert_modification_energies,
                                     orthogonal_projection_generator=orthogonal_projection_generator)
        return calc

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    @classmethod
    def setup_calc(cls, *, calc, **settings):
        return calc

    property_units = 'ElectronVolts'
    batched_orders = True
    analytic_derivative_order = 2
    def prep_ase(self, coords):
        from McUtils.ExternalPrograms import ASEMolecule
        atoms = ASEMolecule.from_coords(
            self.atoms,
            coords.reshape((-1,) + coords.shape[-2:])[0],
            calculator=self.eval,
            charge=self.charge,
            spin=self.multiplicity
        )
        return atoms

    def evaluate_term(self, coords, order, **opts):
        coords = np.asanyarray(coords)
        n = len(self.atoms)
        if coords.shape[-1] == (n * 3):
            base_shape = coords.shape[:-1]
        else:
            base_shape = coords.shape[:-2]
        coords = coords.reshape((-1, n, 3))
        mol = self.prep_ase(coords)
        return [
            r.reshape(base_shape + r.shape[1:])
            for r in mol.calculate_energy(coords, order=order)
        ]

    # def optimize(self,
    #              coords,
    #              mode=None,
    #              initialization_function=None,
    #              gradient_modification_function=None,
    #              return_trajectory=False,
    #              **opts
    #              ):
    #     tol, max_iterations, max_displacement = [
    #         opts.pop(k, self.optimizer_defaults[k])
    #         for k in ["tol", "max_iterations", "max_displacement"]
    #     ]
    #     if mode == 'ase':
    #         coords = np.asanyarray(coords)
    #         mol = self.prep_ase(coords.reshape((-1,) + coords.shape[-2:])[0])
    #         if initialization_function is not None:
    #             d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
    #             coords = initialization_function(d_conv * coords)
    #             coords = coords / d_conv
    #
    #         if gradient_modification_function is not None:
    #             calc = mol.mol.calc
    #             import scipy.optimize as opt
    #             opt.minimize()
    #             old_get_forces = calc.get_forces
    #             d_conv = UnitsData.convert(self.distance_units, "BohrRadius")
    #             e_conv = UnitsData.convert(self.property_units, "Hartrees")
    #             def jacobian(atoms=None):
    #                 res = old_get_forces(atoms)
    #                 if atoms is None: atoms = mol.mol
    #                 crds = atoms.positions
    #                 return -gradient_modification_function(crds * d_conv, -res * e_conv) * d_conv / e_conv
    #             calc.get_forces = jacobian
    #
    #         if return_trajectory:
    #             traj = tempfile.NamedTemporaryFile().name
    #         else:
    #             traj = None
    #
    #         try:
    #             new_opt = mol.optimize_structure(coords,
    #                                              fmax=tol, steps=max_iterations,
    #                                              maxstep=max_displacement,
    #                                              trajectory=traj,
    #                                              method=method)
    #         finally:
    #             if gradient_modification_function is not None:
    #                 calc.get_forces = old_get_forces
    #         if return_trajectory:
    #             from ase.io.trajectory import Trajectory
    #             try:
    #                 with Trajectory(traj, mode='r') as reader:
    #                     traj_data = [
    #                         a.positions for a in reader
    #                     ]
    #             except FileNotFoundError:
    #                 ...
    #             else:
    #                 os.remove(traj)
    #             cond, coords, d = new_opt
    #             new_opt = (cond, (coords, traj_data), d)
    #         return new_opt
    #     else:
    #         return super().optimize(
    #             coords,
    #             mode=mode,
    #             tol=tol,
    #             max_iterations=max_iterations,
    #             max_displacement=max_displacement,
    #             return_trajectory=return_trajectory,
    #             initialization_function=initialization_function,
    #             gradient_modification_function=gradient_modification_function,
    #             **opts
    #         )

@EnergyEvaluator.register('mace')
class MACEEnergyEvaluator(ASECalcEnergyEvaluator):
    analytic_derivative_order = 2
    @classmethod
    def handle_specialization(cls, tag):
        return {'model':tag}

    @classmethod
    @contextlib.contextmanager
    def _overload_mace_modeldir(cls, root_dir):
        cur_xdg = os.environ.get("XDG_CACHE_HOME")
        if root_dir is None:
            root_dir = os.environ.get("MODEL_CACHE_DIR")

        if root_dir is not None:
            try:
                os.environ["XDG_CACHE_HOME"] = root_dir
                yield root_dir
            finally:
                if cur_xdg is not None:
                    os.environ["XDG_CACHE_HOME"] = cur_xdg
                else:
                    os.environ.pop("XDG_CACHE_HOME")
        else:
            yield None

    model_types = {
        'extra_large':'omol/extra_large',
        "small":'off/small',
        "medium":'off/medium',
        "large":'off/large',
        'anicc':'anicc/none',
        'off':'off/large',
        'mp':'mp/medium-mpa-0',
        'polar':'polar/polar-1-l'
    }
    @classmethod
    def setup_calc(cls, model='extra_large', model_type=None, device=None, model_dir=None, **settings):
        model = model.lower()
        model = cls.model_types.get(model, model)
        if model_type is None:
            model_type, model = model.split('/', 1)
        if isinstance(model_type, str):
            with cls.quiet_mode():
                if model_type == 'omol':
                    from mace.calculators import mace_omol as model_type
                elif model_type == 'mp':
                    from mace.calculators import mace_mp as model_type
                elif model_type == 'anicc':
                    from mace.calculators import mace_anicc as model_type
                elif model_type == 'polar':
                    from mace.calculators import mace_anicc as model_type
                elif model_type == 'off':
                    from mace.calculators import mace_off as model_type
                else:
                    import mace.calculators
                    model_type = getattr(mace.calculators, model_type)

        with cls.quiet_mode(), cls._overload_mace_modeldir(model_dir):
            device = cls._resolve_torch_device(device)
            if model == 'none':
                calc = model_type(model=model, device=device, **settings)
            else:
                calc = model_type(model=model, device=device, **settings)

        return calc

@EnergyEvaluator.register('uma')
class UMAEnergyEvaluator(ASECalcEnergyEvaluator):
    @classmethod
    def handle_specialization(cls, tag):
        return {'model':tag}

    model_types = {
        'medium':'uma-m-1p1/omol',
        'small':'uma-s-1p2/omol',
        'crystals':'uma-s-1p2/omc',
        'crystals/small':'uma-s-1p2/omc',
        'crystals/medium':'uma-m-1p1/omc',
        'inorganics/small':'uma-s-1p2/omat',
        'inorganics/medium':'uma-m-1p1/omat',
    }
    @classmethod
    def setup_calc(cls, model="uma-s-1p2",
                   device=None, task_name=None, model_dir=None,
                   login=None,
                   hf_token=None,
                   predict_unit_options=None, **settings):
        from fairchem.core import pretrained_mlip, FAIRChemCalculator

        if login is None:
            if hf_token is None:
                hf_token = os.environ.get("HF_TOKEN")
            login = hf_token is not None
        if login:
            if hf_token is None:
                hf_token = os.environ.get("HF_TOKEN")
            from huggingface_hub import login
            login(token=hf_token)

        model = cls.model_types.get(model, model)
        if task_name is None:
            if '/' in model:
                model, task_name = model.split('/', 1)
            else:
                task_name = 'omol'

        if predict_unit_options is None:
            predict_unit_options = {}
        predict_unit_options = predict_unit_options.copy()
        if model_dir is None:
            model_dir = os.environ.get("MODEL_CACHE_DIR")
        if model_dir is not None and 'cache_dir' not in predict_unit_options:
            predict_unit_options['cache_dir'] = os.path.join(model_dir, 'fairchem')

        model = model.lower()
        with cls.quiet_mode():
            if 'device' not in predict_unit_options:
                device = cls._resolve_torch_device(device)
                predict_unit_options['device'] = device
            predictor = pretrained_mlip.get_predict_unit(model, **predict_unit_options)
            calc = FAIRChemCalculator(predictor, task_name=task_name, **settings)

        return calc

@EnergyEvaluator.register('upet')
class UPETEnergyEvaluator(ASECalcEnergyEvaluator):
    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    model_types = {
        # --- Current / recommended ---
        "mad": "pet-mad-s",
        "mad-xs": "pet-mad-xs",

        "oam": "pet-oam-xl",
        "oam-l": "pet-oam-l",

        "omat-xs": "pet-omat-xs",
        "omat-s": "pet-omat-s",
        "omat-m": "pet-omat-m",
        "omat-l": "pet-omat-l",
        "omat-xl": "pet-omat-xl",

        "omatpes": "pet-omatpes-l",

        "spice": "pet-spice-s",
        "spice-l": "pet-spice-l",

        # --- Legacy ---
        "omad-xs": "pet-omad-xs",
        "omad-s": "pet-omad-s",
        "omad-l": "pet-omad-l",
    }

    @classmethod
    def setup_calc(cls, model="pet-mad-s", device=None, version=None, login=None,
                   hf_token=None,
                   model_dir=None,
                   **settings):
        from upet.calculator import UPETCalculator

        if model_dir is None:
            model_dir = os.environ.get("MODEL_CACHE_DIR")

        model = cls.model_types.get(model, model)
        if version is None:
            if '/' in model:
                model, version = model.split('/', 1)
            else:
                version = 'lastest'

        if model_dir is not None:
            cur_cache = os.environ.get('HF_HUB_CACHE')
            os.environ['HF_HUB_CACHE'] = model_dir
        else:
            cur_cache = False

        try:
            if login is None:
                if hf_token is None:
                    hf_token = os.environ.get("HF_TOKEN")
                login = hf_token is not None
            if login:
                if hf_token is None:
                    hf_token = os.environ.get("HF_TOKEN")
                from huggingface_hub import login
                login(token=hf_token)

            model = model.lower()
            with cls.quiet_mode():
                if device is None:
                    device = cls._resolve_torch_device(device)
                calc = UPETCalculator(model=model, version=version, device=device, **settings)
        finally:
            if cur_cache is not False:
                if cur_cache is None:
                    del os.environ['HF_HUB_CACHE']
                else:
                    os.environ['HF_HUB_CACHE'] = cur_cache

        return calc

@EnergyEvaluator.register('orb')
class OrbModelEnergyEvaluator(ASECalcEnergyEvaluator):
    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    model_types = {
        'v1':"orbmol-v1-conservative",
        'v2':"orbmol-v2",
        'omol':"orbmol-v3-conservative-omol",
        'v3':"orbmol-v3-conservative-omol",
        'v3-direct':"orbmol-v3-direct-omol",
        'v1-direct':"orbmol-v1-direct",
        'omat':"orb-v3-conservative-inf-omat",
        'v3-omat':"orb-v3-conservative-inf-omat",
        'v3-omat-direct':"orb-v3-direct-inf-omat",
        'v3-omat20':"orb-v3-conservative-20-omat",
        'v3-omat20-direct':"orb-v3-direct-20-omat",
        'd3':"separate-d3-5layer",
        'd3-3layer':"separate-d3-3layer",
        'd4':"separate-d3-5layer",
        'd4-3layer':"separate-d4-3layer"
    }
    @classmethod
    def setup_calc(cls, model='omol', device=None, version=None, login=None, **model_opts):
        from orb_models.forcefield import pretrained, ORB_PRETRAINED_MODELS
        from orb_models.forcefield.inference.calculator import ORBCalculator

        model = cls.model_types.get(model, model)
        model = ORB_PRETRAINED_MODELS.get(cls.model_types.get(model, model), model)

        device = "cpu"  # or device="cuda"
        orbff, atoms_adapter = ORB_PRETRAINED_MODELS[model](
            device=device,
            **model_opts,
        )

        calc = ORBCalculator(orbff, atoms_adapter=atoms_adapter, device=device)
        return calc

@EnergyEvaluator.register('hipnn')
class HIPNNEnergyEvaluator(EnergyEvaluator):
    property_key = 'energies'
    gradient_key = 'Grad'
    def __init__(self, atoms, model_dir,
                 property_key=dev.default,
                 gradient_key=dev.default,
                 charge=0, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.eval = self.setup_model(model_dir)
        self.atoms = atoms
        self.numbers = [AtomData[atom, "Number"] for atom in atoms]
        self.charge = charge
        self.multiplicity = multiplicity
        self.quiet = quiet
        if dev.is_default(property_key, allow_None=False):
            property_key = self.property_key
        self.property_key = property_key
        if dev.is_default(gradient_key, allow_None=False):
            gradient_key = self.gradient_key
        self.gradient_key = gradient_key

    @classmethod
    def handle_specialization(cls, tag):
        return {'model_dir':tag}

    @classmethod
    def from_mol(cls, mol, *, model_dir, **opts):
        return cls(mol.atoms, model_dir=model_dir, **cls.prep_mol_opts(mol, **opts))

    @staticmethod
    def autodiff(expr, coord_list, pad_dim, create_graph=False, retain_graph=False):
        import torch
        # here forces have shape (N, 3) and coord has shape (N+1, 3)
        # return hessian with shape (N, 3, N, 3)
        shape = coord_list.shape + expr.shape[1:]
        grads = []
        # N = coord.shape[0] - 1
        # for p in itertools.combinations(coord.shape, expr.ndim):
        npad = len(pad_dim)
        rav_shape = tuple(
            (expr.shape[1+2*i]*expr.shape[2+2*i])
            for i in range((expr.ndim - 1 - npad) // 2)
        )
        cache = {}
        if len(rav_shape) > 0:
            inds = np.array(np.unravel_index(np.arange(np.prod(expr.shape[1:-npad], dtype=int)), rav_shape)).T
            shinds = np.sort(inds, axis=1)
            for n,subexpr in enumerate(expr):
                for i,_f in enumerate(subexpr.flatten()):
                    ravioli = tuple(inds[i])
                    key = tuple(shinds[i])
                    if ravioli == key:
                        g = torch.autograd.grad(_f,
                                                coord_list,
                                                # allow_unused=True,
                                                retain_graph=True,
                                                create_graph=create_graph)
                        g = g[0][n]
                        grads.append(g)
                        cache[key] = g
                    else:
                        grads.append(cache[key])
        else:
            for n,_f in enumerate(expr):
                g = torch.autograd.grad(_f,
                                        coord_list,
                                        # allow_unused=True,
                                        retain_graph=True,
                                        create_graph=create_graph)[0]
                grads.append(g[n])

        deriv = torch.stack(grads).view(shape)
        shape_tuple = (slice(None),) #+ (slice(None, None), slice(None)) * (len(expr.shape - 1) // 2)
        return deriv, shape_tuple

    def process_derivatives(self, expr, coords, order):
        terms = [(expr, slice(None,))]
        pad_dim = expr.shape[1:]
        for o in range(1, order + 1):
            terms.append(
                self.autodiff(terms[-1][0],
                              coords,
                              pad_dim,
                              # create_graph=(o==1),
                              create_graph=(o < order),
                              retain_graph=(o < order)
                              )
            )
        return [
            deriv.__getitem__(subshape)
            for deriv, subshape in terms
        ]

    def setup_model(cls, model_dir, *, device=None, **opts):
        import torch
        import hippynn.graphs
        import hippynn.experiment.serialization
        cur_dir = os.getcwd()

        device = cls._resolve_torch_device(device)
        try:
            os.chdir(model_dir)
            with torch.no_grad():
                model = hippynn.experiment.serialization.load_model_from_cwd(model_device=device)
        finally:
            os.chdir(cur_dir)

        predictor = hippynn.graphs.Predictor.from_graph(model)
        return predictor

    def print_model_properties(self):
        self.eval.graph.print_structure()

    property_units = 'Kilocalories/Mole'
    batched_orders = True
    analytic_derivative_order = 1
    def process_results(self, base_shape, coords, outputs, order,
                        property_key=None,
                        gradient_key=None,
                        property_shape=(),
                        output_shape=None,
                        **etc):
        if property_key is None:
            property_key = self.property_key
        if gradient_key is None:
            gradient_key = self.gradient_key
        property_shape = tuple(property_shape)
        terms = [outputs[property_key]]
        if order > 0:
            if gradient_key is not None:
                g = outputs[gradient_key]
                terms.append(g)
                order = order - 1
            if order > 0:
                terms.extend(
                    self.process_derivatives(
                        terms[-1],
                        coords,
                        order
                    )[1:]
                )

        n = coords.shape[-1] * coords.shape[-2]
        if output_shape is None:
            output_shape = property_shape
        else:
            output_shape = tuple(output_shape)
        return [
            e.detach().numpy().reshape(base_shape + (n,)*i + output_shape)
            for i,e in enumerate(terms)
        ]

        g.reshape(base_shape + (g.shape[-(1 + nshape)] * g.shape[-(2 + nshape)],) + property_shape)
        # # total molecular energy [kcal/mol]
        # outputs['Grad']  # forces [kcal/mol/Å]
        # outputs['dipole']  # dipole moment [e∙Å]
    def prep_eval(self, coords, order, **opts):
        import torch

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])
        # numbers = np.broadcast_to(
        #     np.array(self.numbers)[np.newaxis],
        #     (len(coords), len(self.numbers))
        # )

        torch_coords = torch.tensor(coords, dtype=torch.float32, device=self.eval.model_device)
        arg_dict = {
            'coordinates': torch_coords,
            'species': torch.tensor(
                np.repeat(np.array([self.numbers]), coords.shape[0], axis=0),
                dtype=torch.int64, device=self.eval.model_device
            )
            # 'charge': torch.tensor(self.charge, dtype=torch.float32, device=self.eval.device),
            # 'cell': None,
            # 'mol_idx': torch.tensor(
            #     np.repeat(np.arange(coords.shape[0]), coords.shape[1]).flatten(),
            #     dtype=torch.int64,
            #     device=self.eval.device
            # )
        }
        # if self.multiplicity is not None:
        #     arg_dict['mult'] = torch.tensor(self.multiplicity, dtype=torch.float32, device=self.eval.device)
        # data = self.eval.prepare_input(arg_dict)
        # if hessian and data['mol_idx'][-1] > 0:
        #     raise NotImplementedError('higher-derivative calculation is not supported for multiple molecules')
        # data = self.eval.set_grad_tensors(data,
        #                                   forces=forces,
        #                                   stress=False,
        #                                   hessian=hessian
        #                                   )

        req_grad = self.eval.requires_grad
        prop_key = self.property_key
        grad_key = self.gradient_key
        outputs = None
        output_names = None
        output_dbnames = None
        try:
            self.eval.requires_grad = order > 0
            if prop_key is not None and prop_key not in self.eval.out_names:
                prop = self.eval.graph.node_from_name(prop_key)
                prop.bname = prop.db_name
                if outputs is None:
                    outputs = list(self.eval.outputs)
                    output_names = list(self.eval.out_names)
                    output_dbnames = list(self.eval.out_dbnames)
                self.eval.add_output(prop)
            if grad_key is not None and grad_key not in self.eval.out_names:
                prop = self.eval.graph.node_from_name(grad_key)
                prop.bname = prop.db_name
                if outputs is None:
                    outputs = list(self.eval.outputs)
                    output_names = list(self.eval.out_names)
                    output_dbnames = list(self.eval.out_dbnames)
                self.eval.add_output(grad_key)

            if self.quiet:
                with self.quiet_mode():
                    with torch.jit.optimized_execution(False):
                        data = self.eval(**arg_dict, **opts)
            else:
                with torch.jit.optimized_execution(False):
                    data = self.eval(**arg_dict, **opts)
        finally:
            self.eval.requires_grad = req_grad
            if outputs is not None:
                self.eval.outputs[:] = outputs
                self.eval.out_names[:] = output_names
                self.eval.out_dbnames[:] = output_dbnames

        return base_shape, torch_coords, data
    def evaluate_term(self, coords, order, **opts):
        if coords.ndim == 2 or order < 2:
            base_shape, coords, data = self.prep_eval(coords, order, **opts)
            data = self.process_results(
                base_shape,
                coords,
                data,
                order,
                property_shape=()
                # hessian=hessian,
                # **opts
            )
            return data
        else:
            base_shape = coords.shape[:-2]
            coords = coords.reshape((-1,) + coords.shape[-2:])
            expansions = [
                self.evaluate_term(c, order, **opts)
                for c in coords
            ]
            return [
                np.array(e).reshape(base_shape + e[0].shape)
                for e in zip(*expansions)
            ]


@EnergyEvaluator.register('xtb')
class XTBEnergyEvaluator(EnergyEvaluator):
    """
    Uses XTB to calculate
    """
    def __init__(self, atoms, method="GFN2-xTB", charge=0, quiet=True, **defaults):
        super().__init__(**defaults)

        # self.eval = self.setup_aimnet(model)
        # self.model = model
        self.atoms = atoms
        self.numbers = np.array([AtomData[atom, "Number"] for atom in atoms])
        # self.charge = charge
        # self.multiplicity = multiplicity
        self.quiet = quiet
        self.method = method
        self.charge = charge

    @classmethod
    def handle_specialization(cls, tag):
        return {'method':tag}

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    def get_single_point(self, coords):
        from xtb.interface import Calculator
        from xtb.utils import get_method
        calc = Calculator(get_method(self.method), self.numbers, coords,
                          charge=self.charge)
        if self.quiet:
            calc.set_verbosity('muted')
        return calc.singlepoint()


    distance_units = "BohrRadius"
    property_units = 'Hartrees'
    batched_orders = True
    analytic_derivative_order = 1
    def evaluate_term(self, coords, order, **opts):
        # import torch

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])

        if order > 1:
            raise NotImplementedError('`xtb` only implements up through gradients')

        energies = np.empty(coords.shape[0], dtype=float)
        if order > 0:
            grads = np.empty(coords.shape, dtype=float)
        else:
            grads = None
        for i,coord in enumerate(coords):
            res = self.get_single_point(coord)
            energies[i] = res.get_energy()
            if order > 0:
                grads[i] = res.get_gradient()

        energies = energies.reshape(base_shape)
        if order > 0:
            grads = grads.reshape(base_shape +(-1,))
            return [energies, grads]
        else:
            return [energies]

        # if order > 2:
        # if order > (forces=False, stress=False, hessian=False)
        # if order == 0:
        #     return self.rdmol.calculate_energy(coords, **opts)
        # elif order == 2:
        #     return self.rdmol.calculate_gradient(coords, **opts)
        # else:
        #     raise ValueError(f"order {order} not supported")

        # @staticmethod
        # def calculate_hessian(forces, coord):
        #     # here forces have shape (N, 3) and coord has shape (N+1, 3)
        #     # return hessian with shape (N, 3, N, 3)
        #     hessian = - torch.stack([
        #         torch.autograd.grad(_f, coord, retain_graph=True)[0]
        #         for _f in forces.flatten().unbind()
        #     ]).view(-1, 3, coord.shape[0], 3)[:-1, :, :-1, :]
        #     return hessian

    # def optimize(self,
    #              coords,
    #              method='quasi-newton',
    #              tol=1e-8,
    #              max_iterations=25,
    #              damping_parameter=None,
    #              damping_exponent=None,
    #              restart_interval=None,
    #              max_displacement=None,
    #              **opts
    #              ):
    #     if method == 'xtb-base':
    #         ...
    #     else:
    #         return super().optimize(
    #             coords,
    #             method=method,
    #             tol=tol,
    #             max_iterations=max_iterations,
    #             damping_parameter=damping_parameter,
    #             damping_exponent=damping_exponent,
    #             restart_interval=restart_interval,
    #             **opts
    #         )

@EnergyEvaluator.register('pyscf')
class PySCFEnergyEvaluator(EnergyEvaluator):
    """
    """
    def __init__(self,
                 atoms, *,
                 level_of_theory,
                 basis_set,
                 caller=None,
                 charge=None,
                 multiplicity=None,
                 quiet=True,
                 molecule_options=None,
                 level_of_theory_options=None,
                 disp=None,
                 **defaults):
        super().__init__(**defaults)

        # self.eval = self.setup_aimnet(model)
        # self.model = model
        self.atoms = atoms
        self.basis_set = basis_set
        self.level_of_theory = level_of_theory
        if level_of_theory_options is None:
            level_of_theory_options = {}
        self.level_of_theory_options = level_of_theory_options.copy()
        if disp is not None:
            self.level_of_theory_options['disp'] = disp
        self._pyscf_caller_dispatch = dev.uninitialized
        self.caller = self.resolve_caller(caller, level_of_theory)
        # self.charge = charge
        # self.multiplicity = multiplicity
        if molecule_options is None:
            molecule_options = {}
        if charge is not None:
            molecule_options['charge'] = molecule_options.get('charge', charge)
        if multiplicity is not None:
            molecule_options['spin'] = molecule_options.get('spin', multiplicity)
        if quiet:
            molecule_options['verbose'] = molecule_options.get('verbose', 0)
        self.molecule_options = molecule_options

    @classmethod
    def handle_specialization(cls, tag):
        lot, basis = tag.split("/", 1)
        return {'level_of_theory':lot, 'basis':basis}

    _default_caller = 'dft'
    @property
    def caller_dispatch(self) -> dev.OptionsMethodDispatch:
        self._pyscf_caller_dispatch = dev.handle_uninitialized(
            self._pyscf_caller_dispatch,
            dev.OptionsMethodDispatch,
            args=(self.get_callers,),
            kwargs=dict(
                default_method=self._default_caller,
                # attributes_map=cls.get_evaluators_by_attributes()
            )
        )
        return self._pyscf_caller_dispatch

    def get_callers(self):
        return {
            'mp2': self._call_mp2,
            'dft': self._call_dft,
            'hf': self._call_hf,
            'rhf': self._call_rhf,
            'uhf': self._call_uhf,
            # 'cc': self._call_cc
        }
    def resolve_caller(self, caller, lot):
        if caller is None:
            if callable(lot):
                caller = lot
            else:
                caller, opts = self.caller_dispatch.resolve(lot)
                self.level_of_theory_options.update(opts)
        return caller

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    def get_molecule(self, coords):
        from pyscf import gto

        molecule = gto.Mole(
            basis=self.basis_set,
            **self.molecule_options
        )
        molecule.atom = [
            (a, crd)
            for a,crd in zip(self.atoms, coords)
        ]
        molecule.build()

        return molecule

    def _call_dft(self, molecule, restricted=True, functional_options=None,  disp=None, **opts):
        from pyscf import dft

        if functional_options is None:
            functional_options = {}
        functional_options = functional_options.copy()
        functional_options['xc'] = self.level_of_theory
        d = functional_options.pop('disp', None)
        if disp is None:
            disp = d
        if restricted:
            calc = dft.RKS(molecule, **functional_options)
        else:
            calc = dft.UKS(molecule, **functional_options)
        if disp is not None:
            calc.disp = disp
        return calc.run(**opts)
        # calc = calc.newton()
        # calc.kernel()
        # return calc

    def _call_hf(self, molecule, restricted=False, **opts):
        from pyscf import scf
        if restricted:
            return scf.HF(molecule).run(**opts)
        else:
            return scf.HF(molecule).run(**opts)
    def _call_rhf(self, molecule, **opts):
        return self._call_hf(molecule, restricted=True, **opts)
    def _call_uhf(self, molecule, **opts):
        return self._call_hf(molecule, restricted=False, **opts)

    def _call_mp2(self, molecule, reference='hf', **etc):
        from pyscf import mp
        hf_caller, opts = self.caller_dispatch.resolve(reference)
        mf = hf_caller(molecule, **opts)
        return mp.MP2(mf, **etc).run()

    def get_energy_from_calc(self, calc) -> float:
        if hasattr(calc, 'tot_energy'):
            return calc.tot_energy()
        elif hasattr(calc, 'e_tot'):
            return calc.e_tot
        else:
            return calc.kernel()

    def get_gradient_from_calc(self, calc):
        return calc.nuc_grad_method().run()

    def run_calculation(self, coords, **opts):
        molecule = self.get_molecule(coords)
        return self.caller(molecule, **(self.level_of_theory_options | opts))

    distance_units = "BohrRadius"
    property_units = 'Hartrees'
    batched_orders = True
    analytic_derivative_order = 1
    def evaluate_term(self, coords, order, **opts):
        # import torch

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])

        if order > 1:
            raise NotImplementedError('`xtb` only implements up through gradients')

        energies = np.empty(coords.shape[0], dtype=float)
        if order > 0:
            grads = np.empty(coords.shape, dtype=float)
        else:
            grads = None
        for i,coord in enumerate(coords):
            res = self.run_calculation(coord, **opts)
            energies[i] = self.get_energy_from_calc(res)
            if order > 0:
                grads[i] = self.get_gradient_from_calc(res)

        energies = energies.reshape(base_shape)
        if order > 0:
            grads = grads.reshape(base_shape +(-1,))
            return [energies, grads]
        else:
            return [energies]

        # if order > 2:
        # if order > (forces=False, stress=False, hessian=False)
        # if order == 0:
        #     return self.rdmol.calculate_energy(coords, **opts)
        # elif order == 2:
        #     return self.rdmol.calculate_gradient(coords, **opts)
        # else:
        #     raise ValueError(f"order {order} not supported")

        # @staticmethod
        # def calculate_hessian(forces, coord):
        #     # here forces have shape (N, 3) and coord has shape (N+1, 3)
        #     # return hessian with shape (N, 3, N, 3)
        #     hessian = - torch.stack([
        #         torch.autograd.grad(_f, coord, retain_graph=True)[0]
        #         for _f in forces.flatten().unbind()
        #     ]).view(-1, 3, coord.shape[0], 3)[:-1, :, :-1, :]
        #     return hessian

    # def optimize(self,
    #              coords,
    #              method='quasi-newton',
    #              tol=1e-8,
    #              max_iterations=25,
    #              damping_parameter=None,
    #              damping_exponent=None,
    #              restart_interval=None,
    #              max_displacement=None,
    #              **opts
    #              ):
    #     if method == 'xtb-base':
    #         ...
    #     else:
    #         return super().optimize(
    #             coords,
    #             method=method,
    #             tol=tol,
    #             max_iterations=max_iterations,
    #             damping_parameter=damping_parameter,
    #             damping_exponent=damping_exponent,
    #             restart_interval=restart_interval,
    #             **opts
    #         )

# @EnergyEvaluator.register('job')
class EvaluationServerEnergyEvaluator(EnergyEvaluator):
    def __init__(self, *, request_handler, **defaults):
        super().__init__(**defaults)
        self.request_handler = request_handler

    @abc.abstractmethod
    def setup_request(self, coords, order=None, tasks='energy', **opts) -> (list[Any], dict[Any, Any], Any):
        ...

    @abc.abstractmethod
    def process_response(self, response, state):
        ...

    @abc.abstractmethod
    def cleanup_state(self, response, state):
        ...

    def evaluate_term(self, coords, order, tasks='energy', **opts):
        args, kwargs, state = self.setup_request(coords, tasks=tasks, order=order, **opts)
        response = None
        try:
            response = self.request_handler(*args, **kwargs)
            return self.process_response(response, state)
        finally:
            self.cleanup_state(response, state)

    def optimize(self,
                 coords,
                 method='server',
                 **opts
                 ):

        if method == 'server':
            args, kwargs, state = self.setup_request(coords, tasks='optimize', **opts)
            response = None
            try:
                response = self.request_handler(*args, **kwargs)
                return self.process_response(response, state)
            finally:
                self.cleanup_state(response, state)
        else:
            return super().optimize(
                coords,
                method=method,
                **opts
            )

@EnergyEvaluator.register('mlipserver')
class MLIPServerEnergyEvaluator(EvaluationServerEnergyEvaluator):
    distance_units = 'BohrRadius'
    def __init__(self, atoms, *, container_path, energy_evaluator,
                 session:subprocess.Popen=None,
                 port=None,
                 temp_dir=None,
                 launcher_options=None,
                 bind_sources=True,
                 container_env=None,
                 container_mode=None,
                 initialization_delay_time=1,
                 charge=None,
                 multiplicity=None,
                 model_dir=dev.default,
                 pass_tokens=True,
                 **defaults):
        self.atoms = atoms
        self.charge = charge
        self.multiplicity = multiplicity
        self.port = port
        self.container_path = container_path
        self._managed_session = False
        self._launcher:SingularityLauncher = None
        self.session = session
        self.evaluator = energy_evaluator
        self.temp_dir = None
        self._tmp = temp_dir
        self.launcher_options = launcher_options
        self.bind_sources = bind_sources
        self.container_env = container_env
        self.initialization_delay_time = initialization_delay_time
        if container_mode is None:
            container_mode = 'exec' if self.container_env is not None else 'run'
        self.container_mode = container_mode
        if dev.is_default(model_dir, allow_None=False):
            model_dir = os.path.join(os.path.expanduser("~"), ".cache")
        self.model_dir = model_dir
        self.pass_tokens = pass_tokens
        super().__init__(request_handler=self._run_mlip_request, **defaults)

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    @classmethod
    def pick_port(cls):
        import socket
        with socket.socket() as s:
            s.bind(("", 0))
            return s.getsockname()[1]

    def get_container_env_command(self):
        if self.container_env is None:
            return []
        env = self.container_env
        if isinstance(self.container_env, str):
            env = ['conda', 'run', '--no-capture-output', '-n', self.container_env, 'python', '-u', '-m', 'mlipenv']
        return env

    #bind sources to make a container runnable
    dynamic_sources = ['McUtils', 'Psience', 'mlipenv']
    token_list = ['HF_TOKEN']
    def get_launcher(self):
        lopts = self.launcher_options
        if lopts is None:
            lopts = {}

        bind_paths = lopts.pop('bind', {})
        if isinstance(bind_paths, dict):
            bind_paths[self.temp_dir] = '/tmp/io'
        else:
            bind_paths = list(bind_paths) + [self.temp_dir + ':' + "/tmp/io"]
        if self.bind_sources:
            cur_sources = lopts.pop('bind_sources', None)
            if cur_sources is None:
                cur_sources = []
            new_sources = (self.dynamic_sources if self.bind_sources is True else self.bind_sources)
            lopts['bind_sources'] = list(cur_sources) + list(new_sources)
        if self.pass_tokens:
            tokens = self.pass_tokens
            if tokens is True:
                tokens = self.token_list
            tokens = {
                k:os.environ.get(k)
                for k in tokens
            }
            tokens = {k:v for k,v in tokens.items() if v is not None}
            if len(tokens) > 0:
                lopts['env'] = tokens | lopts.get('env', {})
        return SingularityLauncher(
            self.container_path,
            *self.get_container_env_command(),
            mode=self.container_mode,
            port=self.port,
            bind=bind_paths,
            **lopts
            # forward bind paths and other conveniences
        )

    def resolve_address(self):
        return ('localhost', self.port) #TODO make this handle UDP or other communication modes too

    def launch_session(self):
        if self.temp_dir is None:
            if self._tmp is None:
                self._tmp = tempfile.TemporaryDirectory()
                self._tmp.__enter__()
                self.temp_dir = self._tmp.name
            else:
                self.temp_dir = self._tmp

        if self.session is None:
            self._managed_session = True
            if self.port is None:
                self.port = self.pick_port()
            self._launcher = self.get_launcher()
            self.session = self._launcher.launch()
            is_serving = self.session.stdout.readline()
            if self.initialization_delay_time is not None:
                time.sleep(self.initialization_delay_time)

        return self.resolve_address()

    def cleanup_session(self):
        if self._managed_session:
            if self._launcher is not None:
                self._launcher.terminate()
            else:
                self.session.kill()
            self._managed_session = False

    def __del__(self):
        try:
            self.cleanup_session()
        except:
            ...

    def get_io_files(self):
        id = str(uuid.uuid4())
        return (
            os.path.join(self.temp_dir, 'structures-' + id + '.npz'),
            os.path.join(self.temp_dir, 'config-' + id + '.json')
        )
    def setup_request(self, coords, order=None, tasks='energy', **opts) -> (list[Any], dict[Any, Any], Any):
        address = self.launch_session()
        structures, config = self.get_io_files()
        np.savez(structures,
                 atoms=self.atoms,
                 coordinates=coords * UnitsData.convert("BohrRadius", "Angstroms"),
                 charge=self.charge if self.charge is not None else 0,
                 spin=self.multiplicity if self.multiplicity is not None else 1)
        output_dir, config_name = os.path.split(config)
        output_file = config_name.replace('config-', 'results-').replace('.json', '.npz')
        state = {
            'method':'psience',
            'structures':structures,
            'tasks':tasks,
            'order':order,
            'energy_evaluator':self.evaluator,
            'output_dir':output_dir,
            'output_file':os.path.join(output_dir, output_file),
            'model_dir': self.model_dir
        } | opts

        dev.write_json(config, state)
        return ['evaluate', [config]], {'connection':address}, [config, structures]

    def process_response(self, response, state):
        # response from MLIP env is either an error or
        # outputs as specified in the config
        target_dir = os.path.dirname(state[0])
        if 'output_file' not in response:
            msg = response.get('stdout', '')
            err = response.get('stderr', '')
            if len(err) > 0:
                if len(msg) == 0:
                    msg = err
                else:
                    msg = msg + "\n" + err
            raise ValueError(f"no output from server, got message: {msg}")
        output_file = os.path.basename(response['output_file'])
        res = read_flat_tree(os.path.join(target_dir, output_file))
        return res['results'][0]

    def cleanup_state(self, response, state:list[str]):
        for file in state:
            try:
                os.remove(file)
            except IOError:
                ...

    def _run_mlip_request(self, *args, **kwargs):
        from mlipenv.client import MLIPHandler

        rc = self.session.poll()
        if rc is not None:
            output = (
                    self.session.stdout.read()
                    + self.session.stderr.read()
            )
            raise ValueError(f"container process terminated with status code {rc} and output {output}")

        return MLIPHandler.client_request(*args, **kwargs, print_response=None)

@EnergyEvaluator.register('function')
class PotentialFunctionEnergyEvaluator(EnergyEvaluator):
    # slightly weird, want single inheritance, but using `PropertyFunctionEvaluator` to define the
    # interface--I'm sure there's a better way to do this, but multiple inheritance always feels
    # fragile and python doesn't implement traits...

    def __init__(self,
                 potential_function,
                 property_units=None,
                 energy_units='Hartrees',
                 distance_units='BohrRadius',
                 **opts
                 ):
        """
        **LLM Docstring**

        Wrap a bare potential-energy function as an `EnergyEvaluator`, delegating to `PropertyFunctionEvaluator.__init__` for the actual state setup.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)`
        :type potential_function: callable
        :param property_units: the units the function's output is given in; defaults to `energy_units`
        :type property_units: str | None
        :param energy_units: the default energy units if `property_units` isn't given
        :type energy_units: str
        :param distance_units: the units the function expects its input coordinates in
        :type distance_units: str
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        if property_units is None:
            property_units = energy_units
        PropertyFunctionEvaluator.__init__(
            self,
            potential_function,
            property_units=property_units,
            distance_units=distance_units,
            **opts
        )

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped potential function, delegating to `PropertyFunctionEvaluator.evaluate_term`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        return PropertyFunctionEvaluator.evaluate_term(
            self,
            coords,
            order,
            **opts
        )

    def use_internal_coordinate_handlers(self):
        """
        **LLM Docstring**

        Whether this evaluator should route coordinate handling through the internal-coordinate machinery, delegating to `PropertyFunctionEvaluator.use_internal_coordinate_handlers`.

        :return: whether internal-coordinate handling applies
        :rtype: bool
        """
        return PropertyFunctionEvaluator.use_internal_coordinate_handlers(self)

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 fd_handler=None,
                 **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped potential function's expansion, delegating to `PropertyFunctionEvaluator.evaluate` (defaulting `fd_handler` to `internal_finite_difference_derivs` when working in internal coordinates).

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the highest derivative order to compute
        :type order: int
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param fd_handler: an explicit finite-difference handler to use
        :type fd_handler: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.evaluate`
        :type opts: dict
        :return: the evaluated energy expansion
        :rtype: list
        """
        if fd_handler is None and self.use_internal_coordinate_handlers():
            fd_handler = self.internal_finite_difference_derivs
        return PropertyFunctionEvaluator.evaluate(
            self,
            coords,
            order,
            fd_handler=fd_handler,
            **opts
        )

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        """
        **LLM Docstring**

        Resolve the potential function to wrap; on this class, simply returns `prop_func` unchanged alongside the given options.

        :param prop_func: the candidate potential function
        :type prop_func: callable
        :param mol: the molecule the evaluator is being built for (unused directly)
        :type mol: AbstractMolecule
        :param opts: extra options, passed through unchanged
        :type opts: dict
        :return: `(prop_func, opts)`
        :rtype: tuple
        """
        return prop_func, opts

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 potential_function=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build a `PotentialFunctionEnergyEvaluator` for a molecule, resolving the potential function from either `property_function` or `potential_function` and delegating to `PropertyFunctionEvaluator.initialize_from_mol`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the potential function to wrap
        :type property_function: callable | None
        :param potential_function: an alias for `property_function`, used if that isn't given
        :type potential_function: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.initialize_from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PotentialFunctionEnergyEvaluator
        """
        if property_function is None:
            property_function = potential_function
        return PropertyFunctionEvaluator.initialize_from_mol(
            cls,
            mol,
            property_function=property_function,
            **opts
        )

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default potential function to be used the next time `from_mol` is called without an explicit potential function.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        cls.default_property_function = potential

@EnergyEvaluator.register("expansion")
class PotentialExpansionEnergyEvaluator(PotentialFunctionEnergyEvaluator):

    def __init__(self,
                 expansion,
                 center=None,
                 ref=None,
                 batched_orders=None,
                 transforms=None,
                 transformed_derivatives=False,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build an energy evaluator from a precomputed Taylor-series potential expansion (rather than a live function), wrapping it in a `PotentialSurface` if it isn't already callable.

        :param expansion: the potential expansion, either a callable surface object or raw derivative tensors to build a `PotentialSurface` from
        :type expansion: callable | list[np.ndarray]
        :param center: the reference/expansion center for the surface, if building one from raw derivatives
        :type center: np.ndarray | None
        :param ref: the reference energy value for the surface
        :type ref: float | None
        :param batched_orders: accepted but overridden to `True` (expansions always return all requested orders in one call)
        :type batched_orders: bool | None
        :param transforms: coordinate transformation(s) to apply when building the surface from raw derivatives
        :type transforms: list | None
        :param transformed_derivatives: whether the raw derivatives are already expressed in the transformed coordinate system
        :type transformed_derivatives: bool
        :param opts: extra options forwarded to `PotentialFunctionEnergyEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        if not callable(expansion):
            expansion = PotentialSurface.from_derivatives(expansion,
                                                          center=center, ref=ref,
                                                          transforms=transforms,
                                                          transformed_derivatives=transformed_derivatives
                                                          )
        super().__init__(
            expansion,
            batched_orders=True,
            **opts
        )

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 expansion=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build a `PotentialExpansionEnergyEvaluator` for a molecule, resolving the expansion from either `property_function` or `expansion`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the potential expansion to use
        :type property_function: object | None
        :param expansion: an alias for `property_function`, used if that isn't given
        :type expansion: object | None
        :param opts: extra options forwarded to the base class's `from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PotentialExpansionEnergyEvaluator
        """
        return super().from_mol(
            mol,
            property_function=expansion if property_function is None else property_function,
            **opts
        )

    @classmethod
    def get_property_function(cls, expansion, mol, transforms=None, **ignored):
        """
        **LLM Docstring**

        Resolve the potential expansion to use: falls back to the molecule's own `potential_derivatives` if none is given, re-expresses the derivatives in a normal-mode basis (via the molecule's normal modes) if the cubic term appears to already be in a reduced mode basis, and wraps the result in a `PotentialSurface`.

        :param expansion: the candidate expansion (raw derivative tensors, `None` to use the molecule's own, or an already-callable surface)
        :type expansion: object | None
        :param mol: the molecule the evaluator is being built for
        :type mol: AbstractMolecule
        :param transforms: an explicit coordinate transformation to build the surface with; derived from the molecule's normal modes if needed and not given
        :type transforms: list | None
        :param ignored: extra options, passed through unchanged as the returned options
        :type ignored: dict
        :return: `(expansion, ignored)` -- the resolved (possibly newly built) `PotentialSurface` and the untouched extra options
        :rtype: tuple
        """
        if not callable(expansion):
            transforms = transforms
            transformed_derivatives = False
            if expansion is None:
                expansion = mol.potential_derivatives
                # handle partial quatics
            if len(expansion) > 2 and expansion[2].shape[0] < expansion[1].shape[0]:
                if transforms is None:
                    modes = mol.get_normal_modes(use_internals=False, project_transrot=False).remove_mass_weighting()
                    transforms = [[modes.coords_by_modes], [modes.modes_by_coords]]
                if nput.is_numeric_array_like(transforms[0]):
                    tf = np.asanyarray(transforms[0])
                    if tf.ndim > 2:
                        tf = tf[0]
                else:
                    tf = np.asanyarray(transforms[0][0])
                _ = []
                for i,e in enumerate(expansion):
                    for j in range(min([i+1, 2])):
                        e = np.tensordot(tf, e, axes=[1, -1])
                    _.append(e)
                expansion = _
                transformed_derivatives = True
            expansion = PotentialSurface.from_mol(mol,
                                                  expansion=expansion,
                                                  transforms=transforms,
                                                  transformed_derivatives=transformed_derivatives
                                                  )
        return expansion, ignored

class DipoleEvaluator(PropertyEvaluator):
    evaluator_registry = {}

    target_property_units = ("ElementaryCharge", "BohrRadius")
    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The built-in named dipole-evaluator constructors: an expansion-based evaluator, plus adapters for HIPNN, AIMNet2, and RDKit.

        :return: the evaluator-name-to-constructor mapping
        :rtype: dict
        """
        return {
            'expansion': DipoleExpansionDipoleEvaluator,
            'hipnn': HIPNNDipoleEvaluator,
            'aimnet2': AIMNet2DipoleEvaluator,
            'rdkit': RDKitDipoleEvaluator,
        }
    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as a dipole-evaluator specification.

        :return: `DipoleFunctionDipoleEvaluator`
        :rtype: type
        """
        return DipoleFunctionDipoleEvaluator

class DipoleFunctionDipoleEvaluator(DipoleEvaluator):

    def __init__(self,
                 potential_function,
                 property_units=('ElementaryCharge', 'BohrRadius'),
                 distance_units='BohrRadius',
                 **opts
                 ):
        """
        **LLM Docstring**

        Wrap a bare dipole-evaluating function as a `DipoleEvaluator`, delegating to `PropertyFunctionEvaluator.__init__` for the actual state setup.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)` and return dipole values/derivatives
        :type potential_function: callable
        :param property_units: the units the function's output is given in
        :type property_units: tuple | str
        :param distance_units: the units the function expects its input coordinates in
        :type distance_units: str
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        PropertyFunctionEvaluator.__init__(
            self,
            potential_function,
            property_units=property_units,
            distance_units=distance_units,
            **opts
        )

    def embed_coords(self, coords, **kwargs):
        """
        **LLM Docstring**

        Embed coordinates into this evaluator's working representation, delegating to `PropertyFunctionEvaluator.embed_coords`.

        :param coords: the coordinates to embed
        :type coords: np.ndarray
        :param kwargs: extra options forwarded to `PropertyFunctionEvaluator.embed_coords`
        :type kwargs: dict
        :return: the embedded coordinates
        :rtype: np.ndarray
        """
        return PropertyFunctionEvaluator.embed_coords(
            self,
            coords,
            **kwargs
        )

    def unembed_coords(self, coords, **kwargs):
        """
        **LLM Docstring**

        Convert coordinates back to plain Cartesians, delegating to `PropertyFunctionEvaluator.unembed_coords`.

        :param coords: the coordinates to unembed
        :type coords: np.ndarray
        :param kwargs: extra options forwarded to `PropertyFunctionEvaluator.unembed_coords`
        :type kwargs: dict
        :return: the unembedded (Cartesian) coordinates
        :rtype: np.ndarray
        """
        return PropertyFunctionEvaluator.unembed_coords(
            self,
            coords,
            **kwargs
        )

    def unembed_derivs(self, base_coords, coords, derivs, **kwargs):
        """
        **LLM Docstring**

        Convert dipole derivative tensors back to plain Cartesian coordinates, delegating to `PropertyFunctionEvaluator.unembed_derivs`.

        :param base_coords: the reference coordinates the returned derivatives should be expressed relative to
        :type base_coords: np.ndarray
        :param coords: the coordinates the input derivatives were computed at
        :type coords: np.ndarray
        :param derivs: the derivative tensors to convert
        :type derivs: list[np.ndarray]
        :param kwargs: extra options forwarded to `PropertyFunctionEvaluator.unembed_derivs`
        :type kwargs: dict
        :return: the converted derivative tensors
        :rtype: list[np.ndarray]
        """
        return PropertyFunctionEvaluator.unembed_derivs(
            self,
            base_coords,
            coords,
            derivs,
            **kwargs
        )

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped dipole function, delegating to `PropertyFunctionEvaluator.evaluate_term`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        return PropertyFunctionEvaluator.evaluate_term(
            self,
            coords,
            order,
            **opts
        )

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped dipole function's expansion, delegating to `PropertyFunctionEvaluator.evaluate`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the highest derivative order to compute
        :type order: int
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.evaluate`
        :type opts: dict
        :return: the evaluated dipole expansion
        :rtype: list
        """
        return PropertyFunctionEvaluator.evaluate(
            self,
            coords,
            order,
            **opts
        )

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        """
        **LLM Docstring**

        Resolve the dipole function to wrap; on this class, simply returns `prop_func` unchanged alongside the given options.

        :param prop_func: the candidate dipole function
        :type prop_func: callable
        :param mol: the molecule the evaluator is being built for (unused directly)
        :type mol: AbstractMolecule
        :param opts: extra options, passed through unchanged
        :type opts: dict
        :return: `(prop_func, opts)`
        :rtype: tuple
        """
        return prop_func, opts

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 dipole_function=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build a `DipoleFunctionDipoleEvaluator` for a molecule, resolving the dipole function from either `property_function` or `dipole_function` and delegating to `PropertyFunctionEvaluator.initialize_from_mol`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the dipole function to wrap
        :type property_function: callable | None
        :param dipole_function: an alias for `property_function`, used if that isn't given
        :type dipole_function: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.initialize_from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: DipoleFunctionDipoleEvaluator
        """
        if property_function is None:
            property_function = dipole_function
        return PropertyFunctionEvaluator.initialize_from_mol(
            cls,
            mol,
            property_function=property_function,
            **opts
        )

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default dipole function to be used the next time `from_mol` is called without an explicit dipole function.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        cls.default_property_function = potential

class DipoleExpansionDipoleEvaluator(DipoleFunctionDipoleEvaluator):

    def __init__(self,
                 expansion,
                 center=None,
                 ref=None,
                 property_units=('ElementaryCharge', "BohrRadius"),
                 distance_units='BohrRadius',
                 analytic_derivative_order=None,
                 batched_orders=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build a dipole evaluator from a precomputed Taylor-series dipole expansion (rather than a live function), wrapping it in a `DipoleSurface` if it isn't already callable and defaulting the analytic derivative order to the number of terms in the expansion.

        :param expansion: the dipole expansion, either a callable surface object or raw derivative tensors to build a `DipoleSurface` from
        :type expansion: callable | list[np.ndarray]
        :param center: the reference/expansion center for the surface
        :type center: np.ndarray | None
        :param ref: the reference dipole value for the surface
        :type ref: np.ndarray | None
        :param property_units: the units the dipole values are given in
        :type property_units: tuple
        :param distance_units: the units the coordinates are given in
        :type distance_units: str
        :param analytic_derivative_order: the highest order available analytically; defaults to the number of terms in `expansion`
        :type analytic_derivative_order: int | None
        :param batched_orders: accepted but overridden to `True`
        :type batched_orders: bool | None
        :param opts: extra options forwarded to `DipoleFunctionDipoleEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        if not callable(expansion):
            expansion = DipoleSurface.from_derivatives(expansion, center=center, ref=ref)
        if analytic_derivative_order is None:
            analytic_derivative_order = len(expansion.expansion_tensors)
        super().__init__(
            expansion,
            property_units=property_units,
            distance_units=distance_units,
            analytic_derivative_order=analytic_derivative_order,
            batched_orders=True,
            **opts
        )

    @classmethod
    def expansion_from_mol_charges(cls, mol, charges=None, origin=None):
        """
        **LLM Docstring**

        Build a first-order dipole expansion (a constant term and its Cartesian gradient) directly from a set of atomic partial charges, treating the dipole as the charge-weighted sum of displacements from a reference origin.

        :param mol: the molecule to compute the dipole expansion for
        :type mol: AbstractMolecule
        :param charges: the per-atom partial charges; computed via `mol.calculate_charges()` if not given
        :type charges: np.ndarray | None
        :param origin: the reference origin for the dipole; defaults to the molecule's center of mass
        :type origin: np.ndarray | None
        :return: `[dipole_value, gradient]` -- the zeroth- and first-order dipole expansion terms
        :rtype: list[np.ndarray]
        """
        if charges is None:
            charges = mol.calculate_charges()
        charges = np.asanyarray(charges)
        if origin is None:
            origin = mol.center_of_mass
        coords = mol.coords
        disps = coords - np.asanyarray(origin)[np.newaxis]
        dip_contribs = charges[:, np.newaxis] * disps
        deriv = np.zeros(dip_contribs.shape + (3,), dtype=dip_contribs.dtype)
        for i in range(3):
            deriv[:, i, i] = charges
        return [np.sum(dip_contribs, axis=0), deriv.reshape(-1, 3)]

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 expansion=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build a `DipoleExpansionDipoleEvaluator` for a molecule, resolving the dipole expansion from either `property_function` or `expansion`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the dipole expansion to use
        :type property_function: object | None
        :param expansion: an alias for `property_function`, used if that isn't given
        :type expansion: object | None
        :param opts: extra options forwarded to the base class's `from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: DipoleExpansionDipoleEvaluator
        """
        return super().from_mol(
            mol,
            property_function=expansion
            if property_function is None else property_function,
            **opts
        )
    @classmethod
    def get_property_function(cls, expansion, mol, center=None, transforms=None, use_modes=True, **ignored):
        """
        **LLM Docstring**

        Resolve the dipole expansion to use: falls back to the molecule's own `dipole_derivatives` (or a charge-based estimate via `expansion_from_mol_charges`, if neither is available) if none is given, optionally re-expresses the expansion through the molecule's normal modes (via its Hamiltonian's `dipole_expansion`), and wraps the result in a `DipoleSurface`.

        :param expansion: the candidate dipole expansion (raw derivative tensors, `None` to use the molecule's own/charge-derived one, or an already-callable surface)
        :type expansion: object | None
        :param mol: the molecule the evaluator is being built for
        :type mol: AbstractMolecule
        :param center: the reference/expansion center to build the surface with
        :type center: np.ndarray | None
        :param transforms: an explicit coordinate transformation to use; derived from the molecule's normal-mode embedding if `use_modes` is set and none is given
        :type transforms: list | None
        :param use_modes: whether to re-express the expansion in the molecule's normal-mode basis before building the surface
        :type use_modes: bool
        :param ignored: extra options, passed through unchanged as the returned options
        :type ignored: dict
        :return: `(expansion, ignored)` -- the resolved (possibly newly built) `DipoleSurface` and the untouched extra options
        :rtype: tuple
        """
        if not callable(expansion):
            transformed_derivatives = False
            if expansion is None:
                expansion = mol.dipole_derivatives
                if expansion is None:
                    expansion = cls.expansion_from_mol_charges(mol)
                # handle partial quatics

            if use_modes:
                dip_exp = mol.hamiltonian.dipole_expansion(expansion=expansion)
                expansion = dip_exp.get_terms(transformation=transforms)
                transformed_derivatives = True
                if transforms is None:
                    exp = dip_exp.embedding.get_coords_by_modes(mass_weighted=False)
                    if exp is not None:
                        transforms = [
                            [exp],
                            [dip_exp.embedding.get_modes_by_coords(mass_weighted=False)]
                        ]
                    # print("FWD:", [f.shape for f in transforms[0]])
                    # print("REV:", [f.shape for f in transforms[1]])


                # if len(expansion) > 2 and expansion[2].shape[0] < expansion[1].shape[0]:
                #     raise Exception("?")
                #     if transforms is None:
                #         modes = mol.get_normal_modes(use_internals=False, project_transrot=False).remove_mass_weighting()
                #         transforms = [[modes.coords_by_modes], [modes.modes_by_coords]]
                #     if nput.is_numeric_array_like(transforms[0]):
                #         tf = np.asanyarray(transforms[0])
                #         if tf.ndim > 2:
                #             tf = tf[0]
                #     else:
                #         tf = np.asanyarray(transforms[0][0])
                #     _ = []
                #     for i, e in enumerate(expansion[1:]):
                #         e = np.tensordot(tf, e, axes=[1, -2])
                #         _.append(e)
                #     expansion = [expansion[0]] + _
                #     transformed_derivatives = True

            expansion = DipoleSurface.from_mol(mol,
                                               center=center,
                                               expansion=expansion,
                                               transforms=transforms,
                                               transformed_derivatives=transformed_derivatives
                                               )
        return expansion, ignored

class HIPNNDipoleEvaluator(DipoleEvaluator):

    def __init__(self, atoms, model_dir, charge=0, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.base = HIPNNEnergyEvaluator(atoms, model_dir, charge=charge, multiplicity=multiplicity, quiet=quiet,
                                         property_key='dipole')
        self.base.gradient_key = None

    @classmethod
    def from_mol(cls, mol, *, model_dir, **opts):
        return cls(mol.atoms, model_dir=model_dir, **cls.prep_mol_opts(mol, **opts))

    property_units = ("ElementaryCharge", "Angstroms")
    distance_units = 'Angstroms'
    batched_orders = True
    analytic_derivative_order = 0
    def evaluate_term(self, coords, order, **opts):
        base_shape, coords, data = self.base.prep_eval(coords, order, **opts)

        expansions = self.base.process_results(
            base_shape,
            coords,
            data,
            order,
            forces=order > 0,
            property_key='dipole',
            gradient_key='dipole_grad',
            property_shape=(coords.shape[-1],)
        )

        return expansions

class ChargeEvaluator(PropertyEvaluator):
    evaluator_registry = {}

    target_property_units = "ElementaryCharge"
    default_evaluator_type = 'rdkit'
    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The built-in named charge-evaluator constructors: adapters for RDKit and AIMNet2.

        :return: the evaluator-name-to-constructor mapping
        :rtype: dict
        """
        return {
            'rdkit': RDKitChargeEvaluator,
            'aimnet2': AIMNet2ChargeEvaluator
        }
    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as a charge-evaluator specification.

        :return: `ChargeFunctionChargeEvaluator`
        :rtype: type
        """
        return ChargeFunctionChargeEvaluator

class ChargeFunctionChargeEvaluator(ChargeEvaluator):
    def __init__(self,
                 potential_function,
                 property_units='ElementaryCharge',
                 distance_units='BohrRadius',
                 **opts
                 ):
        """
        **LLM Docstring**

        Wrap a bare charge-evaluating function as a `ChargeEvaluator`, delegating to `PropertyFunctionEvaluator.__init__` for the actual state setup.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)` and return partial-charge values/derivatives
        :type potential_function: callable
        :param property_units: the units the function's output is given in
        :type property_units: str
        :param distance_units: the units the function expects its input coordinates in
        :type distance_units: str
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        PropertyFunctionEvaluator.__init__(
            self,
            potential_function,
            property_units=property_units,
            distance_units=distance_units,
            **opts
        )

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped charge function, delegating to `PropertyFunctionEvaluator.evaluate_term`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        return PropertyFunctionEvaluator.evaluate_term(
            self,
            coords,
            order,
            **opts
        )

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped charge function's expansion, delegating to `PropertyFunctionEvaluator.evaluate`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the highest derivative order to compute
        :type order: int
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.evaluate`
        :type opts: dict
        :return: the evaluated charge expansion
        :rtype: list
        """
        return PropertyFunctionEvaluator.evaluate(
            self,
            coords,
            order,
            **opts
        )

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 **opts
                 ):
        """
        **LLM Docstring**

        Build a `ChargeFunctionChargeEvaluator` for a molecule, delegating to `PropertyFunctionEvaluator.initialize_from_mol`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the charge function to wrap
        :type property_function: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.initialize_from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: ChargeFunctionChargeEvaluator
        """
        return PropertyFunctionEvaluator.initialize_from_mol(
            cls,
            mol,
            property_function=property_function,
            **opts
        )

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default charge function to be used the next time `from_mol` is called without an explicit charge function.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        cls.default_property_function = potential

class RDKitChargeEvaluator(ChargeEvaluator):
    def __init__(self, rdmol, model='gasteiger', **defaults):
        super().__init__(**defaults)
        self.rdmol = rdmol
        self.model = model

    @classmethod
    def from_mol(cls, mol, charge=None, multiplicity=None, **opts):
        return cls(mol.rdmol, **cls.prep_mol_opts(mol, **opts))

    property_units = 'ElementaryCharge'
    analytic_derivative_order = 0
    def evaluate_term(self, coords, order, model=None, **opts):
        if model is None:
            model = self.model
        return np.array(
            self.rdmol.evaluate_charges(coords, model=model, **opts)
        ) / np.sqrt(2) # TODO: I think my base dipoles are just scaled incorrectly...

class AIMNet2ChargeEvaluator(ChargeEvaluator):

    def __init__(self, atoms, model='aimnet2', charge=0, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.base = AIMNet2EnergyEvaluator(atoms, model=model, charge=charge, multiplicity=multiplicity, quiet=quiet)

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    property_units = 'ElementaryCharge'
    distance_units = 'Angstroms'
    batched_orders = True
    analytic_derivative_order = 6
    scaling_factor = 1#0.59667 # necessary to make water work...but I don't know what the actual units are
    def evaluate_term(self, coords, order, **opts):
        base_shape, coords, data = self.base.prep_eval(coords, order, keep_graph=order>0, forces=False, hessian=False, **opts)

        expansions = self.base.process_aimnet_derivs(
            base_shape, coords, data,
            'charges',
            order,
            property_shape=(coords.shape[-2]+1,),
            output_shape=(coords.shape[-2],)
        )
        s = self.scaling_factor / np.sqrt(2)
        expansions = [e * s for e in expansions]

        return expansions

class ChargeEvaluatorDipoleEvaluator(DipoleEvaluator):
    def __init__(self, charge_evaluator:'ChargeEvaluator', **etc):
        self.evaluator = charge_evaluator
        super().__init__(**etc)
        self.distance_units = self.evaluator.distance_units
        self.property_units = (self.evaluator.property_units, self.evaluator.distance_units)
        self.analytic_derivative_order = self.evaluator.analytic_derivative_order

    batched_orders = True
    def evaluate_term(self, coords, order, **opts):
        coords = np.asanyarray(coords)
        charges = self.evaluator.evaluate(coords, order=order, **opts)
        coord_deriv = nput.identity_tensors(coords.shape[:-2], coords.shape[-2]*coords.shape[-1])
        coord_deriv = coord_deriv.reshape(coord_deriv.shape[:-2] + (coord_deriv.shape[-2], -1, 3))
        coord_expansion = [coords, coord_deriv]
        expansion = nput.tensordot_deriv(
            charges, coord_expansion,
            order=order,
            axes=[-1, -2]
        )
        return expansion

class AIMNet2DipoleEvaluator(ChargeEvaluatorDipoleEvaluator):

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(AIMNet2ChargeEvaluator.from_mol(mol, **opts))

class RDKitDipoleEvaluator(ChargeEvaluatorDipoleEvaluator):

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(RDKitChargeEvaluator.from_mol(mol, **opts))

class DipolePolarizabilityEvaluator(PropertyEvaluator):
    target_property_units = "ElementaryCharge"

    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The built-in named dipole-polarizability-evaluator constructors: adapters for AIMNet2 and HIPNN.

        :return: the evaluator-name-to-constructor mapping
        :rtype: dict
        """
        return {
            'aimnet2': AIMNet2DipolePolarizabilityEvaluator,
            'hipnn': HIPNNDipolePolarizabilityEvaluator,
        }

    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as a dipole-polarizability-evaluator specification.

        :return: `PolarizabilityFunctionDipolePolarizabilityEvaluator`
        :rtype: type
        """
        return PolarizabilityFunctionDipolePolarizabilityEvaluator

    @abc.abstractmethod
    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, **opts):
        """
        **LLM Docstring**

        Abstract hook for evaluating the dipole-polarizability tensor's value (or specific coordinate/electrostatic-field derivative orders) at the given coordinates. Concrete evaluator subclasses must implement this.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param coord_order: the coordinate-derivative order to evaluate
        :type coord_order: int
        :param electrostatic_order: the electrostatic-field-derivative order to evaluate
        :type electrostatic_order: int
        :param opts: extra evaluator-specific options
        :type opts: dict
        :return: the evaluated term(s)
        :rtype: object
        """
        ...

    multi_expansion_order = 2
    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the polarizability's coordinate-derivative expansion at a given coordinate order, delegating to `evaluate_polarizability_term` with a fixed electrostatic order (`self.multi_expansion_order`) unless `order` is already given as an explicit `(coord_order, electrostatic_order)` pair.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the coordinate-derivative order, or an explicit `(coord_order, electrostatic_order)` pair
        :type order: int | tuple[int, int]
        :param opts: extra options forwarded to `evaluate_polarizability_term`
        :type opts: dict
        :return: the evaluated term(s)
        :rtype: object
        """
        if nput.is_int(order):
            order = (order, self.multi_expansion_order)
        return self.evaluate_polarizability_term(coords, *order, **opts)

class PolarizabilityFunctionDipolePolarizabilityEvaluator(DipolePolarizabilityEvaluator):
    ...

class FluctuatingChargeDipolePolarizabilityEvaluator(DipolePolarizabilityEvaluator):

    @classmethod
    def mu_factor(cls, total_charge, X, J_inv):
        return (total_charge + np.sum(np.dot(J_inv, X), axis=0)) / np.sum(J_inv)

    @classmethod
    def dipole_factor_expansion(cls, coords, total_charge, X, J_inv):

        charge_factor = np.dot(J_inv, X) + cls.mu_factor(total_charge, X, J_inv) * np.sum(J_inv, axis=-1)
        return np.tensordot(coords, charge_factor, axes=[-2, 0])

    @classmethod
    def polarizability_factor(cls, coords, J_inv):
        coord_factor = np.tensordot(coords, J_inv, axes=[-2, -1])
        sum_factor = np.sum(coord_factor, axis=0)
        return (
                np.tensordot(coord_factor, coords, axes=[-1, -2])
                + sum_factor[:, np.newaxis] * sum_factor[np.newaxis, :] / np.sum(J_inv)
        )

    @abc.abstractmethod
    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, **opts):
        ...

    batched_orders = True
    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, *, charge, **opts):
        potential_by_charges_by_coords = self.evaluate_electrostatic_potential_derivatives(coords, coord_order, electrostatic_order, **opts)
        return self.expand_electrostatic_potential(
            charge, coords, potential_by_charges_by_coords
        )

    @classmethod
    def expand_electrostatic_potential(cls, total_charge, coords, potential_by_charges_by_coords):
        expansion = []
        es_order = len(potential_by_charges_by_coords) - 1
        order = len(potential_by_charges_by_coords[0]) - 1
        if es_order > 2:
            raise NotImplementedError("hyperpolarizabilities require more electrostatic derivatives than I currently have")
        if es_order > 0:
            # pull electrostatic derivatives
            chi_expansion = [e.astype('float64') for e in potential_by_charges_by_coords[1]]
            J_expansion = [e.astype('float64') for e in potential_by_charges_by_coords[2]]
            J_inv_expansion = nput.matinv_deriv(J_expansion, order=order)
            # prep coordinate expansion
            coords = np.asanyarray(coords)
            coord_deriv = nput.identity_tensors(coords.shape[:-2], coords.shape[-2]*coords.shape[-1])
            coord_deriv = coord_deriv.reshape(coord_deriv.shape[:-2] + (coord_deriv.shape[-2], -1, 3))
            coord_expansion = [coords, coord_deriv]
            # construct terms that are common throughout expansions
            J_inv_sum_1 = [np.sum(e, axis=-1) for e in J_inv_expansion]
            J_inv_tot = [np.sum(e, axis=-1) for e in J_inv_sum_1]
            J_inv_tot_inv = nput.scalarinv_deriv(J_inv_tot, order=order)

            # (total_charge + np.sum(np.dot(J_inv, X), axis=0)) / np.sum(J_inv)
            # charge_factor = np.dot(J_inv, X) + cls.mu_factor(total_charge, X, J_inv) * np.sum(J_inv, axis=-1)
            # return np.tensordot(coords, charge_factor, axes=[-2, 0])
            J_inv_chi_expansion = nput.tensordot_deriv(chi_expansion, J_inv_expansion, order=order, axes=[-1, -1])
            J_inv_chi_sum_expansion = nput.tensordot_deriv(J_inv_sum_1, chi_expansion, order=order, axes=[-1, -1])
            mu_expansion =  nput.scalarprod_deriv(J_inv_chi_sum_expansion, J_inv_tot_inv, order=order)
            mu_expansion[0] +=  total_charge * J_inv_tot_inv[0]

            mu_prod_expansion = nput.scalarprod_deriv(mu_expansion, J_inv_sum_1, order=order)
            dipole_expansion_charge_factor = nput.add_expansions(J_inv_chi_expansion, mu_prod_expansion)
            dipole_expansion = nput.tensordot_deriv(dipole_expansion_charge_factor, coord_expansion, axes=[-1, -2], order=order)
            expansion.append(dipole_expansion)

            if es_order > 1:
                # pad_j = [
                #     np.repeat(np.expand_dims(np.moveaxis(j, -1, 0), 0), 3, axis=0)
                #     for j in J_inv_expansion
                # ]
                # pad_c = [
                #     np.moveaxis(np.repeat(np.expand_dims(c, 0), len(J_expansion[0]), axis=0), -1, 0)
                #     for c in coord_expansion
                # ]
                # print([j.shape for j in pad_j], [c.shape for c in pad_c])
                J_inv_X_expansion = nput.tensordot_deriv(
                    J_inv_expansion,
                    coord_expansion,
                    order=order,
                    axes=[-1, -2],
                    shared=coords.ndim-2
                )
                X_J_inv_X_expansion = nput.tensordot_deriv(
                    J_inv_X_expansion,
                    coord_expansion,
                    order=order,
                    axes=[-2, -2]
                )

                J_inv_X_1_expansion = [np.sum(jj, axis=-2) for jj in J_inv_X_expansion]
                JX_x_JX_expansion = nput.tensorprod_deriv(J_inv_X_1_expansion, J_inv_X_1_expansion, order=order,
                                                          axes=[-1, -1]
                                                          )
                outer_expansion = nput.scalarprod_deriv(JX_x_JX_expansion, J_inv_tot_inv, order=order)
                polarizability_expansion = nput.add_expansions(
                    X_J_inv_X_expansion,
                    outer_expansion
                )
                expansion.append(polarizability_expansion)

        return expansion

    @classmethod
    def _distorted_harmonic_polarizabilities(cls, charges, dipole_expansion, coords, hessians):
        coords = np.asanyarray(coords)
        charges = np.asanyarray(charges)
        hessians = np.asanyarray(hessians)
        x = charges[..., :, np.newaxis] * coords
        nats = charges.shape[-1]
        hessians = hessians.reshape(hessians.shape[:-2] + (nats, 3, nats, 3))
        f_factors = np.array([
            np.sum(x[..., i] ** 2) /
            nput.vec_tensordot(
                nput.vec_tensordot(
                    hessians[..., :, i, :, i],
                    x[..., i],
                    axes=[-1, -1],
                    shared=x.ndim - 2
                ),
                x[..., i],
                axes=[-1, -1],
                shared=x.ndim - 2
            )
            for i in range(3)
        ])
        dx_f = x * np.moveaxis(f_factors, 0, -1)[..., np.newaxis, :]
        pad_derivs = np.zeros(x.shape[:-2] + (x.shape[-1] * x.shape[-2], 3))
        for i in range(3):
            pad_derivs[..., i::3, i] = dx_f[..., i]
        dx_expansion = [pad_derivs]
        return nput.tensordot_deriv(dx_expansion, dipole_expansion[1:], axes=[-2, -2], order=len(dipole_expansion) - 2)

class AIMNet2DipolePolarizabilityEvaluator(FluctuatingChargeDipolePolarizabilityEvaluator):
    def __init__(self, atoms, model='aimnet2', charge=0, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.base = AIMNet2EnergyEvaluator(atoms, model=model, charge=charge, multiplicity=multiplicity, quiet=quiet)

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, **cls.prep_mol_opts(mol, **opts))

    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, *, charge=None, **opts):
        if charge is None:
            charge = self.base.charge
        return super().evaluate_polarizability_term(
            coords, coord_order, electrostatic_order, charge=charge, **opts
        )

    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, **opts):
        self.base.eval.model.outputs.lrcoulomb.key_out = 'coul_energy'
        base_shape, coords, data = self.base.prep_eval(coords, 0,
                                                       keep_graph=electrostatic_order>0 or coord_order>0,
                                                       **opts)
        data['coul_energy'] = data['coul_energy'].sum()

        electrostatic_expansions = self.base.process_aimnet_derivs(
            base_shape, coords, data,
            'coul_energy',
            coord_order,
            property_shape=data['coul_energy'].shape[len(base_shape):],
            extra_deriv_coord=(data['charges'], electrostatic_order)
        )

        return electrostatic_expansions

class HIPNNDipolePolarizabilityEvaluator(FluctuatingChargeDipolePolarizabilityEvaluator):
    def __init__(self, atoms, model_dir, **defaults):
        super().__init__()
        self.base = HIPNNEnergyEvaluator(atoms, model_dir, **defaults)
        self._set_extra_outputs(self.base)

    @classmethod
    def _set_extra_outputs(cls, base):
        import hippynn.graphs
        predictor = hippynn.graphs.Predictor.from_graph(base.eval.graph,
                                                        additional_outputs=[
                                                            'ChEQ.coul_energy'
                                                        ]
                                                        )
        cheq = predictor.graph.node_from_name('ChEQ')
        predictor.outputs.extend(cheq.parents[-2:])
        base.eval = predictor

    @classmethod
    def from_mol(cls, mol, *, model_dir, embedding=None, charge=None, multiplicity=None, **opts):
        return cls(mol.atoms, model_dir=model_dir, **cls.prep_mol_opts(mol, opts))

    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, *, charge=None, **opts):
        if charge is None:
            charge = self.base.charge
        return super().evaluate_polarizability_term(
            coords, coord_order, electrostatic_order, charge=charge, **opts
        )

    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, **opts):
        self.base.eval.model.outputs.lrcoulomb.key_out = 'coul_energy'
        base_shape, coords, data = self.base.prep_eval(coords, 0,
                                                       keep_graph=electrostatic_order>0 or coord_order>0,
                                                       **opts)
        data['coul_energy'] = data['coul_energy'].sum()

        electrostatic_expansions = self.base.process_aimnet_derivs(
            base_shape, coords, data,
            'coul_energy',
            coord_order,
            property_shape=data['coul_energy'].shape[len(base_shape):],
            extra_deriv_coord=(data['charges'], electrostatic_order)
        )

        return electrostatic_expansions

class ReducedDimensionalPotentialHandler:

    def __init__(self, mol):
        """
        **LLM Docstring**

        Store the molecule this handler will build reduced-dimensionality (1D) potential slices for.

        :param mol: the molecule to build potentials for
        :type mol: AbstractMolecule
        :return: None
        :rtype: None
        """
        self.mol = mol

    # @property
    # def cartesian_derivatives(self):
    #     ...

    @classmethod
    def get_potential_params(cls, spec, re, g, local_derivs,
                             quartic_potential_cutoff=5e-2,
                             poly_expansion_order=2,
                             **opts):
        """
        **LLM Docstring**

        Fit a simple analytic 1D potential form (a Morse potential, or a low-order polynomial in a Duschinsky-like `w`-scaled coordinate) to a coordinate's local derivative expansion, choosing Morse when the cubic-to-quartic-derivative ratio exceeds `quartic_potential_cutoff` and a polynomial fit otherwise.

        :param spec: the coordinate specification the fit is for (passed through into the returned options for later use by a potential constructor)
        :type spec: object
        :param re: the coordinate's equilibrium/reference value
        :type re: float
        :param g: the coordinate's G-matrix element (effective inverse mass)
        :type g: float
        :param local_derivs: the local derivative expansion (value, 1st, 2nd, 3rd, ... derivatives) of the potential along this coordinate
        :type local_derivs: list[float]
        :param quartic_potential_cutoff: the cubic/quartic-derivative ratio threshold below which a polynomial (rather than Morse) fit is used
        :type quartic_potential_cutoff: float
        :param poly_expansion_order: the polynomial expansion order to keep when fitting the `'poly'` form
        :type poly_expansion_order: int
        :param opts: extra options passed through into the returned parameter dict
        :type opts: dict
        :return: the fitted potential parameters (a `'method'` key of `'morse'` or `'poly'`, plus the corresponding fit parameters and `re`)
        :rtype: dict
        """
        # if (use_morse and len(local_derivs) > 3) or (use_morse is None and len(spec) == 2):
        f2, f3, f4 = [local_derivs[i + 1] if i < len(local_derivs) else 0 for i in range(3)]
        if np.abs(f3 / f4) < quartic_potential_cutoff:
            opts['method'] = 'poly'
            opts['g'] = g
            sg = np.sqrt(g)
            opts['w_coeffs'] = [
                         sg * np.sign(d) * np.power(np.abs(d), 1 / (k + 2))
                         # * (.25 ** (k + 2))
                         for k, d in enumerate(local_derivs[1:4])
                     ][:poly_expansion_order - 1]
        else:
            opts['method'] = 'morse'
            opts['g'] = g
            sg = np.sqrt(g)
            w = sg * np.sqrt(f2)
            wx = -((g / (4 * w)) ** 2) * (f4 - 5 / 3 * (f3**2) / f2)
            params = [w, wx]
            # if len(local_derivs) > 4:
            #     params = params + [
            #         sg * np.sign(d) * np.power(np.abs(d) / math.factorial(k + 4), 1 / (k + 4))  # * (.25 ** (k + 2))
            #         for k, d in enumerate(local_derivs[4:])
            #     ]
            opts['w'] = w
            opts['wx'] = wx
        opts['re'] = re

        return opts

    @classmethod
    def get_morse_potential(cls, spec, *, w, wx, g, re):
        """
        **LLM Docstring**

        Build a Morse-potential `CoordinateFunction` for a given coordinate specification and fitted parameters.

        :param spec: the coordinate specification to build the potential over
        :type spec: object
        :param w: the Morse frequency parameter
        :type w: float
        :param wx: the Morse anharmonicity parameter
        :type wx: float
        :param g: the coordinate's G-matrix element
        :type g: float
        :param re: the equilibrium/reference coordinate value
        :type re: float
        :return: the constructed Morse potential function
        :rtype: CoordinateFunction
        """
        return CoordinateFunction.morse(
                spec,
                re=re,
                w=w,
                wx=wx,
                g=g
            )
    @classmethod
    def get_poly_potential(cls, spec, *, w_coeffs=None, coeffs=None, g, re):
        """
        **LLM Docstring**

        Build a polynomial `CoordinateFunction` for a given coordinate specification, either from raw Taylor coefficients or from `w`-scaled (Duschinsky-like) coefficients that are first converted into raw coefficients.

        :param spec: the coordinate specification to build the potential over
        :type spec: object
        :param w_coeffs: `w`-scaled polynomial coefficients to convert into raw coefficients, if `coeffs` isn't given directly
        :type w_coeffs: list[float] | None
        :param coeffs: raw Taylor-series coefficients for the polynomial, used directly if given
        :type coeffs: list[float] | None
        :param g: the coordinate's G-matrix element, used to convert `w_coeffs` into raw coefficients
        :type g: float
        :param re: the equilibrium/reference coordinate value (the polynomial's expansion center)
        :type re: float
        :return: the constructed polynomial potential function
        :rtype: CoordinateFunction
        """
        if w_coeffs is not None:
            sg = np.sqrt(g)
            coeffs = [0] + [np.sign(c) * (np.abs(c) / sg) ** (k + 2) / math.factorial(k + 2) for k, c in enumerate(w_coeffs)]
        return CoordinateFunction.polynomial(
            spec,
            center=re,
            coeffs=coeffs,
            ref=0
        )
    @classmethod
    def get_potential_types(cls):
        """
        **LLM Docstring**

        The mapping from potential-form name to the constructor method that builds that form, used by the potential-generator dispatch.

        :return: the name-to-constructor mapping (`'morse'`, `'poly'`)
        :rtype: dict
        """
        return {
            'morse':cls.get_morse_potential,
            'poly':cls.get_poly_potential
        }
    @classmethod
    def get_potentials_by_attributes(cls):
        """
        **LLM Docstring**

        Attribute-based dispatch mapping used to infer which potential form to build from the keys present in a parameter dict (e.g. `w`/`wx` implies Morse; `coeffs`/`w_coeffs` implies polynomial).

        :return: the attribute-tuple-to-form-name mapping
        :rtype: dict
        """
        return {
            ('w', 'wx'): 'morse',
            ('coeffs',): 'poly',
            ('w_coeffs',): 'poly'
        }

    _pot_dispatch = dev.uninitialized
    default_evaluator_type = 'poly'
    @classmethod
    def potential_generator_dispatch(cls) -> 'dev.OptionsMethodDispatch':
        """
        **LLM Docstring**

        Lazily build (and cache) the `OptionsMethodDispatch` object used to resolve a parameter dict into the right potential-form constructor and its arguments.

        :return: the cached dispatch object
        :rtype: dev.OptionsMethodDispatch
        """
        cls._pot_dispatch = dev.handle_uninitialized(
            cls._pot_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_potential_types,),
            kwargs=dict(
                attributes_map=cls.get_potentials_by_attributes()
            )
        )
        return cls._pot_dispatch

    def get_g(self, which, return_tf=True, order=1):
        """
        **LLM Docstring**

        Compute the G-matrix (effective inverse-mass) element(s) for a given internal-coordinate specification at the molecule's current geometry, via the mass-weighted B-matrix.

        :param which: the internal coordinate(s) to compute the G-matrix for
        :type which: object
        :param return_tf: whether to also return the raw coordinate values and coordinate-transformation Jacobians used
        :type return_tf: bool
        :param order: the Jacobian derivative order to compute
        :type order: int
        :return: the G-matrix (submatrix restricted to `which`), or `(G, (r, tf))` including the raw coordinate values/Jacobians if `return_tf` is set
        :rtype: np.ndarray | tuple
        """
        r, tf = nput.internal_coordinate_tensors(self.mol.coords, which,
                                                 return_inverse=True,
                                                 masses=self.mol.atomic_masses,
                                                 order=order)
        b = np.diag(np.repeat(1 / np.sqrt(self.mol.atomic_masses), 3)) @ r[1]
        G = b.T @ b
        if return_tf:
            return G, (r, tf)
        else:
            return G


    def get_anharmonic_parameters(self,
                                  which,
                                  evaluator=None,
                                  energy_expansion=None,
                                  params_handler=None,
                                  # allow_fd=False,
                                  # fd_options=None,
                                  **opts
                                  ):
        """
        **LLM Docstring**

        Fit 1D analytic potential parameters for each of a set of internal coordinates, by extracting each coordinate's diagonal local derivative expansion (from the full Cartesian potential expansion, re-expressed through the coordinate Jacobian) and fitting it via `params_handler` (defaulting to `get_potential_params`).

        :param which: the internal coordinate(s) to fit potentials for
        :type which: object
        :param evaluator: an explicit energy evaluator to use when computing the potential derivatives
        :type evaluator: object | None
        :param energy_expansion: a precomputed Cartesian potential-derivative expansion to use instead of computing a fresh one
        :type energy_expansion: list[np.ndarray] | None
        :param params_handler: the function used to fit each coordinate's local derivative expansion into potential parameters; defaults to `get_potential_params`
        :type params_handler: callable | None
        :param opts: extra options forwarded to `params_handler`
        :type opts: dict
        :return: the list of fitted parameter dicts, one per coordinate in `which`
        :rtype: list[dict]
        """
        if params_handler is None:
            params_handler = self.get_potential_params

        if nput.is_numeric(which[0]):
            which = [which]

        G, (r, tf) = self.get_g(which, return_tf=True, order=4)
        if energy_expansion is None:
            energy_expansion = self.mol.get_cartesian_potential_derivatives(evaluator=evaluator, order=4)
        derivs = nput.tensor_reexpand(tf, energy_expansion, order=4)

        bond_params = []
        for n, spec in enumerate(which):
            diag_derivs = []
            for k, d in enumerate(derivs):
                idx = (n,) * (k + 1)
                diag_derivs.append(d[idx])
            bond_params.append(
                params_handler(spec, r[0][n], G[n, n], diag_derivs, **opts)
            )

        return bond_params

    def get_1d_potentials(self,
                          which,
                          evaluator=None,
                          energy_expansion=None,
                          params_handler=None,
                          potential_params=None,
                          coordinate_potential_handler=None,
                          **opts
                          ):
        """
        **LLM Docstring**

        Build 1D analytic potential-energy functions for a set of internal coordinates, using any explicitly supplied per-coordinate parameters where given and automatically fitting the rest (via `get_anharmonic_parameters`), then constructing each potential via the resolved potential-form dispatch.

        :param which: the internal coordinate(s) to build potentials for, or a dict mapping coordinates to explicit parameter dicts
        :type which: object | dict
        :param evaluator: an explicit energy evaluator to use when auto-fitting missing parameters
        :type evaluator: object | None
        :param energy_expansion: a precomputed Cartesian potential-derivative expansion to reuse when auto-fitting
        :type energy_expansion: list[np.ndarray] | None
        :param params_handler: the function used to auto-fit missing coordinate parameters; forwarded to `get_anharmonic_parameters`
        :type params_handler: callable | None
        :param potential_params: explicit per-coordinate parameter dicts (list or dict form) to use instead of auto-fitting, for coordinates where available
        :type potential_params: list | dict | None
        :param coordinate_potential_handler: an explicit potential-form dispatch object to use instead of `potential_generator_dispatch()`
        :type coordinate_potential_handler: dev.OptionsMethodDispatch | None
        :param opts: extra options forwarded to `get_anharmonic_parameters` for any auto-fitted coordinates
        :type opts: dict
        :return: the list of constructed 1D potential functions, one per coordinate in `which`
        :rtype: list[CoordinateFunction]
        """
        if dev.is_dict_like(which):
            potential_params = which
            which = list(which.keys())

        if potential_params is not None:
            if dev.is_dict_like(potential_params):
                which = [
                    coordops.canonicalize_internal(w)
                    for w in which
                ]
                potential_params = {
                    coordops.canonicalize_internal(k): v
                    for k, v in potential_params.items()
                }
                potential_params = [
                    potential_params.get(w)
                    for w in which
                ]
            else:
                potential_params = list(potential_params)
            subwhich = [
                w for w, p in enumerate(potential_params)
                if p is None
            ]
            remwhich = [
                w for w, p in enumerate(potential_params)
                if p is not None
            ]
            gwhich = [
                w for w in remwhich
                if potential_params[w].get('g') is None or potential_params[w].get('re') is None
            ]
            if len(gwhich) > 0:
                G, (r, tf) = self.get_g([which[i] for i in gwhich], return_tf=True)
                for n,w in enumerate(gwhich):
                    potential_params[w] = potential_params[w].copy()
                    if potential_params[w].get('g') is None:
                        potential_params[w]['g'] = G[n, n]
                    if potential_params[w].get('re') is None:
                        potential_params[w]['re'] = r[0][n]
            if len(subwhich) > 0:
                subparams = self.get_anharmonic_parameters([which[k] for k in subwhich],
                                                           evaluator=evaluator,
                                                           energy_expansion=energy_expansion,
                                                           params_handler=params_handler,
                                                           **opts
                                                           )
                for k, p in zip(subwhich, subparams):
                    potential_params[k] = p

        else:
            potential_params = self.get_anharmonic_parameters(
                which,
                evaluator=evaluator,
                energy_expansion=energy_expansion,
                params_handler=params_handler,
                **opts
            )

        if coordinate_potential_handler is None:
            coordinate_potential_handler = self.potential_generator_dispatch()
        pots = []
        for spec, params in zip(which, potential_params):
            handler, props = coordinate_potential_handler.resolve(params)
            pots.append(
                handler(spec, **props)
            )

        return pots