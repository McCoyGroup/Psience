"""
Provides a `MolecularTopology` class that encapsulates all of the bond-topology-related
functionality that used to live directly on `Molecule`, plus a `MolecularCoordinateChoice`
class that encapsulates the internal-coordinate selection/labeling functionality.

The topology owns the atom labels, the bond list, the (lazily built) `EdgeGraph`, the
bond-guessing configuration, and the (lazily computed) fragmentation, along with all of the
graph-based algorithms (backbone/segment finding, Z-matrix construction from the bonding
graph, neighborhood/path queries, etc.).

A `MolecularTopology` is self-contained: it needs only the atom labels and bonds, not a
reference to a `Molecule`. Bond guessing is delegated to a caller-supplied `bond_guesser`
callback (the parent `Molecule` hooks this up to go through the `RDMolecule` pass), and the
few routines that need geometry (currently just `get_bond_zmatrix`'s inter-fragment
distance matrix) accept `coords` explicitly.

`MolecularCoordinateChoice` gathers the routines that decide *which* internal coordinates to
use (auto/natural/redundant primitive selection, canonicalization of the many accepted
internal-coordinate spec forms, bond-graph coordinate generation, pruning, and
labeling/mode-labeling). To keep its semantic surface small it does not hold a full `Molecule`
reference: instead it holds a `MolecularEmbedding` (for coordinates and masses) and a
`MolecularTopology` (for atoms, bonds, edge graph, fragmentation and the bond-graph Z-matrix).
The two genuinely molecule-level quantities it can't derive on its own -- the G-matrix and the
normal modes -- are passed in by the caller (the `Molecule`) at the call sites that need them.
The pure spec-construction helpers are class/static methods that take what they need
explicitly.
"""

from __future__ import annotations

import numpy as np

from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
import McUtils.Iterators as itut
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
from McUtils.Coordinerds import PrimitiveCoordinatePicker, RedundantCoordinateGenerator
from McUtils.Graphs import EdgeGraph
from McUtils.ExternalPrograms import RDMolecule

__all__ = [
    "MolecularTopology",
    "MolecularCoordinateChoice"
]

__reload_hook__ = ['.Properties']


class MolecularTopology:
    """
    Stores and manipulates the bonding topology of a molecule.

    All of the connectivity-level data (atoms, bonds, edge graph, fragmentation) and the graph
    algorithms that act on it live here. The class holds no reference to a `Molecule`; anything
    geometry-dependent is passed in explicitly, and bond guessing is delegated to the
    `bond_guesser` callback.
    """

    def __init__(self, atoms, bonds=None, guess_bonds=True, bond_guesser=None):
        """
        :param atoms: the element symbols (or atom labels) for the nodes of the graph
        :type atoms: Iterable[str]
        :param bonds: bond specification for the molecule
        :type bonds: Iterable[Iterable[int]] | None
        :param guess_bonds: whether or not to guess the bonding arrangement when `bonds` is
            `None` and the bonds are actually needed
        :type guess_bonds: bool
        :param bond_guesser: a callable that returns a guessed bond list; called lazily when
            `bonds` is `None` and `guess_bonds` is `True`. Typically hooked up by the parent
            `Molecule` to run the `RDMolecule` bond-guessing pass.
        :type bond_guesser: callable | None
        """
        self.atoms = tuple(atoms)
        self._bonds = bonds
        self._edge_graph = None
        self.guess_bonds = guess_bonds
        self.bond_guesser = bond_guesser
        self._fragment_indices = None
        # `MolecularProperties.edge_graph`/`fragment_indices` only need `.atoms` and `.bonds`,
        # but the latter still reads (and discards) a `.coords` attribute; expose a harmless
        # `None` so a topology can be passed straight through without a `Molecule`.
        self.coords = None

    def __repr__(self):
        return "{}({} atoms, {} bonds)".format(
            type(self).__name__,
            len(self.atoms),
            "?" if self._bonds is None else len(self._bonds)
        )

    #region Bonds
    @property
    def bonds(self):
        """
        Property getter/setter for the bonds. The getter lazily guesses the bonds (via the
        `bond_guesser` callback) if none are set and `self.guess_bonds` is enabled.

        :param b: (setter only) the new bonds
        :type b: list[tuple] | None
        :return: (getter) the bonds, or `None` if unset and bond-guessing is disabled
        :rtype: list[tuple] | None
        """
        if self._bonds is None and self.guess_bonds:
            self._bonds = self.get_guessed_bonds()
        return self._bonds
    @bonds.setter
    def bonds(self, b):
        """
        Property getter/setter for the bonds. The getter lazily guesses the bonds (via the
        `bond_guesser` callback) if none are set and `self.guess_bonds` is enabled.

        :param b: (setter only) the new bonds
        :type b: list[tuple] | None
        :return: (getter) the bonds, or `None` if unset and bond-guessing is disabled
        :rtype: list[tuple] | None
        """
        self._bonds = b
        self._edge_graph = None
        self._fragment_indices = None

    def get_guessed_bonds(self, **opts):
        """
        Guess the bonding arrangement by invoking the caller-supplied `bond_guesser` callback.

        :param opts: extra options forwarded to the `bond_guesser` callback
        :type opts: dict
        :return: the guessed bonds
        :rtype: list[tuple]
        :raises ValueError: if no `bond_guesser` callback was supplied
        """
        if self.bond_guesser is None:
            raise ValueError(
                "{} can't guess bonds without a `bond_guesser` callback".format(
                    type(self).__name__
                )
            )
        return self.bond_guesser(**opts)

    def break_bonds(self, bonds):
        """
        Build a copy of this topology with the specified bonds removed, filtering the bond list
        directly.

        :param bonds: the bonds to remove, each an atom-index pair
        :type bonds: Iterable[tuple]
        :return: the new topology with the given bonds removed
        :rtype: MolecularTopology
        """
        bond_sets = [{b[0], b[1]} for b in bonds]
        return self.modify(
            bonds=[b for b in self.bonds if
                   all(b[0] not in bs or b[1] not in bs for bs in bond_sets)]
        )
    #endregion

    #region Construction helpers
    def modify(self, atoms=None, bonds=None, guess_bonds=None, bond_guesser=None):
        """
        Build a new `MolecularTopology` that copies this one with the given fields overridden.

        :param atoms: replacement atoms, or `None` to keep the current ones
        :type atoms: Iterable[str] | None
        :param bonds: replacement bonds, or `None` to keep the current ones
        :type bonds: Iterable[Iterable[int]] | None
        :param guess_bonds: replacement bond-guessing flag, or `None` to keep the current one
        :type guess_bonds: bool | None
        :param bond_guesser: replacement bond-guessing callback, or `None` to keep the current one
        :type bond_guesser: callable | None
        :return: the new topology
        :rtype: MolecularTopology
        """
        return type(self)(
            self.atoms if atoms is None else atoms,
            bonds=self._bonds if bonds is None else bonds,
            guess_bonds=self.guess_bonds if guess_bonds is None else guess_bonds,
            bond_guesser=self.bond_guesser if bond_guesser is None else bond_guesser
        )

    def take_subgraph(self, pos):
        """
        Build a new `MolecularTopology` containing only the graph nodes at the given positions
        and the bonds internal to that subset, remapping bonds to the new (sub)indexing and
        dropping any bonds that reference atoms outside the subset.

        :param pos: the atom indices to keep, in the desired order for the subgraph
        :type pos: Iterable[int]
        :return: the constructed subgraph topology
        :rtype: MolecularTopology
        """
        pos = list(pos)
        atoms = [self.atoms[i] for i in pos]
        remap = {p: i for i, p in enumerate(pos)}
        if self._bonds is not None:
            bonds = [
                (
                    [remap[b[0]], remap[b[1]]]
                        if len(b) == 2 else
                    [remap[b[0]], remap[b[1]], b[2]]
                )
                for b in self._bonds
                if b[0] in remap and b[1] in remap
            ]
        else:
            bonds = None
        return type(self)(
            atoms,
            bonds=bonds,
            guess_bonds=self.guess_bonds,
            bond_guesser=None
        )
    #endregion

    #region Graph
    @property
    def edge_graph(self) -> EdgeGraph:
        """
        The (cached) `EdgeGraph` representation of the bonding structure, built lazily via
        `MolecularProperties.edge_graph` (which needs only the atoms and bonds this topology
        already exposes).

        :return: the edge graph
        :rtype: EdgeGraph
        """
        if self._edge_graph is None:
            from .Properties import MolecularProperties
            self._edge_graph = MolecularProperties.edge_graph(self)
        return self._edge_graph

    def find_path(self, atom1, atom2):
        """
        Find a path between two atoms through the bonding graph.

        :param atom1: the starting atom index
        :type atom1: int
        :param atom2: the ending atom index
        :type atom2: int
        :return: the path between the two atoms
        :rtype: list[int]
        """
        return self.edge_graph.get_path(atom1, atom2)

    def neighborhood(self, loc, size=1, heavy_only=False):
        """
        Find the atoms within a given graph-distance of a location in the bonding graph.

        :param loc: the atom index to center the neighborhood on
        :type loc: int
        :param size: the neighborhood radius (in bond-graph steps)
        :type size: int
        :return: the neighboring atom indices
        :rtype: tuple[int]
        """
        nb = tuple(l for l in self.edge_graph.neighbor_iterator(loc, num=size))
        if heavy_only:
            nb = tuple(n for n in nb if self.atoms[n] not in ("H", "D", "T"))
        return nb
    #endregion

    #region Fragments
    @property
    def fragment_indices(self):
        """
        The (cached) grouping of atom indices into connected molecular fragments, computed
        lazily via `MolecularProperties.fragment_indices` (which needs only the atoms and bonds
        this topology already exposes).

        :return: the list of per-fragment atom-index groups
        :rtype: list[np.ndarray]
        """
        if self._fragment_indices is None:
            from .Properties import MolecularProperties
            self._fragment_indices = MolecularProperties.fragment_indices(self)
        return self._fragment_indices

    @property
    def fragments(self):
        """
        The topology split into its connected fragments (as separate `MolecularTopology`
        objects), by taking the subgraph for each connected component.

        :return: the list of fragment topologies
        :rtype: list[MolecularTopology]
        """
        return [self.take_subgraph(c) for c in self.fragment_indices]
    #endregion

    #region Backbone / Z-matrix
    def find_heavy_atom_backbone(self, root=None):
        """
        Find the longest chain of atoms in the bonding graph (the heavy-atom backbone), via
        `edge_graph.find_longest_chain`.

        :param root: an atom index to force as one end of the chain
        :type root: int | None
        :return: the backbone atom indices, in chain order
        :rtype: list[int]
        """
        return self.edge_graph.find_longest_chain(root=root)

    def find_backbone_segments(self, root=None, initial_backbone=None):
        """
        Split the bonding graph into backbone-connected segments, via
        `edge_graph.segment_by_chains`.

        :param root: an atom index to anchor the segmentation at
        :type root: int | None
        :param initial_backbone: an initial backbone chain to build the segmentation around
        :type initial_backbone: Iterable[int] | None
        :return: the resulting chain segments
        :rtype: list
        """
        return self.edge_graph.segment_by_chains(root=root, backbone=initial_backbone)

    def get_backbone_zmatrix(self, root=None,
                             segments=None,
                             return_remainder=False,
                             return_segments=False,
                             required_coordinates=None,
                             isolated_coordinates=None,
                             root_coordinates=None,
                             initial_backbone=None,
                             validate=True
                             ):
        """
        Build a Z-matrix for a (typically single-fragment) molecule by first segmenting its
        bonding graph into backbone chains (via `find_backbone_segments`, validating there are
        no duplicate atoms across segments), constructing the base Z-matrix graph from the
        bonds and segments, filling in any bonds missing from the initial graph, and (if
        requested) enforcing required/isolated/root coordinate constraints.

        :param root: an atom index to anchor the backbone segmentation at
        :type root: int | None
        :param segments: precomputed backbone segments to use instead of calling
            `find_backbone_segments`
        :type segments: list | None
        :param return_remainder: whether to also return the bonds that had to be added beyond
            the base backbone graph
        :type return_remainder: bool
        :param return_segments: whether to also return the backbone segments used
        :type return_segments: bool
        :param required_coordinates: internal coordinates that must appear in the resulting
            Z-matrix
        :type required_coordinates: Iterable | None
        :param isolated_coordinates: coordinates that must be built in isolation from others
        :type isolated_coordinates: Iterable | None
        :param root_coordinates: coordinates to anchor at the root of the Z-matrix
        :type root_coordinates: Iterable | None
        :param initial_backbone: an initial backbone chain to seed the segmentation with
        :type initial_backbone: Iterable[int] | None
        :param validate: whether to validate the Z-matrix construction at each step (duplicate
            atoms, valid additions)
        :type validate: bool
        :return: the Z-matrix, or a tuple additionally including the segments and/or remainder
            bonds depending on the `return_*` flags
        :rtype: np.ndarray | tuple
        :raises ValueError: if `validate` is set and duplicate atoms are found across backbone
            segments
        """
        if segments is None:
            segments = self.find_backbone_segments(root=root, initial_backbone=initial_backbone)
            if validate:
                flat_frags = list(itut.flatten(segments))
                frag_counts = itut.counts(flat_frags)
                bad_frags = {k: v for k, v in frag_counts.items() if v > 1}
                if len(bad_frags) > 0:
                    raise ValueError(f"diplicate atoms {list(bad_frags.keys())} encountered in {segments}")

        bond_list = [b[:2] for b in self.bonds]
        base_graph = coordops.bond_graph_zmatrix(
            bond_list,
            segments,
            validate_additions=validate,
            required_coordinates=required_coordinates,
            isolated_coordinates=isolated_coordinates,
            root_coordinates=root_coordinates,
            enforce_requirements=False
        )
        zmat, new_bonds = coordops.add_missing_zmatrix_bonds(
            base_graph,
            bond_list,
            validate_additions=validate
        )
        if (
                required_coordinates is not None
                or isolated_coordinates is not None
                or root_coordinates is not None
        ):
            zmat = coordops.enforce_required_zmatrix_coordinates(zmat,
                                                                 required_coordinates,
                                                                 isolated_coordinates=isolated_coordinates,
                                                                 root_coordinates=root_coordinates,
                                                                 validate=validate)

        if return_segments or return_remainder:
            res = (zmat,)
            if return_segments:
                res = res + (segments,)
            if return_remainder:
                res = res + (new_bonds,)

            return res
        else:
            return zmat

    def get_canonical_zmatrix(self, ordering=None, validate=True):
        """
        Build a canonical Z-matrix ordering for the molecule from a canonical fragmentation of
        the bonding graph.

        :param ordering: the atom ordering to use as the basis for canonicalization; defaults
            to the natural `0..N` ordering
        :type ordering: np.ndarray | None
        :param validate: whether to validate each Z-matrix addition
        :type validate: bool
        :return: the canonical Z-matrix
        :rtype: np.ndarray
        """
        if ordering is None: ordering = np.arange(len(self.atoms))
        frags = self.edge_graph.get_canonical_fragments(ordering)
        return coordops.canonical_fragment_zmatrix(frags, validate_additions=validate)

    @classmethod
    def _filter_coordinates_by_fragments(cls, inds, frags, required_coordinates):
        """
        Split a list of required internal-coordinate specs into those fully contained within a
        single fragment (reindexed to that fragment's local atom numbering) versus those that
        span multiple fragments and must be merged in afterward.

        :param inds: the atom-index lists defining each fragment
        :type inds: list[Iterable[int]]
        :param frags: the corresponding sub-topologies for each fragment (unused directly in the
            body but kept for interface/length consistency)
        :type frags: list[MolecularTopology]
        :param required_coordinates: the coordinate specs to classify, each a tuple of atom
            indices in the full-molecule numbering
        :type required_coordinates: Iterable[tuple] | None
        :return: `(merge_coordinates, fragment_requireds)` -- the coordinates that couldn't be
            assigned to a single fragment, and a per-fragment list of the (locally reindexed)
            coordinates assigned to each fragment
        :rtype: tuple[list, list]
        """
        merge_coordinates = []
        if required_coordinates is not None:
            fragment_requireds = [[] for _ in range(len(frags))]
            for c in required_coordinates:
                for i, f in enumerate(inds):
                    ff = list(f)
                    sub = []
                    for j in c:
                        try:
                            x = ff.index(j)
                        except ValueError:
                            break
                        else:
                            sub.append(x)
                    if len(sub) == len(c):
                        fragment_requireds[i].append(tuple(sub))
                        break
                else:
                    merge_coordinates.append(c)
        else:
            fragment_requireds = [None] * len(inds)
        return merge_coordinates, fragment_requireds

    def get_bond_zmatrix(self,
                         coords=None,
                         fragments=None,
                         segments=None,
                         root=None,
                         required_coordinates=None,
                         isolated_coordinates=None,
                         root_coordinates=None,
                         attachment_points=None,
                         check_attachment_points=True,
                         validate=True,
                         for_fragment=None,
                         fragment_ordering=None,
                         connect_fragments=True,
                         initial_backbone=None,
                         distance_matrix=None
                         ):
        """
        Build a full Z-matrix for the molecule from its bonding graph, handling the
        single-fragment case via `get_backbone_zmatrix` directly and the multi-fragment case by
        building a per-fragment Z-matrix for each fragment (optionally
        reordering/rooting/filtering required coordinates per fragment) and then splicing them
        together into one connected Z-matrix (via `coordops.complex_zmatrix`) unless
        `connect_fragments` is `False`; can also be restricted to build the Z-matrix for just
        one fragment (`for_fragment`), in which case it recurses on the corresponding
        subgraph and reindexes the result back to the full atom numbering.

        The inter-fragment connection step needs a distance matrix; if `distance_matrix` is not
        supplied it is computed from `coords` when those are given, and left as `None`
        otherwise.

        :param coords: the Cartesian coordinates, used only to build the inter-fragment distance
            matrix when connecting fragments; may be `None` if no distance information is needed
        :type coords: np.ndarray | None
        :param fragments: explicit fragment atom-index groups to use instead of
            `self.fragment_indices`
        :type fragments: list[Iterable[int]] | None
        :param segments: precomputed backbone segments for the single-fragment case
        :type segments: list | None
        :param root: root atom(s) to anchor the Z-matrix construction at (per fragment, in the
            multi-fragment case)
        :type root: int | Iterable | None
        :param required_coordinates: internal coordinates that must appear in the resulting
            Z-matrix
        :type required_coordinates: Iterable | None
        :param isolated_coordinates: coordinates that must be built in isolation from others
        :type isolated_coordinates: Iterable | None
        :param root_coordinates: coordinates to anchor at the root of the Z-matrix
        :type root_coordinates: Iterable | None
        :param attachment_points: explicit inter-fragment attachment points to use when
            connecting fragments
        :type attachment_points: dict | Iterable | None
        :param check_attachment_points: whether to validate the attachment points used when
            connecting fragments
        :type check_attachment_points: bool
        :param validate: whether to validate each Z-matrix addition
        :type validate: bool
        :param for_fragment: restrict the Z-matrix construction to just this fragment (an index
            into `self.fragment_indices`, or an explicit list of atom indices), returning the
            result reindexed to the full molecule
        :type for_fragment: int | Iterable[int] | None
        :param fragment_ordering: explicit ordering to apply to the fragments before connecting
            them
        :type fragment_ordering: Iterable[int] | None
        :param connect_fragments: whether to splice the per-fragment Z-matrices into one
            connected Z-matrix, or return them as a list of separate Z-matrices
        :type connect_fragments: bool
        :param initial_backbone: an initial backbone chain to seed the segmentation of each
            fragment with
        :type initial_backbone: Iterable[int] | None
        :param distance_matrix: an explicit inter-atom distance matrix to use when connecting
            fragments; computed from `coords` if not given (and left `None` if `coords` is
            `None`)
        :type distance_matrix: np.ndarray | None
        :return: the (connected) Z-matrix, or a list of per-fragment Z-matrices if
            `connect_fragments` is `False`
        :rtype: np.ndarray | list[np.ndarray]
        """
        if for_fragment is not None:
            if nput.is_int(for_fragment):
                for_fragment = self.fragment_indices[for_fragment]
            if attachment_points is not None:
                frag_map = dict(zip(for_fragment, np.arange(len(for_fragment))))
                if hasattr(attachment_points, 'items'):
                    attachment_points = {
                        frag_map[i]: (frag_map[k] if nput.is_numeric(k) else tuple(frag_map[kk] for kk in k))
                        for i,k in attachment_points.items()
                    }
                else:
                    attachment_points = [
                        frag_map[i] for i in attachment_points
                    ]

            if initial_backbone is not None:
                ff = list(for_fragment)
                initial_backbone = [ff.index(i) for i in initial_backbone]
            sub_coords = None if coords is None else np.asanyarray(coords)[..., for_fragment, :]
            base_ints = self.take_subgraph(for_fragment).get_bond_zmatrix(
                coords=sub_coords,
                fragments=fragments, segments=segments, root=root,
                attachment_points=attachment_points,
                check_attachment_points=check_attachment_points,
                fragment_ordering=fragment_ordering,
                required_coordinates=required_coordinates,
                isolated_coordinates=isolated_coordinates,
                root_coordinates=root_coordinates,
                initial_backbone=initial_backbone,
                validate=validate
            )
            zm = coordops.reindex_zmatrix(base_ints, for_fragment)
            return np.asarray(zm)
        else:
            no_frag = fragments is None
            if no_frag:
                fragments = self.fragment_indices

            if len(fragments) == 1:
                if segments is not None and len(segments) == 1:
                    segments = segments[0]
                zm = self.get_backbone_zmatrix(
                    root=root, segments=segments,
                    required_coordinates=required_coordinates,
                    isolated_coordinates=isolated_coordinates,
                    root_coordinates=root_coordinates,
                    validate=validate,
                    initial_backbone=initial_backbone
                )
                zm = np.asarray(zm)
                if connect_fragments:
                    return zm
                else:
                    return [zm]
            else:
                inds = fragments
                if no_frag:
                    if fragment_ordering is None:
                        fragment_ordering = np.argsort([-len(x) for x in inds])
                    inds = [inds[i] for i in fragment_ordering]
                if root is not None:
                    if nput.is_numeric(root):
                        inds = list(sorted(inds, key=lambda x:root not in x))
                    else:
                        inds = list(
                            sorted(inds,
                                   key=lambda x:sum(i if r is not None and r in x else len(inds) for i,r in enumerate(root))
                                   )
                        )

                sort_attch = isinstance(attachment_points, dict)
                if sort_attch:
                    check_attachment_points = False
                    inds, attachment_points = coordops.sort_complex_attachment_points(
                        inds,
                        attachment_points
                    )

                frags = [self.take_subgraph(ix) for ix in inds]
                if root is None and sort_attch:
                    root = [ix[0] for ix in inds]
                if root is None:
                    root = [root]


                if initial_backbone is not None:
                    initial_backbones = []
                    for frag in inds:
                        ff = list(frag)
                        sub = []
                        for i in initial_backbone:
                            try:
                                x = ff.index(i)
                            except ValueError:
                                continue
                            else:
                                sub.append(x)
                        if len(sub) > 0:
                            initial_backbones.append(sub)
                        else:
                            initial_backbones.append(None)
                else:
                    initial_backbones = [None] * len(inds)

                root = list(root) + [None] * (len(inds) - len(root))
                merge_coordinates, fragment_requireds = self._filter_coordinates_by_fragments(inds, frags, required_coordinates)
                merge_isolated, fragment_isolated = self._filter_coordinates_by_fragments(inds, frags, isolated_coordinates)
                merge_root, fragment_root = self._filter_coordinates_by_fragments(inds, frags, root_coordinates)
                if len(merge_coordinates) == 0:
                    merge_coordinates = None
                if len(merge_isolated) == 0:
                    merge_isolated = None
                if len(merge_root) == 0:
                    merge_root = None
                zmats = [
                    f.get_backbone_zmatrix(root=r, initial_backbone=bb,
                                           required_coordinates=rq,
                                           isolated_coordinates=iso,
                                           root_coordinates=rot
                                           )
                    for r,f,bb,rq,iso,rot in zip(root, frags, initial_backbones,
                                                 fragment_requireds, fragment_isolated, fragment_root)
                ]

                if connect_fragments:

                    # inds = [inds[i] for i in ordering]
                    # zmats = [zmats[i] for i in ordering]

                    if distance_matrix is None and coords is not None:
                        distance_matrix = nput.distance_matrix(coords)
                    dm = distance_matrix
                    if dm is not None:
                        dm = np.array(dm, copy=True)
                        h_pos = [i for i, a in enumerate(self.atoms) if a in {'H', 'D'}]
                        dm[:, h_pos] = 1e8
                        dm[h_pos, :] = 1e8

                    zm = coordops.complex_zmatrix(
                        [b[:2] for b in self.bonds],
                        inds,
                        zmats,
                        distance_matrix=dm,
                        attachment_points=attachment_points,
                        check_attachment_points=check_attachment_points,
                        required_coordinates=required_coordinates,
                        isolated_coordinates=isolated_coordinates,
                        root_coordinates=root_coordinates,
                        validate_additions=validate
                    )
                    return np.asarray(zm)
                else:
                    shift_mats = []
                    offsets = 0
                    for zmat in zmats:
                        shift_mats.append([
                            [(z + offsets) if z >= 0 else z for z in zm]
                            for zm in zmat
                        ])
                        offsets += len(zmat)
                    return [np.asarray(zm) for zm in shift_mats]
    #endregion


class MolecularCoordinateChoice:
    """
    Encapsulates the internal-coordinate *choice* functionality: selecting primitive internal
    coordinates from the bonding graph (auto/natural/redundant), canonicalizing the many
    accepted internal-coordinate specification forms, generating and pruning bond-graph
    coordinates, and labeling coordinates / normal modes.

    The pure spec-construction helpers (`_generate_auto_spec`, `_generate_stretch_spec`,
    `_auto_auto_spec`, `_auto_spec`, `_stretch_spec`, `_check_label`, `get_coordinate_filer`)
    are class/static methods that take everything they need explicitly.

    The instance-level entry points (`canonicalize_internals`, `prep_internal_spec`,
    `get_bond_graph_internals`, `prune_internals`, `get_labeled_internals`, `get_mode_labels`)
    need geometry (coordinates, masses) and connectivity (atoms, bonds, edge graph,
    fragmentation, bond-graph Z-matrix). Rather than holding a full `Molecule` reference, they
    hold just a `MolecularEmbedding` (for coordinates/masses and Z-matrix geometry) and a
    `MolecularTopology` (for atoms/bonds/graph/fragments) -- keeping the semantic surface small.
    The two genuinely molecule-level quantities that this object can't derive on its own --
    the G-matrix and the normal modes -- are *passed in* by the caller (the `Molecule`) at the
    call sites that need them (`prune_internals`'s B-matrix generator and `get_mode_labels`).
    """

    default_coordinate_pruning = 'graph'

    def __init__(self, embedding, topology):
        """
        :param embedding: the coordinate embedding supplying Cartesian coordinates and
            (atomic-unit) masses
        :type embedding: MolecularEmbedding
        :param topology: the bonding topology supplying atoms, bonds, edge graph, fragmentation
            and the bond-graph Z-matrix
        :type topology: MolecularTopology
        """
        self.embedding = embedding
        self.topology = topology

    def set_embedding(self, embedding):
        """
        Rebind this coordinate-choice helper to a (new) coordinate embedding.

        :param embedding: the coordinate embedding this helper should use
        :type embedding: MolecularEmbedding
        :return: None
        :rtype: None
        """
        self.embedding = embedding

    def set_topology(self, topology):
        """
        Rebind this coordinate-choice helper to a (new) bonding topology.

        :param topology: the bonding topology this helper should use
        :type topology: MolecularTopology
        :return: None
        :rtype: None
        """
        self.topology = topology

    def _take_subchoice(self, fragment):
        """
        Build a `MolecularCoordinateChoice` restricted to a fragment (subset of atom indices),
        by slicing this object's embedding down to those atoms and taking the corresponding
        subgraph of the topology.

        :param fragment: the atom indices defining the fragment
        :type fragment: Iterable[int]
        :param fragment: the atom indices to keep
        :return: the sub-fragment coordinate-choice helper
        :rtype: MolecularCoordinateChoice
        """
        fragment = list(fragment)
        emb = self.embedding
        sub_coords = np.asanyarray(emb.coords)[..., fragment, :]
        sub_masses = np.asanyarray(emb.masses)[fragment,]
        sub_embedding = type(emb)(sub_masses, sub_coords, None)
        sub_topology = self.topology.take_subgraph(fragment)
        return type(self)(sub_embedding, sub_topology)

    #region Primitive spec generation
    @classmethod
    def _generate_auto_spec(cls, atoms, bonds, base_coords=None, **opts):
        """
        Automatically pick a set of primitive internal coordinates (bond stretches, angles,
        dihedrals) from the bonding graph, via `PrimitiveCoordinatePicker`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param bonds: the bonds to build coordinates from
        :type bonds: Iterable[Iterable[int]]
        :param base_coords: coordinates to seed/prioritize the picker with
        :type base_coords: Iterable | None
        :param opts: extra options forwarded to `PrimitiveCoordinatePicker`
        :type opts: dict
        :return: the picked coordinate specs
        :rtype: list
        """
        return PrimitiveCoordinatePicker(
            atoms,
            [b[:2] for b in bonds],
            base_coords=base_coords,
            **opts
        ).coords

    @classmethod
    def _generate_stretch_spec(cls, atoms, bonds, **opts):
        """
        Build a "natural"-coordinate specification consisting only of the bond-stretch
        coordinates implied by `bonds`.

        :param atoms: the atom labels (unused directly, kept for interface consistency with
            `_generate_auto_spec`)
        :type atoms: Iterable[str]
        :param bonds: the bonds to build stretch coordinates from
        :type bonds: Iterable[Iterable[int]]
        :param opts: extra options, unused
        :type opts: dict
        :return: the list of bond-stretch coordinate specs
        :rtype: list
        """
        return sum(coordops.get_stretch_coordinate_system([tuple(s[:2]) for s in bonds]), [])

    @classmethod
    def _auto_auto_spec(cls, spec_generator, atoms, coords, bonds, redundant=False, base_coordinates=None,
                        masses=None,
                        untransformed_coordinates=None,
                        prune_coordinates=True,
                        pruning_options=None,
                        formal_charges=None,
                        **opts):
        """
        Shared driver behind `_auto_spec`/`_stretch_spec`: guesses bonds if not given,
        optionally sets up a redundant-coordinate specification (folding in any
        `untransformed_coordinates`), generates the primitive coordinate specs via
        `spec_generator`, and (optionally) prunes them down to a well-conditioned,
        non-redundant subset via `RedundantCoordinateGenerator.prune_coordinate_specs`.

        :param spec_generator: the coordinate-generating function to use (`_generate_auto_spec`
            or `_generate_stretch_spec`)
        :type spec_generator: callable
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use; guessed via RDKit if `None`
        :type bonds: Iterable[Iterable[int]] | None
        :param redundant: whether to build a redundant coordinate specification
        :type redundant: bool
        :param base_coordinates: seed coordinates to prioritize/include
        :type base_coordinates: Iterable | None
        :param masses: atomic masses, used for pruning; computed from `atoms` if not given
        :type masses: np.ndarray | None
        :param untransformed_coordinates: coordinates that should remain untransformed under the
            redundant transformation
        :type untransformed_coordinates: Iterable | None
        :param prune_coordinates: whether to prune the generated coordinate specs down to a
            well-conditioned subset
        :type prune_coordinates: bool
        :param pruning_options: extra options forwarded to
            `RedundantCoordinateGenerator.prune_coordinate_specs`
        :type pruning_options: dict | None
        :param formal_charges: formal charges used when guessing bonds
        :type formal_charges: Iterable | None
        :param opts: extra options forwarded to `spec_generator`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict (with `'specs'`, and
            `'redundant'`/`'untransformed_coordinates'` if applicable)
        :rtype: dict
        """
        base_coords = base_coordinates
        if bonds is None:
            bonds = RDMolecule.from_coords(
                                           atoms,
                                           coords * UnitsData.convert("BohrRadius", "Angstroms"),
                                           bonds,
                                           formal_charges=formal_charges,
                                           guess_bonds=True
                                           ).bonds
        if redundant and untransformed_coordinates is None:
            untransformed_coordinates = base_coords
            base_coords = None
        if redundant and untransformed_coordinates is not None:
            untransformed_coordinates = PrimitiveCoordinatePicker.prep_unique_coords(untransformed_coordinates)
            base_coords = untransformed_coordinates + ([] if base_coords is None else list(base_coords))
        if base_coords is not None:
            base_coords = PrimitiveCoordinatePicker.prep_unique_coords(base_coords)
        specs = spec_generator(atoms, bonds, base_coords=base_coords, **opts)
        if prune_coordinates:
            if pruning_options is None:
                pruning_options = {}
            expansion = nput.internal_coordinate_tensors(coords, specs, order=1)[1:]
            if masses is None:
                ats = [AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms]
                masses = np.array([a["Mass"] for a in ats])
            prune_pos = RedundantCoordinateGenerator.prune_coordinate_specs(
                expansion,
                masses=masses,
                untransformed_coordinates=np.arange(len(base_coords)) if base_coords is not None else None,
                **pruning_options
            )
            specs = [specs[i] for i in prune_pos]
        spec = {'specs':specs}
        if redundant:
            spec['redundant'] = True
            if untransformed_coordinates is not None:
                spec['untransformed_coordinates'] = np.arange(len(untransformed_coordinates))
        return spec

    @classmethod
    def _auto_spec(cls, atoms, coords, bonds, **opts):
        """
        Build an automatically-chosen internal-coordinate specification from the bonding graph,
        via `_auto_auto_spec` with `_generate_auto_spec`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use
        :type bonds: Iterable[Iterable[int]] | None
        :param opts: extra options forwarded to `_auto_auto_spec`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict
        :rtype: dict
        """
        return cls._auto_auto_spec(cls._generate_auto_spec, atoms, coords, bonds, **opts)

    @classmethod
    def _stretch_spec(cls, atoms, coords, bonds, **opts):
        """
        Build a "natural"/stretch-only internal-coordinate specification from the bonding graph,
        via `_auto_auto_spec` with `_generate_stretch_spec`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use
        :type bonds: Iterable[Iterable[int]] | None
        :param opts: extra options forwarded to `_auto_auto_spec`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict
        :rtype: dict
        """
        return cls._auto_auto_spec(cls._generate_stretch_spec, atoms, coords, bonds, **opts)
    #endregion

    #region Canonicalization
    def canonicalize_internals(self, spec, atoms, coords, bonds, relocalize=True, masses=None):
        """
        Normalize the many accepted forms of an internal-coordinate specification (the strings
        `'auto'`/`'zmatrix'`, a dict with `'primitives'`/`'specs'`/`'zmatrix'` keys where
        `'specs'` may itself be `'auto'`/`'natural'`, a bare Z-matrix-like array, or a bare list
        of primitive specs) down into the canonical dict form expected by `MolecularEmbedding`,
        recursively re-dispatching as needed.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: str | dict | Iterable | None
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use when auto-generating coordinates
        :type bonds: Iterable[Iterable[int]] | None
        :param relocalize: whether redundant coordinates should be relocalized by default
        :type relocalize: bool
        :param masses: atomic masses, forwarded to the auto-generation routines
        :type masses: np.ndarray | None
        :return: the canonicalized specification
        :rtype: dict | None
        :raises ValueError: if `spec` is a string that isn't recognized
            (`'auto'`/`'zmatrix'`/`'natural'`)
        """
        if isinstance(spec, str) and spec.lower() == 'auto':
            spec = {
                'primitives': 'auto'
            }
        elif dev.str_is(spec, 'zmatrix'):
            spec = self.topology.get_bond_zmatrix(coords=self.embedding.coords)

        if isinstance(spec, str):
            # if spec.lower() == 'auto':
            #     spec = cls._auto_spec(atoms, coords, bonds)
            # else:
            raise ValueError(f"can't understand internal spec '{spec}'")
        elif isinstance(spec, dict):
            if 'zmatrix' in spec: return spec
            prims = spec.pop('primitives', None)
            if prims is not None:
                spec = spec.copy()
                spec['specs'] = prims
                spec['redundant'] = True
            subspec = spec.get('specs', '')
            if isinstance(subspec, str):
                if subspec.lower() == 'auto':
                    opts = spec.copy()
                    del opts['specs']
                    if 'relocalize' in opts:
                        relocalize = spec.get('relocalize', relocalize)
                        del opts['relocalize']
                    spec = self._auto_spec(atoms, coords, bonds, masses=masses, **opts)
                elif subspec.lower() == 'natural':
                    opts = spec.copy()
                    del opts['specs']
                    if 'relocalize' in opts:
                        relocalize = spec.get('relocalize', relocalize)
                        del opts['relocalize']
                    spec = self._stretch_spec(atoms, coords, bonds, masses=masses, **opts)
                else:
                    raise ValueError(f"can't understand internal spec '{spec}'")
            else:
                untransformed_coordinates = spec.get('untransformed_coordinates')
                if untransformed_coordinates is not None:
                    if not nput.is_int(untransformed_coordinates[0]):
                        prims = spec.get('specs')
                        untransformed_coordinates = [
                            prims.index(u)
                            for u in untransformed_coordinates
                        ]
                    spec['untransformed_coordinates'] = untransformed_coordinates
            if spec.get('redundant'):
                spec['relocalize'] = spec.get('relocalize', relocalize)
        elif not isinstance(spec, dict) and spec is not None:
            if all(not isinstance(x, dict) and len(x) == 4 for x in spec):
                spec = {'zmatrix': spec}
            else:
                spec = {'primitives':spec}
            spec = self.canonicalize_internals(spec, atoms, coords, bonds, relocalize=relocalize, masses=masses)
        return spec

    def prep_internal_spec(self, spec, relocalize=True, masses=None):
        """
        Canonicalize an internal-coordinate specification against this object's own atoms,
        coordinates, and bonds (from the topology and embedding), via `canonicalize_internals`.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: Any
        :param relocalize: whether redundant coordinates should be relocalized by default
        :type relocalize: bool
        :param masses: atomic masses to use instead of this object's own
        :type masses: np.ndarray | None
        :return: the canonicalized specification
        :rtype: dict | None
        """
        return self.canonicalize_internals(
            spec,
            self.topology.atoms,
            self.embedding.coords,
            self.topology._bonds,
            relocalize=relocalize,
            masses=masses
        )
    #endregion

    #region Labeling helpers
    @classmethod
    def _check_label(cls, label,
                     allowed_coordinate_types=None,
                     excluded_coordinate_types=None,
                     allowed_ring_types=None,
                     excluded_ring_types=None,
                     allowed_group_types=None,
                     excluded_group_types=None,
                     ):
        """
        Test whether a coordinate label passes a set of allow/exclude filters on its atom types,
        ring membership, and functional-group membership.

        :param label: the coordinate label to test (exposing `.atoms`, `.ring`, `.group`
            attributes)
        :type label: Any
        :param allowed_coordinate_types: if given, `label.atoms` must be among these to pass
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: if given, `label.atoms` must not be among these to pass
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: if given, `label.ring` must be among these to pass
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: if given, `label.ring` must not be among these to pass
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: if given, `label.group` must be among these to pass
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: if given, `label.group` must not be among these to pass
        :type excluded_group_types: Iterable | None
        :return: whether the label passes every specified filter
        :rtype: bool
        """
        if allowed_coordinate_types is not None:
            if label.atoms not in allowed_coordinate_types: return False
        if excluded_coordinate_types is not None:
            if label.atoms in excluded_coordinate_types: return False
        if allowed_ring_types is not None:
            if label.ring not in allowed_ring_types: return False
        if excluded_ring_types is not None:
            if label.ring in excluded_ring_types: return False
        if allowed_group_types is not None:
            if label.group not in allowed_group_types: return False
        if excluded_group_types is not None:
            if label.group in excluded_group_types: return False
        return True

    @classmethod
    def get_coordinate_filer(cls,
                             allowed_coordinate_types=None,
                             excluded_coordinate_types=None,
                             allowed_ring_types=None,
                             excluded_ring_types=None,
                             allowed_group_types=None,
                             excluded_group_types=None
                             ):
        """
        Build a filter function (closing over the given allow/exclude criteria) that, given a
        dict of coordinate-to-label mappings, returns only the entries whose label passes
        `_check_label`.

        :param allowed_coordinate_types: forwarded to `_check_label`
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: forwarded to `_check_label`
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: forwarded to `_check_label`
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: forwarded to `_check_label`
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: forwarded to `_check_label`
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: forwarded to `_check_label`
        :type excluded_group_types: Iterable | None
        :return: the constructed coordinate-filtering function
        :rtype: callable
        """
        def coordinate_filter(coords):
            """
            Filter a dict of coordinate-to-label mappings down to just the entries whose label
            satisfies the enclosing allow/exclude criteria (via `_check_label`).

            :param coords: the coordinate-to-label mapping to filter
            :type coords: dict
            :return: the filtered mapping
            :rtype: dict
            """
            return {
                c: l
                for c, l in coords.items()
                if cls._check_label(l,
                                    allowed_coordinate_types=allowed_coordinate_types,
                                    excluded_coordinate_types=excluded_coordinate_types,
                                    allowed_ring_types=allowed_ring_types,
                                    excluded_ring_types=excluded_ring_types,
                                    allowed_group_types=allowed_group_types,
                                    excluded_group_types=excluded_group_types
                                    )
            }

        return coordinate_filter
    #endregion

    #region Bond-graph internals
    def get_bond_graph_internals(self,
                                 include_stretches=True,
                                 include_bends=True,
                                 include_dihedrals=True,
                                 include_fragments=True,
                                 pruning=None,
                                 fragment=None,
                                 base_internals=None,
                                 use_distance_matrix=True,
                                 concatenate=True,
                                 gmatrix=None
                                 ):
        """
        Build a set of internal coordinates (bond stretches, bends, dihedrals, and/or
        inter-fragment coordinates) directly from the bonding graph, optionally restricted to a
        single fragment (recursively, with the result permuted back into the full atom indexing)
        and/or pruned down to a well-conditioned subset.

        :param include_stretches: whether to include bond-stretch coordinates
        :type include_stretches: bool
        :param include_bends: whether to include bond-angle coordinates
        :type include_bends: bool
        :param include_dihedrals: whether to include dihedral-angle coordinates
        :type include_dihedrals: bool
        :param include_fragments: whether to include coordinates connecting separate molecular
            fragments
        :type include_fragments: bool
        :param pruning: whether/how to prune the resulting coordinates (`True` for the default
            method, or an explicit method spec), forwarded to `prune_internals`
        :type pruning: bool | str | dict | None
        :param fragment: restrict to a single fragment, given as a fragment index or an explicit
            list of atom indices
        :type fragment: int | Iterable[int] | None
        :param base_internals: accepted and forwarded when recursing on a fragment, but not
            otherwise used directly in this method's own body
        :type base_internals: Any | None
        :param use_distance_matrix: whether to precompute a distance matrix for the
            fragment-coordinate generation
        :type use_distance_matrix: bool
        :param concatenate: whether to concatenate the different coordinate categories
            (stretches/bends/dihedrals/fragments) into a single list, or return them as separate
            groups
        :type concatenate: bool
        :param gmatrix: the (fractional-power) G-matrix to hand to `prune_internals` if
            B-matrix pruning is requested; must be supplied by the caller (the `Molecule`) when
            such pruning is used, since this object can't compute a G-matrix on its own
        :type gmatrix: np.ndarray | None
        :return: the generated internal coordinates, as a single concatenated list or a list of
            category groups depending on `concatenate`
        :rtype: list
        :raises ValueError: if `pruning` is requested while `concatenate` is `False`
        """
        topology = self.topology
        embedding = self.embedding
        if fragment is not None:
            if nput.is_int(fragment):
                fragment = topology.fragment_indices[fragment]
            base_ints = self._take_subchoice(fragment).get_bond_graph_internals(
                include_stretches=include_stretches,
                include_bends=include_bends,
                include_dihedrals=include_dihedrals,
                include_fragments=include_fragments,
                base_internals=base_internals,
                pruning=pruning,
                concatenate=concatenate,
                gmatrix=gmatrix
            )
            if concatenate:
                return coordops.permute_internals(base_ints, fragment)
            else:
                return [
                    coordops.permute_internals(b, fragment)
                    for b in base_ints
                ]
        else:
            st, bo, di = coordops.get_stretch_coordinate_system(
                [tuple(b[:2]) for b in topology.bonds],
                include_bends=include_bends,
                include_dihedrals=include_dihedrals
            )
            bits = []
            if include_fragments:
                if use_distance_matrix:
                    dm = nput.distance_matrix(embedding.coords)
                else:
                    dm = None
                frag_bits = coordops.get_fragment_coordinate_system(
                    topology.edge_graph,
                    masses=embedding.masses,
                    distance_matrix=dm
                )
                bits.append(frag_bits)
            if include_stretches:
                bits.append(st)
            if include_bends:
                bits.append(bo)
            if include_dihedrals:
                bits.append(di)

            if concatenate:
                internals = bits[0]
                for b in bits[1:]:
                    internals = internals + b

                if pruning:
                    if pruning is True:
                        pruning = self.default_coordinate_pruning
                    internals = self.prune_internals(internals, method=pruning, gmatrix=gmatrix)
            else:
                if pruning:
                    raise ValueError("can't prune without concatenating")
                internals = bits

            return internals

    def prune_internals(self, coords, method='b_matrix', check_rigidity=True, gmatrix=None):
        """
        Reduce a set of internal coordinates down to a non-redundant, well-conditioned subset,
        defaulting to a B-matrix-rank-based method (building the necessary
        translation/rotation-projected B-matrix generator and a sensible `max_coords` cap) if no
        custom method is supplied.

        The B-matrix generator needs the (square-root) G-matrix. Since this object doesn't hold
        a `Molecule` (and so can't compute a G-matrix itself), the caller must pass it in as
        `gmatrix` whenever the B-matrix method is used and no explicit `'b_matrix'` generator is
        already present in `method`.

        :param coords: the internal-coordinate specs to prune
        :type coords: list
        :param method: the pruning method: a method-name string, or a dict of method options
            (with a `'method'` key defaulting to `'b_matrix'`)
        :type method: str | dict
        :param check_rigidity: whether to check that the pruned coordinate set spans a rigid
            (non-redundant) representation
        :type check_rigidity: bool
        :param gmatrix: the fractional-power (`G^{1/2}`) G-matrix to build the B-matrix
            generator from; supplied by the caller (the `Molecule`, via
            `get_gmatrix(power=1/2)`) and required when the B-matrix method is used
        :type gmatrix: np.ndarray | None
        :return: the pruned coordinate specs
        :rtype: list
        :raises ValueError: if the B-matrix pruning method is used but no `gmatrix` (nor an
            explicit `'b_matrix'` generator) was provided
        """
        embedding = self.embedding
        topology = self.topology
        if isinstance(method, str):
            method = {'method':method}
        if hasattr(method, 'items'):
            meth = method.get('method')
            if meth is None:
                method = method.copy()
                method['method'] = 'b_matrix'
                meth = 'b_matrix'
            if dev.str_is(meth, 'b_matrix'):
                if 'b_matrix' not in method:
                    if gmatrix is None:
                        raise ValueError(
                            "b_matrix pruning needs the G^1/2 matrix; pass it in as `gmatrix`"
                            " (e.g. from `mol.get_gmatrix(power=1/2)`)"
                        )
                    g12 = gmatrix
                    coords_cart = embedding.coords
                    proj = nput.translation_rotation_projector(coords_cart, embedding.masses, mass_weighted=True)
                    def b_gen(pos, crds):
                        """
                        Compute the (translation/rotation-projected, mass-weighted) B-matrix for
                        a candidate set of internal coordinates at this geometry, used by the
                        default `'b_matrix'` pruning method to assess rank/conditioning.

                        :param pos: the coordinate index/indices under consideration (unused
                            directly in the body, but part of the callback signature expected by
                            the pruning routine)
                        :type pos: Any
                        :param crds: the candidate coordinate specs to build the B-matrix for
                        :type crds: list
                        :return: the projected, mass-weighted B-matrix
                        :rtype: np.ndarray
                        """
                        return proj @ g12 @ nput.internal_coordinate_tensors(coords_cart, crds, order=1)[1]
                    method = method.copy()
                    method['b_matrix'] = b_gen
                if 'max_coords' not in method:
                    method = method.copy()
                    method['max_coords'] =  min([3 * len(topology.atoms) - 6, len(coords)])
        return coordops.prune_internal_coordinates(
            coords,
            method=method,
            check_rigidity=check_rigidity,
        )

    def get_labeled_internals(self,
                              coordinate_filter=None,
                              allowed_coordinate_types=None,
                              excluded_coordinate_types=None,
                              allowed_ring_types=None,
                              excluded_ring_types=None,
                              allowed_group_types=None,
                              excluded_group_types=None,
                              include_stretches=True,
                              include_bends=True,
                              include_dihedrals=True,
                              include_fragments=True,
                              coordinate_sorting=None,
                              pruning=False,
                              gmatrix=None
                              ):
        """
        Build the internal coordinates from the bonding graph (via `get_bond_graph_internals`)
        and label each one by its atom types/ring/functional-group membership (via
        `edge_graph.get_label_types` and `coordops.get_coordinate_label`), then filter and sort
        them.

        :param coordinate_filter: an explicit filter function to apply instead of building one
            from the allow/exclude arguments
        :type coordinate_filter: callable | None
        :param allowed_coordinate_types: forwarded to `get_coordinate_filer` if
            `coordinate_filter` is not given
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: forwarded to `get_coordinate_filer`
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: forwarded to `get_coordinate_filer`
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: forwarded to `get_coordinate_filer`
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: forwarded to `get_coordinate_filer`
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: forwarded to `get_coordinate_filer`
        :type excluded_group_types: Iterable | None
        :param include_stretches: whether to include bond-stretch coordinates
        :type include_stretches: bool
        :param include_bends: whether to include bond-angle coordinates
        :type include_bends: bool
        :param include_dihedrals: whether to include dihedral-angle coordinates
        :type include_dihedrals: bool
        :param include_fragments: whether to include inter-fragment coordinates
        :type include_fragments: bool
        :param coordinate_sorting: a custom sorting function to apply to the labeled coordinates
            instead of the default `coordops.sort_internal_coordinates`; pass a falsy value to
            skip sorting
        :type coordinate_sorting: callable | bool | None
        :param pruning: whether/how to prune the coordinates, forwarded to
            `get_bond_graph_internals`
        :type pruning: bool | str | dict
        :param gmatrix: the fractional-power G-matrix passed through to
            `get_bond_graph_internals` (needed only if B-matrix pruning is requested)
        :type gmatrix: np.ndarray | None
        :return: a mapping from coordinate spec to its label, filtered and sorted
        :rtype: dict
        """
        internals = self.get_bond_graph_internals(
            include_stretches=include_stretches,
            include_bends=include_bends,
            include_dihedrals=include_dihedrals,
            include_fragments=include_fragments,
            pruning=pruning,
            gmatrix=gmatrix
        )

        labels = self.topology.edge_graph.get_label_types()
        internals = {
            (c if isinstance(c, tuple) else coordops.InternalCoordinateType.resolve(c)): coordops.get_coordinate_label(
                c,
                labels
            )
            for c in internals
        }

        if coordinate_filter is None:
            coordinate_filter = self.get_coordinate_filer(
                allowed_coordinate_types=allowed_coordinate_types,
                excluded_coordinate_types=excluded_coordinate_types,
                allowed_ring_types=allowed_ring_types,
                excluded_ring_types=excluded_ring_types,
                allowed_group_types=allowed_group_types,
                excluded_group_types=excluded_group_types,
            )

        if coordinate_filter:
            internals = coordinate_filter(internals)

        if coordinate_sorting is None:
            coordinate_sorting = coordops.sort_internal_coordinates

        if coordinate_sorting:
            internals = coordinate_sorting(internals)

        return internals

    def get_mode_labels(self,
                        internals=None,
                        modes=None,
                        use_redundants=True,
                        expansions=None,
                        return_modes=False,
                        gmatrix=None,
                        **internals_opts
                        ):
        """
        Assign human-readable labels (e.g. "C-H stretch") to a set of normal modes by projecting
        them onto labeled internal coordinates, handling both redundant and non-redundant
        internal-coordinate expansions and both Cartesian- and internal-coordinate-basis modes.

        Since this object doesn't hold a `Molecule`, the normal `modes` and (for the
        Cartesian-mode branch) the `gmatrix` must be supplied by the caller (the `Molecule`,
        via `get_normal_modes()` and `get_gmatrix()`).

        :param internals: the labeled internal coordinates to project onto; computed via
            `get_labeled_internals` if not given
        :type internals: dict | None
        :param modes: the normal modes to label; must be provided by the caller
        :type modes: Any
        :param use_redundants: whether to build a redundant-coordinate expansion (with
            relocalization) for the projection, rather than using the internal coordinates
            directly
        :type use_redundants: bool
        :param expansions: precomputed `(expansions, inverse_expansion)` internal-coordinate
            Jacobian data to reuse instead of recomputing it
        :type expansions: tuple | None
        :param return_modes: whether to also return the internal-coordinate-basis mode matrix
            alongside the labels
        :type return_modes: bool
        :param gmatrix: the G-matrix, supplied by the caller (`mol.get_gmatrix()`); required for
            the Cartesian-mode projection branch
        :type gmatrix: np.ndarray | None
        :param internals_opts: extra options forwarded to `get_labeled_internals` if `internals`
            is not given
        :type internals_opts: dict
        :return: the mode labels, or `(internal_modes, labels)` if `return_modes` is set
        :rtype: list | tuple
        :raises ValueError: if `modes` is not provided, or if the Cartesian-mode branch is
            reached without a `gmatrix`
        """
        embedding = self.embedding
        if modes is None:
            raise ValueError(
                "get_mode_labels needs the normal modes; pass them in as `modes`"
                " (e.g. from `mol.get_normal_modes()`)"
            )
        modes = modes.remove_mass_weighting()

        if internals is None:
            internals = self.get_labeled_internals(**internals_opts)

        if modes.is_cartesian:
            if expansions is not None:
                expansions, inv_expansion = expansions
            else:
                expansions = inv_expansion = None

            if use_redundants:
                redundant_tf, expansions = coordops.RedundantCoordinateGenerator(
                    internals,
                    masses=embedding.masses,
                    relocalize=True
                ).compute_redundant_expansions(embedding.coords,
                                               expansions=expansions
                                               )

                redund_labs = coordops.get_mode_labels(
                    internals,
                    redundant_tf,
                    norm_cutoff=.3
                )

                inv_expansion = nput.inverse_internal_coordinate_tensors(
                    expansions,
                    coords=embedding.coords,
                    masses=embedding.masses,
                    order=1,
                    remove_translation_rotation=True
                )

            else:
                redund_labs = internals
                if expansions is None:
                    expansions, inv_expansion = nput.internal_coordinate_tensors(
                        embedding.coords,
                        internals,
                        order=1,
                        masses=embedding.masses,
                        return_inverse=True
                    )
                    expansions = expansions[1:]

            if gmatrix is None:
                raise ValueError(
                    "labeling Cartesian modes needs the G-matrix; pass it in as `gmatrix`"
                    " (e.g. from `mol.get_gmatrix()`)"
                )
            g = expansions[0].T @ gmatrix @ expansions[0]
            g12 = nput.fractional_power(g, 1 / 2)
            internal_modes = g12 @ inv_expansion[0] @ modes.modes_by_coords
        else:
            redund_labs = internals
            internal_modes = modes

        labels = coordops.get_mode_labels(
            redund_labs,
            internal_modes,
            norm_cutoff=.8
        )
        if return_modes:
            return internal_modes, labels
        else:
            return labels
    #endregion
