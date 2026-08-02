"""
Provides a `MolecularTopology` class that encapsulates all of the bond-topology-related
functionality that used to live directly on `Molecule`.

The topology owns the atom labels, the bond list, the (lazily built) `EdgeGraph`, the
bond-guessing configuration, and the (lazily computed) fragmentation, along with all of the
graph-based algorithms (backbone/segment finding, Z-matrix construction from the bonding
graph, neighborhood/path queries, etc.).

A `MolecularTopology` is self-contained: it needs only the atom labels and bonds, not a
reference to a `Molecule`. Bond guessing is delegated to a caller-supplied `bond_guesser`
callback (the parent `Molecule` hooks this up to go through the `RDMolecule` pass), and the
few routines that need geometry (currently just `get_bond_zmatrix`'s inter-fragment
distance matrix) accept `coords` explicitly.
"""

from __future__ import annotations

import numpy as np

import McUtils.Numputils as nput
import McUtils.Iterators as itut
import McUtils.Coordinerds as coordops
from McUtils.Graphs import EdgeGraph

__all__ = [
    "MolecularTopology"
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

    def neighborhood(self, loc, size=1):
        """
        Find the atoms within a given graph-distance of a location in the bonding graph.

        :param loc: the atom index to center the neighborhood on
        :type loc: int
        :param size: the neighborhood radius (in bond-graph steps)
        :type size: int
        :return: the neighboring atom indices
        :rtype: tuple[int]
        """
        return tuple(l for l in self.edge_graph.neighbor_iterator(loc, num=size))
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
