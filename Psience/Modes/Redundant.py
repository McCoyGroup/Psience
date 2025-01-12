import itertools

import numpy as np
import collections
import enum
import scipy.sparse as sparse
import McUtils.Numputils as nput
from McUtils.ExternalPrograms import RDMolecule

__all__ = [
    "EdgeGraph",
    "PrimitiveCoordinatePicker",
    "RedundantCoordinateGenerator"
]

class EdgeGraph:
    __slots__ = ["labels", "edges", "graph", "map"]
    def __init__(self, labels, edges, graph=None, edge_map=None):
        self.labels = labels
        self.edges = np.asanyarray(edges)
        if graph is None:
            graph = self.adj_mat(len(labels), self.edges)
        self.graph = graph
        if edge_map is None:
            edge_map = self.build_edge_map(self.edges)
        self.map = edge_map

    @classmethod
    def adj_mat(cls, num_nodes, edges):
        adj = np.zeros((num_nodes, num_nodes), dtype=bool)
        rows,cols = edges.T
        adj[rows, cols] = 1
        adj[cols, rows] = 1

        return sparse.csr_matrix(adj)

    @classmethod
    def build_edge_map(cls, edge_list):
        map = {}
        for e1,e2 in edge_list:
            if e1 not in map: map[e1] = set()
            map[e1].add(e2)
            if e2 not in map: map[e2] = set()
            map[e2].add(e1)
        return map

    @classmethod
    def _remap(cls, labels, pos, rows, cols):
        if len(rows) == 0:
            edge_list = np.array([], dtype=int).reshape(-1, 2)
        else:
            new_mapping = np.zeros(len(labels), dtype=int)
            new_mapping[pos,] = np.arange(len(pos))
            new_row = new_mapping[rows,]
            new_col = new_mapping[cols,]
            edge_list = np.array([new_row, new_col]).T

        return [labels[p] for p in pos], edge_list
    @classmethod
    def _take(cls, pos, labels, adj_mat:sparse.compressed):
        (rows, cols), _ = sparse.find(adj_mat)
        utri = cols >= rows
        rows = rows[utri]
        cols = cols[utri]
        row_cont, _ = nput.contained(rows, pos)
        col_cont, _ = nput.contained(cols, pos)
        cont = np.logical_and(row_cont, col_cont)

        labels, edge_list = cls._remap(labels, pos, rows[cont], cols[cont])
        return cls(labels, edge_list)

    def take(self, pos):
        return self._take(pos, self.labels, self.graph)

    def split(self, backbone_pos):
        new_adj = self.graph.copy()
        for n in backbone_pos:
            for i,j in self.map[n]:
                new_adj[i,j] = 0
        ncomp, labels = sparse.csgraph.connected_components(new_adj, directed=False, return_labels=True)
        groups, _ = nput.group_by(np.arange(len(labels)), labels)
        return [
            self._take(pos, self.labels, new_adj)
            for _, pos in groups
        ]

    @classmethod
    def _bfs(cls, test, root, edge_map, visited:set):
        if root in visited:
            return
        queue = collections.deque([root])
        while queue:
            head = queue.pop()
            visited.add(head)
            test(head)
            queue.extend(h for h in edge_map.get(head, []) if h not in visited)
    @classmethod
    def _subgraph_match(cls,
                        root1, labels1, edge_map1,
                        root2, labels2, edge_map2,
                        visited=None
                        ):
        if labels1[root1] != labels2[root2]:
            return False
        if len(edge_map1.get(root1, [])) != len(edge_map2.get(root2, [])):
            return False

        if visited is None:
            visited = (set(), set())
        visit1, visit2 = visited

        # queue1 = collections.deque([root1])
        tests1 = set(edge_map1.get(root1, [])) - visit1
        tests2 = set(edge_map2.get(root2, [])) - visit2
        for r1 in tests1:
            visit1.add(r1)
            for r2 in tests2:
                visit2.add(r2)
                # check if all the subgraphs match via DFS
                # we could do this without recursion, but it'd be
                # more annoying
                if cls._subgraph_match(
                    r1, labels1, edge_map1,
                    r2, labels2, edge_map2,
                    visited=(visit1, visit2)
                ):
                    tests2.remove(r2)
                    break
                else:
                    visit2.remove(r2)

        return len(tests2) == 0

    @classmethod
    def graph_match(cls, graph1:'EdgeGraph', graph2:'EdgeGraph'):
        # we do some quick prunes
        atoms1 = graph1.labels
        atoms2 = graph2.labels
        if (
                len(atoms1) != len(atoms2)
                or atoms1[0] != atoms2[0]
                or len(graph1.edges) != len(graph2.edges)
                or list(sorted(atoms1)) != list(sorted(atoms2))
                or list(sorted(len(v) for v in graph1.map.values())) != list(sorted(len(v) for v in graph2.map.values()))
        ):
            return False

        return cls._subgraph_match(
            0, graph1.labels, graph1.map,
            0, graph2.labels, graph2.map
        )

    def __eq__(self, other):
        return self.graph_match(self, other)

    @classmethod
    def build_neighborhood_graph(cls, node, labels, edge_map, ignored=None, num=1):
        edges = []
        visited = set()
        if ignored is None: ignored = []
        ignored = set(ignored)
        queue = [node]
        for i in range(num):
            new_queue = []
            for node in queue:
                visited.add(node)
                new_nodes = set(edge_map[node]) - visited - ignored
                edges.extend((node, e) for e in new_nodes)
                new_queue.extend(new_nodes)
            queue = new_queue

        edges = np.array(edges, dtype=int)
        if len(edges) == 0:
            edges = np.reshape(edges, (-1, 2))
        labels, edges = cls._remap(labels, list(visited), edges[:, 0], edges[:, 1])
        return cls(labels, edges)

    def neighbor_graph(self, root, ignored=None, num=1):
        return self.build_neighborhood_graph(root, self.labels, self.map, ignored=ignored, num=num)

    def get_rings(self):
        return RDMolecule.from_coords(
            ["C"]*len(self.labels),
            coords=np.zeros((len(self.labels), 3)),
            bonds=[[int(i), int(j), 1] for i,j in self.edges]
        ).rings

class PrimitiveCoordinatePicker:

    light_atom_types = {"H", "D"}
    def __init__(self, atoms, bonds, base_coords=None, rings=None, light_atoms=None, backbone=None, neighbor_count=3):
        self.graph = EdgeGraph(atoms, bonds)
        if rings is None:
            rings = self.graph.get_rings()
        self.rings = rings
        self.ring_atoms = {k:n for n,rats in enumerate(rings) for k in rats}
        self.backbone = backbone
        self.light_atoms = (
            [
                i for i, l in enumerate(self.graph.labels)
                if l in self.light_atom_types
            ]
                if light_atoms is None else
            light_atoms
        )
        self.base_coords = list(base_coords) if base_coords is not None else []
        self.neighbors = neighbor_count
        self._coords = None
    @property
    def coords(self):
        if self._coords is None:
            self._coords = tuple(self.generate_coords())
        return self._coords
    def generate_coords(self):
        coords = []
        for ring in self.rings:
            coords.extend(self.ring_coordinates(ring))
        for i,r1 in enumerate(self.rings):
            for r2 in self.rings[i+1:]:
                coords.extend(self.fused_ring_coordinates(r1, r2))
        for a in range(len(self.graph.labels)):
            if a not in self.ring_atoms:
                symm_c = self.symmetry_coords(a, neighborhood=self.neighbors, backbone=self.backbone)
                coords.extend(symm_c)
        coords = self.base_coords + coords
        return self.prune_excess_coords(coords)

    @classmethod
    def canonicalize_coord(cls, coord):
        dupes = len(np.unique(coord)) < len(coord)
        if dupes: return None
        if len(coord) == 2:
            i,j = coord
            if i > j:
                coord = (j, i)
        elif len(coord) == 3:
            i,j,k = coord
            if i > k:
                coord = (k, j, i)
        elif len(coord) == 4:
            i, j, k, l = coord
            if i > l:
                coord = (l, k, j, i)
        elif coord[0] > coord[-1]:
            coord = tuple(reversed(coord))
        return coord

    @classmethod
    def prune_excess_coords(cls, coord_set, canonicalized=False):
        if not canonicalized:
            coord_set = [cls.canonicalize_coord(c) for c in coord_set]
            coord_set = [c for c in coord_set if c is not None]
        dupe_set = set()
        coords = []
        coord_counts = {}
        for coord in coord_set:
            if coord in dupe_set: continue
            dupe_set.add(coord)
            if len(coord) == 4:
                i,j,k,l = coord
                if (i, k, j, l) in dupe_set: continue
                # only one choice of ordering for shared bond, even though implies different exterior bonds
            elif len(coord) == 3:
                i,j,k = coord
                if ( # can't have all three angles, dealing with the ordering manually
                    (k < j and (i, k, j) in dupe_set and (k, i, j) in dupe_set)
                    or (k > j and (j, i, k) in dupe_set) and ((i, k, j) in dupe_set or (j, k, i) in dupe_set)
                ): continue
            coords.append(coord)
            # key = tuple(sorted(coord))
            # coord_counts[key] = coord_counts.get(key, 0) + 1
            # if len(coord_counts)

        return coords

    @classmethod
    def ring_coordinates(cls, ring_atoms):
        # ordered as pair-wise bonded

        bonds = list(zip(
            ring_atoms, ring_atoms[1:] + ring_atoms[:1]
        ))
        angles = list(zip(
            ring_atoms, ring_atoms[1:] + ring_atoms[:1], ring_atoms[2:] + ring_atoms[:2]
        ))
        dihedrals = list(zip(
            ring_atoms, ring_atoms[1:] + ring_atoms[:1],
                        ring_atoms[2:] + ring_atoms[:2],
                        ring_atoms[3:] + ring_atoms[:3],
        ))

        return bonds + angles + dihedrals

    fused_ring_dispatch_table = {

    }
    @classmethod
    def _fused_dispatch(cls):
        return dict({
            0:cls.unfused_ring_coordinates,
            1:cls.pivot_fused_ring_coordinates,
            2:cls.simple_fused_ring_coordinates
        }, **cls.fused_ring_dispatch_table)
    @classmethod
    def unfused_ring_coordinates(cls, ring_atoms1, ring_atoms2, shared_atoms, shared_indices1, shared_indices2):
        return []
    @classmethod
    def pivot_fused_ring_coordinates(cls, ring_atoms1, ring_atoms2, shared_atoms, shared_indices1, shared_indices2):
        p = shared_atoms[0]
        i = shared_indices1[0]
        n = len(ring_atoms1)
        j = shared_indices2[0]
        m = len(ring_atoms2)

        # add in all relative angles
        ip1 = (i+1) % n
        jp1 = (j+1) % m
        return [
            (ring_atoms1[i-1], p, ring_atoms2[j-1]),
            (ring_atoms1[ip1], p, ring_atoms2[j-1]),
            (ring_atoms1[i-1], p, ring_atoms2[jp1]),
            (ring_atoms1[ip1], p, ring_atoms2[jp1])
        ]

    @classmethod
    def simple_fused_ring_coordinates(cls, ring_atoms1, ring_atoms2, shared_atoms, shared_indices1, shared_indices2):
        j, k = shared_atoms
        j1, k1 = shared_indices1
        j2, k2 = shared_indices2
        # we want relative orientation indices for both rings
        inds = []
        if j1 > k1:
            j1, k1 = k1, j1
        if j2 > k2:
            j2, k2 = k2, j2
        i1 = ring_atoms1[j1 - 1]; l1 = ring_atoms1[(k1 + 1) % len(ring_atoms1)]
        i2 = ring_atoms2[j2 - 1]; l2 = ring_atoms2[(k2 + 1) % len(ring_atoms2)]

        return [
            (i1, j, k, l2),
            (i2, j, k, l1)
        ]

    @classmethod
    def fused_ring_coordinates(cls, ring_atoms1, ring_atoms2):
        shared_atoms, _, _, r1_indices, r2_indices = nput.intersection(ring_atoms1, ring_atoms2, return_indices=True)
        if len(shared_atoms) == 0:
            return []
        else:
            n = len(shared_atoms)
            coord_func = cls._fused_dispatch().get(n)
            if coord_func is None:
                raise ValueError(f"can't deal with fused rings with {n} shared atoms")
            return coord_func(ring_atoms1, ring_atoms2, shared_atoms, r1_indices, r2_indices)


    def get_neighborhood_symmetries(self, atoms, ignored=None, neighborhood=3):
        graphs = [
            self.graph.neighbor_graph(a, ignored=ignored, num=neighborhood)
            for a in atoms
        ]
        rows, cols = np.triu_indices(len(graphs), k=1)
        return [graphs[r] == graphs[c] for r,c in zip(rows, cols)]

    def chain_coords(self, R, y):
        coords = []
        if len(R) > 0:
            coords.append((y, R[-1]))
        if len(R) > 1:
            coords.append((y, R[-1], R[-2]))
        if len(R) > 2:
            coords.append((y, R[-1], R[-2], R[-3]))
        return coords

    def RYX2_coords(self, R, y, X):
        coords = []
        coords.extend((y, x) for x in X)
        coords.extend(
            (X[i], y, X[j])
            for i, j in itertools.combinations(range(3), 2)
        )
        if len(R) > 0:
            coords.append((R[-1], y))
            # add in RYX angles
            coords.extend(
                (R[-1], y, x)
                for x in X
            )
        if len(R) > 1:
            # add in RRY angle
            coords.append(
                (R[-2], R[-1], y)
            )
        if len(R) > 2:
            # add in RRRY dihedral
            coords.append(
                (R[-3], R[-2], R[-1], y)
            )

        return coords

    def RYX3_coords(self, R, y, X):
        coords = []
        coords.extend((y, x) for x in X)
        coords.extend(
            (X[i], y, X[j])
            for i, j in itertools.combinations(range(3), 2)
        )
        if len(R) > 0:
            coords.append((R[-1], y))
            # add in RYX angles
            coords.extend(
                (R[-1], y, x) for x in X
            )
        if len(R) > 1:
            # add in RRY angle
            coords.append(
                (R[-2], R[-1], y)
            )
        if len(R) > 2:
            # add in RRRY dihedral
            coords.append(
                (R[-3], R[-2], R[-1], y)
            )

        return coords

    def get_precedent_chains(self, atom, num_precs=2, ring_atoms=None, light_atoms=None, ignored=None, backbone=None):
        chains = []
        visited = set([] if ignored is None else ignored)
        ring_atoms = set(self.ring_atoms if ring_atoms is None else ring_atoms)
        light_atoms = set(self.light_atoms if light_atoms is None else light_atoms)
        if backbone is not None:
            backbone = set(backbone)

        visited = visited #| ring_atoms # | light_atoms

        # do a dfs exploration up to the given depth over non-ring, heavy atoms
        queue = collections.deque([[[], 0, atom]])
        while queue:
            chain, depth, root = queue.pop()
            neighbors = self.graph.map[root]
            visited.add(root)
            branches = neighbors - visited
            if len(branches) == 0:
                chains.append(chain)
            elif depth == num_precs - 1:
                chains.extend(chain + [n] for n in branches)
            else:
                queue.extend([chain + [n], depth+1, n] for n in branches)

            #
            # queue.extend(rem)
            # while len(rem) == 0 and queue:
            #     # walk up dfs tree until we find a branch with nodes that work
            #     root = queue.pop()
            #     if root not in visited:
            #         rem = {root}
            # else:
            #     if len(rem) == 0: rem = neighbors - visited - light_atoms
            #     # if len(rem) == 0: rem = neighbors - visited
            #     if len(rem) == 0: break
            #
            #     if backbone is not None:
            #         bb_chain = rem & backbone
            #     else:
            #         bb_chain = []
            #     if len(bb_chain) > 0:
            #         atom = min(bb_chain)
            #     else:
            #         atom = min(rem)
            #     chain.append(atom)

        # chain = list(reversed(chain))

        return [list(reversed(c)) for c in chains]

    symmetry_type_dispatch = {}
    def _symmetry_dispatch(self):
        return dict({
            (3,):self._3_coords,
            (2,1):self._2_1_coords,
            (3,1):self._3_1_coords
        })
    def _2_1_coords(self, atom, neighbors, X, R, backbone=None):
        coords = []
        R = R[0]
        chains = self.get_precedent_chains(R, 1, ignored=[atom] + neighbors, backbone=backbone)
        if len(chains) == 0: chains = [[]]
        for c in chains:
            coords.extend(self.RYX2_coords(c + [R], atom, X))
        return coords
    def _3_coords(self, atom, neighbors, X, backbone=None):
        return self.RYX3_coords([], atom, X)
    def _3_1_coords(self, atom, neighbors, X, R, backbone=None):
        coords = []
        R = R[0]
        chains = self.get_precedent_chains(R, 1, ignored=[atom] + neighbors, backbone=backbone)
        if len(chains) == 0: chains = [[]]
        for c in chains:
            coords.extend(self.RYX3_coords(c + [R], atom, X))
        return coords

    # def _2_2_coords(self, atom, neighbors, X, R, backbone=None):
    #     coords = []
    #     R = R[0]
    #     # chains = self.get_precedent_chains(R, 1, ignored=[atom] + neighbors, backbone=backbone)
    #     if len(chains) == 0: chains = [[]]
    #     for c in chains:
    #         coords.extend(self.R2YX2_coords(c + [R], atom, X))
    #     return coords


    @classmethod
    def get_symmetry_groups(cls, neighbors, matches):
        groups = {}
        n = len(neighbors)
        k = 0
        for i in range(n):
            if i not in groups:
                groups[i] = {i}
            for j in range(i+1, n):
                eq = matches[k]
                k += 1
                if eq:
                    groups[j] = groups[i]
                    groups[i].add(j)
        groups = {id(g):g for g in groups.values()}
        groups = list(reversed(sorted(groups.values(), key=lambda g:len(g))))
        return [
            [neighbors[i] for i in g]
            for g in groups
        ]

    def symmetry_coords(self, atom, neighborhood=3, backbone=None):
        # neighbors = list(self.graph.map[atom])
        coords = []
        # dispatch = self._symmetry_dispatch()
        # neighbor_counts = {sum(k) for k in dispatch.keys()}
        # if len(neighbors) in neighbor_counts:
        #     symms = self.get_neighborhood_symmetries(neighbors, ignored=[atom], neighborhood=neighborhood)
        #     groups = self.get_symmetry_groups(neighbors, symms)
        #     key = tuple(len(g) for g in groups)
        #     dfunc = dispatch.get(key)
        #     if dfunc is not None:
        #         coords = dfunc(atom, neighbors, *groups, backbone=backbone)
        chains = self.get_precedent_chains(atom, 3, backbone=backbone)
        if len(chains) == 0: chains = [[]]
        for R in chains:
            coords.extend(self.chain_coords(R, atom))

        # if coords is None:
        #     R = self.get_precedent_chain(atom, 3, backbone=backbone)
        #     coords = self.chain_coords(R, atom)

        return coords


class RedundantCoordinateGenerator:

    def __init__(self, coordinate_specs, angle_ordering='ijk', untransformed_coordinates=None, masses=None,
                 relocalize=False,
                 **opts):
        self.specs = coordinate_specs
        self.untransformed_coordinates = untransformed_coordinates
        self.relocalize = relocalize
        self.masses = masses
        self.opts = dict(
            opts,
            angle_ordering=angle_ordering
        )

    @classmethod
    def _pad_redund_tf(cls, redund_tf, n):
        ndim = redund_tf.ndim - 2
        eye = nput.identity_tensors(redund_tf.shape[:-2], n)
        return np.concatenate(
            [
                np.pad(eye, ([[0, 0]] * ndim) + [[0, redund_tf.shape[-2]], [0, 0]]),
                np.pad(redund_tf, ([[0, 0]] * ndim) + [[n, 0], [0, 0]])
            ],
            axis=-1
        )

    @classmethod
    def _relocalize_tf(cls, redund_tf):
        n = redund_tf.shape[-1]
        # target = np.pad(np.eye(n), [[0, 173 - 108], [0, 0]])
        ndim = redund_tf.ndim - 2
        eye = nput.identity_tensors(redund_tf.shape[:-2], n)
        target = np.pad(eye, ([[0, 0]] * ndim) + [[0, redund_tf.shape[-2] - n], [0, 0]])
        loc = np.linalg.lstsq(redund_tf, target, rcond=None)
        U, s, V = np.linalg.svd(loc[0])
        R = U @ V
        return redund_tf @ R


    @classmethod
    def base_redundant_transformation(cls, expansion,
                                      untransformed_coordinates=None,
                                      masses=None,
                                      relocalize=False
                                      ):
        conv = np.asanyarray(expansion[0])
        if untransformed_coordinates is not None:
            transformed_coords = np.setdiff1d(np.arange(conv.shape[-1]), untransformed_coordinates)
            ut_conv = conv[..., untransformed_coordinates]
            conv = conv[..., transformed_coords]

            # project out contributions along untransformed coordinates to ensure
            # dimension of space remains unchanged
            ut_conv = nput.vec_normalize(ut_conv)
            proj = ut_conv @ np.moveaxis(ut_conv, -1, -2)
            proj = nput.identity_tensors(proj.shape[:-2], proj.shape[-1]) - proj

            conv = proj @ conv
        else:
            transformed_coords = None

        if masses is None:
            G_internal = np.moveaxis(conv, -1, -2) @ conv
        else:
            masses = np.asanyarray(masses)
            if len(masses) == conv.shape[-2] // 3:
                masses = np.repeat(masses, 3)
            masses = np.diag(1 / masses)
            if conv.ndim > 2:
                masses = np.broadcast_to(
                    np.expand_dims(masses, list(range(conv.ndim - 2))),
                    conv.shape[:-2] + masses.shape
                )
            G_internal = np.moveaxis(conv, -1, -2)  @ masses @ conv


        redund_vals, redund_tf = np.linalg.eigh(G_internal)
        redund_pos = np.where(np.abs(redund_vals) > 1e-10)
        if redund_vals.ndim > 1:
            redund_tf = nput.take_where_groups(redund_tf, redund_pos)
        else:
            redund_tf = redund_tf[:, redund_pos[0]]
        if isinstance(redund_tf, np.ndarray):
            redund_tf = np.flip(redund_tf, axis=-1)
        else:
            redund_tf = [
                np.flip(tf, axis=-1)
                for tf in redund_tf
            ]
        if transformed_coords is not None:
            n = len(untransformed_coordinates)
            if isinstance(redund_tf, np.ndarray):
                redund_tf = cls._pad_redund_tf(redund_tf, n)
            else:
                redund_tf = [
                    cls._pad_redund_tf(tf, n)
                    for tf in redund_tf
                ]

        if relocalize:
            # provide some facility for rearranging coords?
            if isinstance(redund_tf, np.ndarray):
                redund_tf = cls._relocalize_tf(redund_tf)
            else:
                redund_tf = [
                    cls._relocalize_tf(tf)
                    for tf in redund_tf
                ]



        return redund_tf


    @classmethod
    def get_redundant_transformation(cls, base_expansions, untransformed_coordinates=None, masses=None,
                                     relocalize=False):
        if isinstance(base_expansions, np.ndarray) and base_expansions.ndim == 2:
            base_expansions = [base_expansions]
        redund_tf = cls.base_redundant_transformation(base_expansions,
                                                      untransformed_coordinates=untransformed_coordinates,
                                                      masses=masses,
                                                      relocalize=relocalize
                                                      )

        # dQ/dR, which we can transform with dR/dX to get dQ/dX
        if isinstance(redund_tf, np.ndarray):
            redund_expansions = nput.tensor_reexpand(base_expansions, [redund_tf])
        else:
            redund_expansions = [
                nput.tensor_reexpand(base_expansions, [tf])
                for tf in redund_tf
            ]

        return redund_tf, redund_expansions

    def compute_redundant_expansions(self, coords, order=None, untransformed_coordinates=None,
                                     relocalize=None):
        coords = np.asanyarray(coords)
        if order is None:
            opts = dict(dict(order=1), **self.opts)
        else:
            opts = dict(self.opts, order=order)
        if untransformed_coordinates is None:
            untransformed_coordinates = self.untransformed_coordinates
        if relocalize is None:
            relocalize = self.relocalize
        base_expansions = nput.internal_coordinate_tensors(coords, self.specs, **opts)
        return self.get_redundant_transformation(base_expansions,
                                                 untransformed_coordinates=untransformed_coordinates,
                                                 masses=self.masses,
                                                 relocalize=relocalize
                                                 )

    @classmethod
    def _prune_coords_svd(cls, b_mat, svd_cutoff=5e-2, sort=True, fixed_vecs=None):
        # turns out, equivalent to finding maximimum loc in eigenvectors of G
        if fixed_vecs is not None:
            transformed_coords = np.delete(np.arange(b_mat.shape[-1]), fixed_vecs)
            ut_conv = b_mat[..., fixed_vecs]
            conv = b_mat[..., transformed_coords]

            # project out contributions along untransformed coordinates to ensure
            # dimension of space remains unchanged
            ut_conv = nput.vec_normalize(ut_conv)
            proj = ut_conv @ np.moveaxis(ut_conv, -1, -2)
            proj = nput.identity_tensors(proj.shape[:-2], proj.shape[-1]) - proj

            b_mat = proj @ conv

        _, s, Q = np.linalg.svd(b_mat)
        pos = np.where(s > 1e-8)
        loc_val = np.max(Q[pos]**2, axis=0)
        coords = np.where(loc_val > svd_cutoff)[0]
        if sort:
            coords = coords[np.argsort(-loc_val[coords,],)]
        if fixed_vecs is not None:
            coords = transformed_coords[coords,]
            coords = np.concatenate([fixed_vecs, coords])
        return coords

    @classmethod
    def _prune_coords_loc(cls, b_mat, loc_cutoff=.33, sort=True, fixed_vecs=None):
        # if fixed_vecs is not None:
        #     transformed_coords = np.delete(np.arange(b_mat.shape[-1]), fixed_vecs)
        #     ut_conv = b_mat[..., fixed_vecs]
        #     b_mat = b_mat[..., transformed_coords]
        #
        #     # # project out contributions along untransformed coordinates to ensure
        #     # # dimension of space remains unchanged
        #     # ut_conv = nput.vec_normalize(ut_conv)
        #     # proj = ut_conv @ np.moveaxis(ut_conv, -1, -2)
        #     # proj = nput.identity_tensors(proj.shape[:-2], proj.shape[-1]) - proj
        #     #
        #     # b_mat = proj @ conv

        s, Q = np.linalg.eigh(b_mat.T @ b_mat)
        pos = np.where(s > 1e-8)[0]
        Q = Q[:, pos]
        loc_mat = np.linalg.lstsq(Q, np.eye(Q.shape[-2]), rcond=None)
        loc_val = np.diag(Q @ loc_mat[0])
        coords = np.where(loc_val > loc_cutoff)[0]
        if fixed_vecs is not None:
            coords = np.setdiff1d(coords, fixed_vecs)
        if sort:
            coords = coords[np.argsort(-loc_val[coords,], )]
        if fixed_vecs is not None:
            coords = np.concatenate([fixed_vecs, coords])

            # have to run this back again in case these added coords break the localization
            # a third run back isn't work it since that will just repeat what's found here
            b_mat = b_mat[:, coords]
            s, Q = np.linalg.eigh(b_mat.T @ b_mat)
            pos = np.where(s > 1e-8)[0]
            Q = Q[:, pos]
            loc_mat = np.linalg.lstsq(Q, np.eye(Q.shape[-2]), rcond=None)
            loc_val = np.diag(Q @ loc_mat[0])
            subcoords = np.setdiff1d(np.where(loc_val > loc_cutoff)[0], np.arange(len(fixed_vecs)))
            if sort:
                subcoords = subcoords[np.argsort(-loc_val[subcoords,], )]
            coords = np.concatenate([coords[:len(fixed_vecs)], coords[subcoords]])

        return coords

    @classmethod
    def _prune_coords_gs(cls, b_mat, fixed_vecs=None, core_scaling=1e3, max_condition_number=1e8):
        if fixed_vecs is None:
            _, s, Q = np.linalg.svd(b_mat)
            pos = np.where(s > 0)
            fixed_vecs = np.where(np.linalg.norm(Q[pos], axis=0) > .9)[0]
        fixed_vecs = np.array(fixed_vecs)
        core_mask = np.full(b_mat.shape[1], False)
        core_mask[fixed_vecs] = True
        all_mask = core_mask.copy()
        rem_pos = np.setdiff1d(np.arange(b_mat.shape[1]), fixed_vecs)
        base_evals = np.linalg.eigvalsh(b_mat[:, core_mask].T @ b_mat[:, core_mask])
        base_cond = base_evals[-1] / base_evals[0] * core_scaling
        for r in rem_pos:
            core_mask[r] = True
            evals = np.linalg.eigvalsh(b_mat[:, core_mask].T @ b_mat[:, core_mask])
            new_cond = abs(evals[-1] / evals[0])
            if new_cond < max_condition_number:
                all_mask[r] = True
                if new_cond > base_cond:
                    core_mask[r] = False
            else:
                core_mask[r] = False

        return np.where(all_mask)[0]

    @classmethod
    def prune_coordinate_specs(cls, expansion,
                               masses=None,
                               untransformed_coordinates=None,
                               pruning_mode='loc',
                               **opts
                               ):
        conv = np.asanyarray(expansion[0])
        masses = np.asanyarray(masses)
        if len(masses) == conv.shape[-2] // 3:
            masses = np.repeat(masses, 3)
        masses = np.diag(1 / np.sqrt(masses))
        b_mat = masses @ conv

        if pruning_mode == 'svd':
            return cls._prune_coords_svd(b_mat, fixed_vecs=untransformed_coordinates, **opts)
        elif pruning_mode == 'loc':
            return cls._prune_coords_loc(b_mat, fixed_vecs=untransformed_coordinates, **opts)
        elif pruning_mode == 'gs':
            return cls._prune_coords_gs(b_mat, fixed_vecs=untransformed_coordinates, **opts)
        else:
            raise ValueError(f"don't understand pruning mode {pruning_mode}")
