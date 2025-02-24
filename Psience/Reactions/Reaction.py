

__all__ = [
    "Reaction"
]

import numpy as np
import McUtils.Numputils as nput
import McUtils.Devutils as dev
from ..Molecools import Molecule

from .ProfileGenerator import ProfileGenerator


class Reaction:
    def __init__(self,
                 reactants:list[Molecule],
                 products:list[Molecule],
                 align=True,
                 atom_mapping=None,
                 use_product_structure=True,
                 fragment_expansion_method=dev.default,
                 optimize=True,
                 energy_evaluator='rdkit',
                 profile_generator=None,
                 **fragment_initialization_opts
                 ):
        self.reactants = reactants
        self.products = products

        self.product_complex = self.make_reaction_complex(products,
                                                          fragment_expansion_method=fragment_expansion_method,
                                                          **fragment_initialization_opts
                                                          )
        if use_product_structure:
            self.reactant_complex = self.make_reaction_complex(reactants,
                                                               fragment_expansion_method=None)
        else:
            self.reactant_complex = self.make_reaction_complex(reactants,
                                                               fragment_expansion_method=fragment_expansion_method,
                                                               **fragment_initialization_opts
                                                               )

        if align:
            self.reactant_complex = self.product_complex.align_molecule(
                self.reactant_complex,
                align_structures=not use_product_structure,
                permute_atoms=False if use_product_structure else 'all'
            )

        if use_product_structure:
            if fragment_expansion_method is not None:
                method, opts = self.get_expansion_dispatch().resolve(fragment_expansion_method)
                target_coords = method(self.reactant_complex.fragment_indices,
                                       self.product_complex.coords,
                                       **dict(fragment_initialization_opts, **opts)
                                       )
            else:
                target_coords = self.product_complex.coords

            self.reactant_complex = self.reactant_complex.modify(coords=target_coords)

        if optimize:
            self.product_complex = self.product_complex.optimize(evaluator=energy_evaluator)
            self.reactant_complex = self.reactant_complex.optimize(evaluator=energy_evaluator)
            if align:
                self.reactant_complex = self.product_complex.align_molecule(
                    self.reactant_complex,
                    reindex_bonds=False,
                    align_structures=True,
                    permute_atoms=False if use_product_structure else 'all'
                )

        self.profile_generator = profile_generator

    @classmethod
    def expand_fragments_centroid(cls, inds, coords,
                                  min_distance=1,
                                  add_radius=True,
                                  expansion_factor=.1,
                                  max_iterations=5):
        coords = coords.copy()
        centroids = np.array([
            np.average(coords[i, :], axis=0)
            for i in inds
        ])
        rows, cols = np.triu_indices(len(centroids), k=1)
        initial_centroid_distances = np.linalg.norm(centroids[rows,] - centroids[cols,], axis=1)
        initial_centroids = centroids.copy()

        min_distance = [min_distance] * len(centroids)
        if add_radius:
            for i,x in enumerate(inds):
                rx, cx = np.triu_indices(len(x), k=1)
                rx = [x[j] for j in rx]
                cx = [x[j] for j in cx]
                radius = np.max(np.linalg.norm(coords[rx, :] - coords[cx, :], axis=1))
                min_distance[i] += radius

        for it_no in range(max_iterations):
            for k,(i,j) in enumerate(zip(rows, cols)):
                v = centroids[i] - centroids[j]
                v_norm = np.linalg.norm(v)
                if v_norm < 1e-4:
                    v = np.array([min_distance[k]/(1 + expansion_factor/2), 0, 0])
                di = v * (1 + expansion_factor/2)
                dj = -v * (1 + expansion_factor/2)
                centroids[i] += di
                centroids[j] += dj

            centroid_distances = np.linalg.norm(centroids[rows,] - centroids[cols,], axis=1)
            if all(
                    cd > md and cd >= (1+expansion_factor) * id
                    for md,cd,id in zip(min_distance, centroid_distances, initial_centroid_distances)
            ):
                break

        displacements = centroids - initial_centroids
        for idx, d in zip(inds, displacements):
            coords[idx,] += d[np.newaxis]

        return coords

    @classmethod
    def expand_fragments_nearest(cls, inds, coords, min_distance=5, max_iterations=15):
        coords = coords.copy()
        rows, cols = np.triu_indices(len(inds), k=1)
        rows = [inds[i] for i in rows]
        cols = [inds[j] for j in cols]

        for it_no in range(max_iterations):
            for (i, j) in zip(rows, cols):
                vecs = coords[i, np.newaxis, :] - coords[np.newaxis, j, :]
                dist = np.linalg.norm(vecs, axis=2)
                x,y = np.unravel_index(np.argmin(dist), (len(i), len(j)))
                md = dist[x, y]
                print(md)
                if md < min_distance:
                    v = vecs[x, y]
                    disp = v * (min_distance / md / 2)
                    coords[i,] += disp[np.newaxis]
                    coords[j,] -= disp[np.newaxis]

            nearest_distances = np.array([
                np.min(
                    np.linalg.norm(coords[i, np.newaxis, :] - coords[np.newaxis, j, :], axis=2)
                )
                for i,j in zip(rows, cols)
            ])

            if np.all(nearest_distances > min_distance):
                break

        return coords

    @classmethod
    def expand_fragments_attractive(cls,
                                    inds, coords,
                                    min_distance=5,
                                    max_displacement=.2,
                                    max_iterations=10,
                                    tol=1e-3
                                    ):
        coords = coords.copy()
        rows, cols = np.triu_indices(len(inds), k=1)
        rows = [inds[i] for i in rows]
        cols = [inds[j] for j in cols]

        for it_no in range(max_iterations):
            disp_norm = 0
            for (i, j) in zip(rows, cols):
                vecs = coords[i, np.newaxis, :] - coords[np.newaxis, j, :]
                dist = np.linalg.norm(vecs, axis=2)
                vecs = nput.vec_normalize(vecs, norms=dist)
                lj_param = min_distance/dist
                lj_grads = -(lj_param**13 - lj_param**7)[:, :, np.newaxis] * vecs
                disps = np.sign(lj_grads) * np.clip(np.abs(lj_grads), 0, max_displacement)

                d = np.sum(np.sum(disps, axis=0), axis=0)

                coords[i,] += d[np.newaxis]
                coords[j,] -= d[np.newaxis]

                disp_norm += np.max(np.abs(disps))

            if disp_norm < tol:
                break

        return coords

    @classmethod
    def get_expansion_methods(cls):
        return {
            'nearest': cls.expand_fragments_nearest,
            'attractive': cls.expand_fragments_attractive,
            'centroid': cls.expand_fragments_centroid
        }

    _expansion_dispatch = dev.uninitialized
    @classmethod
    def get_expansion_dispatch(cls):
        cls._expansion_dispatch = dev.handle_uninitialized(
            cls._expansion_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_expansion_methods,),
            kwargs=dict(default_method='attractive')
        )
        return cls._expansion_dispatch
    @property
    def expansion_methods(self):
        return self.get_expansion_dispatch()

    @classmethod
    def initialize_reactant_structures(cls,
                                       reactant_complex: Molecule,
                                       product_complex,
                                       fragment_expansion_method=None,
                                       **expansion_method_options
                                       ):
        if fragment_expansion_method is not None:
            fragment_indices = reactant_complex.fragment_indices
            method, opts = cls.get_expansion_dispatch().resolve(fragment_expansion_method)
            target_coords = method(fragment_indices, product_complex.coords, **dict(opts, **expansion_method_options))
        else:
            target_coords = product_complex.coords

        return reactant_complex.modify(coords=target_coords)

    @classmethod
    def concenate_molecules(cls, mols):
        atoms = sum([m.atoms for m in mols], ())
        masses = np.concatenate([m.masses for m in mols])
        coords = np.concatenate(
            [
                m.coords for m in mols
            ],
            axis=0
        )
        mol_shift = np.cumsum([0] + [len(m.atoms) for m in mols])
        bonds = sum(
            [
                [
                    [b[0] + shift, b[1] + shift] + list(b[2:])
                    for b in mol.bonds
                ]
                for mol, shift in zip(mols, mol_shift)
            ],
            []
        )

        return Molecule(
            atoms,
            coords,
            masses=masses,
            bonds=bonds
        )

    @classmethod
    def make_reaction_complex(cls,
                              mols,
                              coords=None,
                              fragment_expansion_method=None,
                              **expansion_method_options
                              ):
        reactant_complex = cls.concenate_molecules(mols)
        if coords is not None:
            reactant_complex = reactant_complex.modify(coords=coords)
        if fragment_expansion_method is not None:
            fragment_indices = reactant_complex.fragment_indices
            method, opts = cls.get_expansion_dispatch().resolve(fragment_expansion_method)
            target_coords = method(fragment_indices, reactant_complex.coords, **dict(opts, **expansion_method_options))
            reactant_complex = reactant_complex.modify(coords=target_coords)
        return reactant_complex

    @classmethod
    def from_smiles(cls, smiles, fragment_expansion_method=dev.default, align=True,
                    add_implicit_hydrogens=True,
                    **opts):
        reactants, products = [s.strip() for s in smiles.split(">>")]
        reactants = [
            Molecule.from_string(
                s.strip(), "smi", add_implicit_hydrogens=add_implicit_hydrogens
            )
            for s in reactants.split(".")
        ]

        products = [
            Molecule.from_string(
                s.strip(), "smi", add_implicit_hydrogens=add_implicit_hydrogens
            )
            for s in products.split(".")
        ]

        return cls(reactants, products,
                   fragment_expansion_method=fragment_expansion_method,
                   align=align,
                   **opts
                   )

    def get_profile_generator(self, profile_generator=None, internals=None, **profile_opts) -> ProfileGenerator:
        if profile_generator is None:
            profile_generator = self.profile_generator
        return ProfileGenerator.resolve(
            self.reactant_complex,
            self.product_complex,
            profile_generator=profile_generator,
            internals=internals,
            **profile_opts
        )
