import collections, functools, numpy as np
# might be worth excising theses since this could be standalone otherwise?
from McUtils.Data import AtomData, UnitsData
from McUtils.Zachary import FiniteDifferenceDerivative, RBFDInterpolator
import McUtils.Numputils as nput

from ..Molecools import Molecule, MolecularZMatrixCoordinateSystem
from ..Molecools.Properties import PropertyManager

__all__ = [
    "AIMDSimulator",
    "PairwisePotential"
]

class AIMDSimulator:
    def __init__(self,
                 atoms,
                 coordinates,
                 force_function,
                 atomic_structures=None,
                 internals=None,
                 velocities=0,
                 track_kinetic_energy=False,
                 track_velocities=False,
                 timestep=.1,
                 sampling_rate=1
                 ):
        self.masses = np.array([
                AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
            if isinstance(a, str) else a for a in atoms
        ])
        coords = np.asanyarray(coordinates)
        if coords.ndim == 1:
            coords = coords[np.newaxis]

        if atomic_structures is None:
            atomic_structures = (
                coords.shape[-1] == 3
                    and
                coords.shape[-2] == len(self.masses)
            )
        self._atomic_structs = atomic_structures
        if self._atomic_structs:
            if coords.ndim == 2: coords = coords[np.newaxis]
            self._mass = np.expand_dims(self.masses, -1)
            for _ in range(coords.ndim - 2):
                self._mass = np.expand_dims(self._mass, 0)
            # self._mass = self.masses[np.newaxis, :,  np.newaxis]
        else: # regular coordinates
            self._mass = self.masses[np.newaxis]
            for _ in range(coords.ndim - 2):
                self._mass = np.expand_dims(self._mass, 0)
        self.coords = coords
        # for _ in range(coordinates.ndim - 2):
        #     self._mass = np.expand_dims(self._mass, -1)
        if isinstance(velocities, (int, float, np.integer, np.floating)):
            velocities = np.full_like(self.coords, velocities)
        velocities = np.asanyarray(velocities)
        self.velocities = velocities
        if internals is not None:
            raise NotImplementedError("haven't fully managed getting embedding right yet...")
            if not isinstance(internals, MolecularZMatrixCoordinateSystem):
                base_mol = Molecule(
                    ['H']*len(self.masses),
                    coords=self.coords,
                    masses=self.masses,
                    internals=internals
                )
                internals = [base_mol.coords.system, base_mol.internal_coordinates.system]
            else:
                internals = [internals.molecule.coords.system, internals]
        self._internals = internals
        self.force_function = force_function
        self._prev_forces = None

        self.trajectory = collections.deque()
        self.trajectory.append(self.coords)
        if track_kinetic_energy:
            self.kinetic_energies = collections.deque()
            self.kinetic_energies.append(
                self._ke(self._mass, self.velocities)
            )
        else:
            self.kinetic_energies = None
        if track_velocities:
            self.velocity_deque = collections.deque()
            self.velocity_deque.append(self.velocities)
        else:
            self.velocity_deque = None
        self.steps = 0
        self.dt = timestep
        self.sampling_rate = sampling_rate

    @classmethod
    def mode_energies_to_velocities(cls, modes, masses, energy_splits, inverse=None):

        masses = np.asanyarray(masses)
        e_part = np.asanyarray(energy_splits)
        smol = e_part.ndim == 1
        if smol: e_part = e_part[np.newaxis]

        # if not np.allclose(np.linalg.norm(modes, axis=0), np.ones(modes.shape[1])):  # i.e. came out of a generalized eigenvalue run
        #     modes = modes.reshape(-1, 3, modes.shape[1]) * np.sqrt(masses[:, np.newaxis, np.newaxis])
        #     modes = modes.reshape(-1, modes.shape[-1])

        if inverse is None:
            if not np.allclose(np.linalg.norm(modes, axis=0), np.ones(modes.shape[1])):
                raise ValueError("modes aren't normalized -> non-unitary transformation")
            inverse = modes.T

        tripmass = np.broadcast_to(masses[:, np.newaxis], (len(masses), 3)).flatten()
        g = inverse @ np.diag(tripmass) @ modes

        gv, Q = np.linalg.eigh(g)
        sorting = np.argsort(np.argmax(Q ** 2, axis=0))
        gv = gv[sorting]  # maximum similarity to OG vectors
        Q = Q[:, sorting]

        axes = (Q.T @ inverse).reshape(modes.shape[1], -1, 3)
        v_part = np.sign(e_part) * np.sqrt(2 * np.abs(e_part)) # / gv[np.newaxis])
        vels = np.sum(
            axes[np.newaxis] *
                v_part[:, :, np.newaxis, np.newaxis],
            axis=1
        )

        return vels

    #TODO: add internal coordinate propagation following the Krimm paper using the B-matrix
    #      or not...maybe I just need a faster Jacobian?
    def get_forces(self, coords):
        if self._internals is not None:
            ccs, ics = self._internals
            icoords = ccs.convert_coords(coords, ics)
            if not isinstance(icoords, np.ndarray):
                icoords = icoords[0] # converter info
            forces = self.force_function(icoords)
            jacs = ccs(coords).jacobian(ics, order=1, all_numerical=True)[0] #dRdX?
            ncs = np.prod(coords.shape[-2:], dtype=int)
            jacs = np.reshape(np.moveaxis(jacs, 1, 0), (len(coords), ncs, ncs))
            forces = jacs@forces.reshape((len(coords), ncs, 1))
            forces = forces.reshape(coords.shape)
            # raise Exception(forces.shape)
        else:
            forces = self.force_function(coords)
        return forces

    def step(self):
        forces = self._prev_forces
        if forces is None:
            forces = self.get_forces(self.coords)


        v = self.velocities
        coords = self.coords + v * self.dt + forces / (2 * self._mass) * self.dt**2
        if self._atomic_structs:
            com = np.tensordot(self.masses, coords, axes=[0, -2]) / np.sum(self._mass)
            coords = coords - com[:, np.newaxis, :] # don't let COM move
        forces_new = self.get_forces(coords)
        vels = v + self.dt * (forces + forces_new) / (2 * self._mass)

        self._prev_forces = forces_new
        self.velocities = vels
        self.coords = coords
        self.steps += 1

        return coords, vels, forces_new

    @staticmethod
    def _ke(_mass, v):
        # _mass and v are 2D b.c. atomic structs
        return 1 / 2 * np.sum(np.sum(_mass * v ** 2, axis=-1), axis=-1)
    def propagate(self, num_steps=1):

        for _ in range(num_steps):
            c, v, f = self.step()
            if self.steps % self.sampling_rate == 0:
                if self.kinetic_energies is not None:
                    self.kinetic_energies.append(self._ke(self._mass, v))
                self.trajectory.append(c)
                if self.velocity_deque is not None:
                    self.velocity_deque.append(v)

        return self.trajectory

    def build_interpolation(self, energy_function, interpolation_order=2,
                            equilibration_steps=None,
                            interpolator_class=None, eckart_embed=True, reference=None, **interpolator_options):

        traj = np.array(self.trajectory)
        if equilibration_steps is not None:
            traj = traj[equilibration_steps:]
        traj = traj.reshape((-1,) + traj.shape[-2:])

        vals = [energy_function(traj)]

        # if clustering_criterion > ...:
        #     ...

        # TODO: add clustering radius + energy cutoff within disks

        if interpolation_order > 0:
            vals.append(-self.force_function(traj))
        if interpolation_order > 1:
            cshape = (-1,) + traj.shape[-2:]
            npts = np.prod(traj.shape[-2:], dtype=int)
            grad = lambda x:-self.force_function(x.reshape(cshape)).reshape(x.shape)
            hess_fun = FiniteDifferenceDerivative(grad, function_shape=(npts, npts))
            hess = np.moveaxis(
                hess_fun.derivatives(traj.reshape(-1, npts)).derivative_tensor([1])[0],
                1, 0
            )
            vals.append(hess)
        if interpolation_order > 2:
            raise ValueError("can only do order 2 interps for now because I am laaazy")

        if eckart_embed:
            if isinstance(reference, Molecule):
                ref = reference
            elif reference is None:
                ref_pos = np.argmin(vals[0])
                ref = Molecule(
                    ["H"] * len(self.masses),
                    traj[ref_pos],
                    masses=self.masses,
                )
            else:
                ref = Molecule(
                    ["H"] * len(self.masses),
                    reference,
                    masses=self.masses,
                )
            ref = ref.get_embedded_molecule(load_properties=False)

            # all_crds = []
            # all_ders = [[] for _ in range(interpolation_order)]

            rots, _, (pax_traj, _, pax_rots) = ref.get_embedding_data(traj)
            traj = pax_traj @ np.swapaxes(rots, -2, -1)

            if interpolation_order > 0:
                npts = np.prod(traj.shape[-2:], dtype=int)
                cshape = traj.shape[-2:]

                # new_derivs = [vals[1] @ pax_rots @ np.swapaxes(rots, -2, -1)]
                # raise Exception(new_derivs[0])
                new_derivs = PropertyManager._transform_derivatives(
                    [v.reshape((-1,) + (npts,)*(i+1)) for i,v in enumerate(vals[1:])],
                    rots @ pax_rots
                )
                vals = [vals[0]] + [v.reshape((-1,) + cshape*(i+1)) for i,v in enumerate(new_derivs)]


            # for crd_ders in zip(traj, *vals[1:]): #TODO: speed this up...
            #     c = crd_ders[0]
            #     mol = Molecule(
            #         ["H"] * len(self.masses),
            #         c,
            #         masses=self.masses,
            #         potential_derivatives=crd_ders[1:]
            #     ).get_embedded_molecule(ref, load_properties=False)
            #     all_crds.append(mol.coords)
            #     pes = mol.potential_derivatives
            #     for i, p in enumerate(pes):
            #         all_ders[i].append(p)

            # traj = np.array(all_crds)

        if interpolator_class is None:
            interpolator_class = RBFDInterpolator

        return interpolator_class(traj, *vals, **interpolator_options)

    def extract_trajectory(self, flatten=True, embed=True, extract_velocities=None):
        ref_struct = None
        if isinstance(embed, (np.ndarray, list, tuple)):
            ref_struct = np.asanyarray(embed)
            embed = True
        base_coords = np.array(self.trajectory)
        if flatten:
            base_coords = base_coords.reshape((-1,)+base_coords.shape[-2:])
        if extract_velocities is None:
            extract_velocities = self.velocity_deque is not None
        if extract_velocities:
            velocities = np.array(self.velocity_deque)
        else:
            velocities = None
        if embed:
            if ref_struct is None:
                if not flatten:
                    ref_crds = base_coords.reshape((-1,)+base_coords.shape[-2:])
                else:
                    ref_crds = base_coords
                if self.kinetic_energies is not None:
                    # assuming symplectic integration...
                    kes = np.array(self.kinetic_energies).flatten()
                    max_pos = np.argmax(kes)
                    ref_struct = base_coords[max_pos]
                else:
                    ref_struct = ref_crds[0]
            mol = Molecule(["H"]*len(self.masses), ref_struct, masses=self.masses)
            if extract_velocities:
                embedding = mol.get_embedding_data(base_coords)
                coords = np.tensordot(
                    nput.vec_tensordot(
                        embedding.coord_data.coords,
                        embedding.rotations,
                        shared=1,
                        axes=[-1, 2]
                    ),
                    embedding.reference_data.axes,
                    axes=[-1, 1]
                )
                velocities = np.tensordot(
                    nput.vec_tensordot(
                        np.reshape(velocities, coords.shape),  # dx/dt has same embedding
                        embedding.rotations,
                        shared=1,
                        axes=[-1, 2]
                    ),
                    embedding.reference_data.axes,
                    axes=[-1, 1]
                )
                coords = np.reshape(coords, base_coords.shape)
                velocities = np.reshape(velocities, base_coords.shape)
            else:
                coords = mol.embed_coords(base_coords)
        else:
            coords = base_coords
        if extract_velocities:
            return coords, velocities
        else:
            return coords



class PairwisePotential:
    def __init__(self, fun, deriv=None):
        self.fun = fun
        if deriv is None:
            if hasattr(deriv, 'deriv'):
                fun.deriv()
            else:
                deriv = self.numerical_deriv(fun)
        self.deriv = deriv
        self._inds = None # to cache indices when applied to same-shape data

    @classmethod
    def numerical_deriv(cls, fun, step_size=.01):
        @functools.wraps(fun)
        def deriv(dists):
            return (fun(dists + step_size) - fun(dists - step_size)) / (2 * step_size)
        return deriv

    def _get_inds_cached(self, coords):
        if self._inds is not None and len(coords) != self._inds[0]:
            self._inds = None
        if self._inds is None:
            row_inds, col_inds = np.triu_indices(len(coords), k=1)
            _, row_map = np.unique(row_inds, return_index=True) # One element less than the coords
            sorting = np.argsort(col_inds)
            _, col_map = np.unique(col_inds, return_index=True) # One element less than the coords
            self._inds = [len(coords), (row_inds, col_inds, row_map, sorting, col_map)]
        return self._inds[1]

    def eval(self, coords):
        row_inds, col_inds, row_map, col_sort, col_map = self._get_inds_cached(coords)
        pot = self.fun

        # compute unsigned force for each pair
        diffs = coords[row_inds] - coords[col_inds]
        dists = np.linalg.norm(diffs, axis=1)
        pot_vals = pot(dists)

        # print(row_map)

        # now reduce over different chunks for the positive contribs...
        pos_chunks = np.add.reduceat(
            pot_vals,
            row_map
        )
        neg_chunks = np.add.reduceat(
            pot_vals[col_sort],
            col_map
        )

        return np.concatenate(
            [
                pos_chunks[:1],
                pos_chunks[1:] + neg_chunks[:-1],
                pos_chunks[-1:]
            ],
            axis=0
        )

    def forces(self, coords):
        row_inds, col_inds, row_map, col_sort, col_map = self._get_inds_cached(coords)
        pot_deriv = self.deriv

        # compute unsigned force for each pair
        diffs = coords[row_inds] - coords[col_inds]
        dists = np.linalg.norm(diffs, axis=1)
        base_force = pot_deriv(dists)
        normals = diffs / dists[:, np.newaxis]
        force_list = -normals * base_force[:, np.newaxis]

        # now reduce over different chunks for the positive contribs...
        pos_chunks = np.add.reduceat(
            force_list,
            row_map
        )
        neg_chunks = np.add.reduceat(
            force_list[col_sort],
            col_map
        )

        return np.concatenate(
            [
                pos_chunks[:1],
                pos_chunks[1:] - neg_chunks[:-1],
                pos_chunks[-1:]
            ],
            axis=0
        )

    def __call__(self, coords):
        return self.eval(coords)
