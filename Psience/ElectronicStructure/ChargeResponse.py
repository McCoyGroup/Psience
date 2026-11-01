"""
Dipoles and polarizabilities of charge-equilibration (fluctuating-charge) models from the
derivatives of their electrostatic energy with respect to the atomic charges
"""

import numpy as np

import McUtils.Numputils as nput

__all__ = [
    "FluctuatingChargeResponse",
]


class FluctuatingChargeResponse:
    """
    For a model whose charges minimize `E(q) = chi . q + q . J q / 2` subject to a fixed total
    charge, the induced dipole and the dipole polarizability follow from the electronegativities
    `chi = dE/dq` and hardness matrix `J = d^2E/dq^2`. Given those (and their nuclear-coordinate
    derivatives), `expand` returns the coordinate expansions of both.
    """

    @classmethod
    def constrained_polarizability(cls, coords, hardness, hardness_derivatives=None,
                                  *, allow_indefinite=False, return_diagnostics=False,
                                  stability_tolerance=1e-10, solver="constrained"):
        """Charge-conserving static response and its first Cartesian derivative.

        ``hardness`` has shape (..., N, N), in energy/charge**2. Coordinate
        derivatives have shape (..., 3*N, N, N). Outputs are in the corresponding
        charge**2 * distance**2 / energy units. No hardness regularization is
        applied. ``solver="inverse"`` instead projects a full inverse and needs
        a nonsingular full hardness. ``allow_indefinite`` permits evaluation of unstable
        response models; singular constrained systems always raise.
        """
        if solver not in ("constrained", "inverse"):
            raise ValueError(f"unknown charge-response solver {solver!r}")
        if not np.isfinite(stability_tolerance) or stability_tolerance < 0:
            raise ValueError("stability_tolerance must be finite and nonnegative")
        coords = np.asarray(coords, dtype=float)
        hardness = np.asarray(hardness, dtype=float)
        n = coords.shape[-2]
        if coords.shape[-1] != 3 or hardness.shape[-2:] != (n, n):
            raise ValueError("expected coordinates (..., N, 3) and hardness (..., N, N)")
        batch = np.broadcast_shapes(coords.shape[:-2], hardness.shape[:-2])
        coords = np.broadcast_to(coords, batch + (n, 3))
        hardness = np.broadcast_to(hardness, batch + (n, n))
        if not np.all(np.isfinite(coords)) or not np.all(np.isfinite(hardness)):
            raise ValueError("coordinates and hardness must be finite")
        if not np.allclose(hardness, np.swapaxes(hardness, -1, -2), atol=1e-10, rtol=1e-8):
            raise ValueError("charge hardness must be symmetric")
        hardness = (hardness + np.swapaxes(hardness, -1, -2))/2
        # This basis removes only changes in total charge, not a spatial direction.
        Z = np.linalg.svd(np.ones((1, n)), full_matrices=True)[2][1:].T
        K = np.einsum('na,...nm,mb->...ab', Z, hardness, Z)
        D = np.einsum('na,...nj->...aj', Z, coords-coords.mean(axis=-2, keepdims=True))
        eigenvalues = np.linalg.eigvalsh(K)
        if n > 1:
            scale = np.maximum(1., np.max(np.abs(eigenvalues), axis=-1))
            if np.any(np.min(np.abs(eigenvalues), axis=-1) <= stability_tolerance*scale):
                raise np.linalg.LinAlgError("singular or near-singular constrained charge hardness")
            if not allow_indefinite and np.any(eigenvalues <= stability_tolerance*scale[..., None]):
                raise ValueError("unstable constrained charge hardness; diagnostic evaluation requires allow_indefinite=True")
        Y = np.linalg.solve(K, D)
        X = np.einsum('na,...aj->...nj', Z, Y)
        if solver == "inverse" and n > 1:
            # Compatibility route requiring a nonsingular full hardness. Its
            # neutral response equals the null-space solve when both exist.
            A = np.linalg.inv(hardness)
            a = A.sum(axis=-1)
            denominator = a.sum(axis=-1)
            if np.any(np.abs(denominator) <= stability_tolerance*np.maximum(1., np.sum(np.abs(A), axis=(-2,-1)))):
                raise np.linalg.LinAlgError("singular full-inverse charge constraint")
            P = A-a[..., :, None]*a[..., None, :]/denominator[..., None, None]
            X = P@(coords-coords.mean(axis=-2,keepdims=True))
            Y = np.einsum('na,...nj->...aj',Z,X)
        alpha = np.einsum('...ai,...aj->...ij', D, Y)
        expansion = [(alpha + np.swapaxes(alpha, -1, -2))/2]
        if hardness_derivatives is not None:
            dJ = np.broadcast_to(np.asarray(hardness_derivatives, dtype=float), batch+(3*n,n,n))
            if not np.all(np.isfinite(dJ)):
                raise ValueError("hardness derivatives must be finite")
            dx = np.eye(3*n).reshape(3*n, n, 3)
            dD = np.einsum('na,knj->kaj', Z, dx)
            dK = np.einsum('na,...knm,mb->...kab', Z, dJ, Z)
            cross = np.einsum('kai,...aj->...kij', dD, Y)
            dalpha = cross + np.swapaxes(cross,-1,-2) - np.einsum('...ai,...kab,...bj->...kij', Y,dK,Y)
            expansion.append((dalpha + np.swapaxes(dalpha,-1,-2))/2)
        if return_diagnostics:
            diagnostics = dict(constrained_eigenvalues=eigenvalues,
                               stable=np.all(eigenvalues>stability_tolerance*np.maximum(1.,np.max(np.abs(eigenvalues),axis=-1,initial=0))[...,None],axis=-1),
                               solver=solver,
                               solve_residual=np.linalg.norm(K@Y-D,axis=(-2,-1)),
                               charge_response=X)
            return expansion, diagnostics
        return expansion


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

    @classmethod
    def expand(cls, total_charge, coords, potential_by_charges_by_coords):
        """
        :param total_charge: the total charge
        :param coords: nuclear coordinates, `(..., natoms, 3)`
        :param potential_by_charges_by_coords: `[[E, dE/dx, ...], [chi, dchi/dx, ...], [J, dJ/dx, ...]]`
            (charge-derivative order outermost, coordinate-derivative order innermost)
        :return: `[dipole_expansion, polarizability_expansion]` (each a list over coordinate orders),
            truncated to the available charge-derivative orders
        :rtype: list[list[np.ndarray]]
        """
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

            J_inv_chi_expansion = nput.tensordot_deriv(chi_expansion, J_inv_expansion, order=order, axes=[-1, -1])
            J_inv_chi_sum_expansion = nput.tensordot_deriv(J_inv_sum_1, chi_expansion, order=order, axes=[-1, -1])
            mu_expansion = nput.scalarprod_deriv(J_inv_chi_sum_expansion, J_inv_tot_inv, order=order)
            mu_expansion[0] += total_charge * J_inv_tot_inv[0]

            mu_prod_expansion = nput.scalarprod_deriv(mu_expansion, J_inv_sum_1, order=order)
            dipole_expansion_charge_factor = nput.add_expansions(J_inv_chi_expansion, mu_prod_expansion)
            dipole_expansion = nput.tensordot_deriv(dipole_expansion_charge_factor, coord_expansion, axes=[-1, -2], order=order)
            expansion.append(dipole_expansion)

            if es_order > 1:
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
    def distorted_harmonic_polarizabilities(cls, charges, dipole_expansion, coords, hessians):
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
