
import numpy as np
from .CoordinateSystems import MolecularEmbedding
from .Properties import NormalModesManager
from McUtils.Zachary import TensorDerivativeConverter

__all__ = [
    "MolecularEvaluator"
]

class MolecularEvaluator:
    def __init__(self, embedding:MolecularEmbedding, normal_modes:NormalModesManager):
        self.embedding = embedding
        self.normal_modes = normal_modes
    @property
    def coords(self):
        return self.embedding.coords

    @property
    def masses(self):
        return self.embedding.masses

    def evaluate(self,
                 func,
                 use_internals=None,
                 deriv_order=None,
                 strip_embedding=False
                 ):
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
            if deriv_order is None:
                return func(coords).view(np.ndarray)

            terms = func(coords, deriv_order=deriv_order)
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
            jacs = self.embedding.get_internals_by_cartesians(deriv_order)

            terms = TensorDerivativeConverter(
                jacs,
                terms[1:],
                jacobians_name='dXdR',
                values_name='f'
            ).convert()  # , check_arrays=True)

            return [const.view(np.ndarray)] + [t.view(np.ndarray) for t in terms]
        else:
            if deriv_order is None:
                return func(self.embedding.coords).view(np.ndarray)
            else:
                return [x.view(np.ndarray) for x in func(self.embedding.coords, deriv_order=deriv_order)]

    def evaluate_at(self,
                    func,
                    coords,
                    use_internals=None,
                    deriv_order=None,
                    strip_embedding=False
                    ):
        return type(self)(
            MolecularEmbedding(self.embedding.masses, coords, self.embedding.internals),
            self.normal_modes
        ).evaluate(func, use_internals=use_internals, deriv_order=deriv_order, strip_embedding=strip_embedding)

    def get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None,
                                  use_internals=False,
                                  strip_embedding=False,
                                  shift=True
                                  ):
        displacements = np.asanyarray(displacements)

        if which is not None:
            which = tuple(
                np.ravel_multi_index(idx, (len(self.embedding.masses), 3))
                    if not isinstance(idx, (int, np.integer)) else
                idx
                for idx in which
            )

        if use_internals and self.embedding.internals is None:
            raise ValueError("can't displace in internals without internal coordinate spec")
        base_coords = self.embedding.coords if not use_internals else self.embedding.internal_coordinates
        if strip_embedding:
            ecs = self.embedding.embedding_coords
            all_coords = np.arange(len(self.embedding.masses) * 3)
            which = np.setdiff1d(all_coords, ecs)[which,]

        if which is not None:
            if displacements.shape[-1] != len(which):  # displacements provided in atom coordinates
                displacements = displacements.reshape(
                    displacements.shape[:-2] +
                    (np.prod(displacements.shape[-2:], dtype=int),)
                )
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
            if displacements.shape[-2:] != base_coords.shape:
                raise ValueError("displacements with shape {} passed but coordinates have shape {}".format(
                    displacements.shape,
                    base_coords.shape
                ))
            base_coords = displacements

        if use_internals:
            # track the embedding info...
            base_coords = self.embedding.internal_coordinates.system(base_coords, **self.embedding.internal_coordinates.converter_options)
            if isinstance(use_internals, str):
                if use_internals == 'convert':
                    base_coords = base_coords.convert(self.embedding.coords.system)
                elif use_internals == 'reembed':
                    base_coords = self.embedding.embed_coords(base_coords.convert(self.embedding.coords.system))
        else:
            base_coords = self.embedding.coords.system(base_coords)
        return base_coords

    def get_scan_coordinates(self,
                             domains,
                             internals=False,
                             which=None, sel=None, axes=None,
                             shift=True
                             ):

        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        return self.get_displaced_coordinates(displacement_mesh, shift=shift,
                                              use_internals=internals, which=which, sel=sel, axes=axes)

    def get_nearest_displacement_atoms(self,
                                       points,
                                       sel=None, axes=None, weighting_function=None,
                                       return_distances=False
                                       ):

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
        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        return self.get_nearest_displacement_coordinates(displacement_mesh, axes=axes, sel=sel)