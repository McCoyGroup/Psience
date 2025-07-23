from Psience.Molecools import Molecule
from Psience.Modes import NormalModes, LocalizedModes
import McUtils.Coordinerds as coordops
import McUtils.Numputils as nput
import McUtils.Formatters as mfmt
from McUtils.Data import UnitsData
import numpy as np

from . import expansions
from . import paths


RPMODE_GRADIENT_CUTOFF = 2e-5
def prep_rp_modes(me_ints, nms=None, internals=None,
                  project_transrot=False,
                  proj_coord=(3, 2, 1, 0),
                  return_status=False
                  ) -> ((LocalizedModes, NormalModes,NormalModes), bool):
    if nms is None:
        nms = me_ints.get_normal_modes(use_internals=False, project_transrot=project_transrot)
    rpnms, stat = me_ints.get_reaction_path_modes(
        use_internals=False,
        return_status=True,
        zero_gradient_cutoff=RPMODE_GRADIENT_CUTOFF
    )
    if stat:
        locs = nms.localize(target_modes=rpnms.modes_by_coords[:, 1:])
    else:
        if me_ints.internals is None:
            me_ints = me_ints.modify(internals=internals)
        nms_int = me_ints.get_normal_modes(use_internals=True, project_transrot=False)
        # print(
        #     mfmt.TableFormatter("{:.3f}").format(nms_int.modes_by_coords)
        # )
        locs = nms_int.localize(
            coordinate_constraints=coordops.zmatrix_indices(
                internals,
                [proj_coord]
                # [(3, 2, 1, 0)]
            ),
            maximum_similarity=True
        )

        print("=" * 75)
        print(
            mfmt.TableFormatter("{:.3f}").format(
                nms_int.make_mass_weighted().modes_by_coords
                @ np.diag(nms_int.freqs ** 2)
                @ nms_int.make_mass_weighted().modes_by_coords.T * 219474.63
            )
        )
        print("=" * 75)
        print(
            mfmt.TableFormatter("{:.3f}").format(
                locs.make_mass_weighted().modes_by_coords
                @ np.diag(locs.freqs ** 2)
                @ locs.make_mass_weighted().modes_by_coords.T * 219474.63
            )
        )
        print("=" * 75)

        # print("="*75)
        # print(
        #     mfmt.TableFormatter("{:.3f}").format(
        #         locs.localizing_transformation[1] @ locs.localizing_transformation[0]

        #         )
        # )
        # print("="*75)
        # print(
        #     mfmt.TableFormatter("{:.3f}").format(
        #         locs.make_mass_weighted().coords_by_modes @ nms_int.make_mass_weighted().modes_by_coords

        #     )
        # )
    if return_status:
        return (locs, nms, rpnms), stat
    else:
        return locs, nms, rpnms

def print_rpnm_hessian(me_ints, nms, rpnms, locs):
    # locs = nms.localize(coordinate_constraints=coordops.zmatrix_indices(
    #     methanol_zmatrix,
    #     [(3, 2, 1, 0)]
    # ))
    # print(mfmt.TableFormatter('{:.0f}').format(
    #     nms.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
    # print(
    #     mfmt.TableFormatter('{:.0f}').format(locs.local_hessian *
    #                                     UnitsData.hartrees_to_wavenumbers)
    # )
    # nms_carts = me_ints.get_normal_modes(project_transrot=True, use_internals=False)
    # locs = nms_carts.localize(internals=[(3, 2, 1, 0)], orthogonal_projection=True)
    # loc_2 = nms_carts.apply_transformation(locs.localizing_transformation).make_dimensionless()

    # cart_udim = nms_carts.make_dimensionless()

    locs_2 = nms.apply_transformation(locs.localizing_transformation).make_dimensionless()
    f_nmw = me_ints.potential_derivatives[1]
    g12 = nput.fractional_power(me_ints.get_gmatrix(use_internals=False), 1 / 2)
    f_mw = g12 @ f_nmw @ g12
    f_cart = nput.tensor_reexpand(
        [locs_2.coords_by_modes],
        [0, f_mw]
    )[-1]
    # g_cart = nms_carts.compute_gmatrix()
    # print(
    #     TableFormatter('{:.3f}').format(locs.localizing_transformation[1] @ locs.localizing_transformation[0])
    # )
    # f_loc = locs.localizing_transformation[1] @ f_cart @ locs.localizing_transformation[1].T
    # print(f_loc.shape)

    # break

    print(mfmt.TableFormatter('{:.0f}').format(
        rpnms.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
    print(
        mfmt.TableFormatter('{:.1f}').format(
            f_cart * (UnitsData.hartrees_to_wavenumbers)
        )
    )

def get_gaussian_modes(lot, subkey, use_internals=False, return_status=True, return_mol=False):
    me_gaussian = Molecule.from_file(
        paths.torsion_scan_path(lot, f'{subkey}.fchk'),
        internals=expansions.methanol_gaussian_zmatrix if use_internals else None
    )

    mode_data = prep_rp_modes(me_gaussian,
                  internals=expansions.methanol_gaussian_zmatrix if use_internals else None,
                  proj_coord=(3, 2, 1, 0),
                  return_status=return_status
                  )

    if return_mol:
        return me_gaussian, mode_data
    else:
        return mode_data