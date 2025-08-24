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
                  try_reaction_path=True,
                  return_status=False
                  ) -> ((LocalizedModes, NormalModes,NormalModes), bool):
    if nms is None:
        nms = me_ints.get_normal_modes(
            use_internals=False,
            project_transrot=project_transrot,
            potential_derivatives=me_ints.potential_derivatives
        )

    # rpnms, stat = nms.get_reaction_path_modes(
    #     me_ints.potential_derivatives[0],
    #     return_status=True,
    #     zero_gradient_cutoff=RPMODE_GRADIENT_CUTOFF,
    #     zero_freq_cutoff=50 * UnitsData.convert("Wavenumbers", "Hartrees")
    # )
    rpnms, stat = me_ints.get_reaction_path_modes(
        use_internals=False,
        return_status=True,
        zero_gradient_cutoff=RPMODE_GRADIENT_CUTOFF,
        project_transrot=project_transrot,
        zero_freq_cutoff=30 * UnitsData.convert("Wavenumbers", "Hartrees")
    )

    if try_reaction_path and stat:
        rp2 = rpnms.make_mass_weighted()
        nm2 = nms.make_mass_weighted()
        tf = nm2.coords_by_modes @ rp2.modes_by_coords[:, 1:]
        inv = rp2.coords_by_modes[1:, :] @ nm2.modes_by_coords
        locs = nms.apply_transformation(
            (tf, inv)
        )
        freqs, m2, i2 = NormalModes.get_normal_modes(
            locs.compute_hessian(),
            locs.compute_gmatrix()
        )
        tf = tf @ m2
        inv = i2 @ inv
        locs = nms.apply_transformation(
            (tf, inv)
        )
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

    # print(
    #     mfmt.TableFormatter("{:.1f}").format(
    #         nput.tensor_reexpand(
    #             [nms.make_frequency_scaled().coords_by_modes],
    #             [0, me_ints.potential_derivatives[1]]
    #         )[1] * UnitsData.hartrees_to_wavenumbers
    #     )
    # )
    # print("-"*30)
    # print(
    #     mfmt.TableFormatter("{:.1f}").format(
    #         nput.tensor_reexpand(
    #             [rpnms.make_frequency_scaled().coords_by_modes],
    #             [0, me_ints.potential_derivatives[1]]
    #         )[1] * UnitsData.hartrees_to_wavenumbers
    #     )
    # )
    #
    # print("-"*30)
    # print(
    #     mfmt.TableFormatter("{:.1f}").format(
    #         nput.tensor_reexpand(
    #             [locs.make_frequency_scaled().coords_by_modes],
    #             [0, me_ints.potential_derivatives[1]]
    #         )[1] * UnitsData.hartrees_to_wavenumbers
    #     )
    # )
    #
    # raise Exception(...)

    # print(np.round(
    #     locs.make_mass_weighted().coords_by_modes @
    #     rpnms.make_mass_weighted().modes_by_coords[:, 1:]
    # ))
    # print(locs.coords_by_modes.shape, rpnms.coords_by_modes.shape)
    # raise ValueError(...)
    # print("=" * 75)
    # print(
    #     mfmt.TableFormatter("{:.3f}").format(
    #         nms_int.make_mass_weighted().modes_by_coords
    #         @ np.diag(nms_int.freqs ** 2)
    #         @ nms_int.make_mass_weighted().modes_by_coords.T * 219474.63
    #     )
    # )
    # print("=" * 75)
    # print(
    #     mfmt.TableFormatter("{:.3f}").format(
    #         locs.make_mass_weighted().modes_by_coords
    #         @ np.diag(locs.freqs ** 2)
    #         @ locs.make_mass_weighted().modes_by_coords.T * 219474.63
    #     )
    # )
    # print("=" * 75)

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
        locs.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
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