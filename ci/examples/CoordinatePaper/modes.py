from Psience.Molecools import Molecule
from Psience.Modes import MixtureModes, NormalModes, LocalizedModes, ObliqueModeGenerator
import McUtils.Coordinerds as coordops
import McUtils.Numputils as nput
import McUtils.Formatters as mfmt
import McUtils.Devutils as dev
from McUtils.Data import UnitsData
import numpy as np

from . import paths

ZERO_FREQ_CUTOFF = 25 * UnitsData.convert("Wavenumbers", "Hartrees")
RPMODE_GRADIENT_CUTOFF = 50 / UnitsData.hartrees_to_wavenumbers
def prep_rp_modes(me_ints,
                  nms=None,
                  potential_derivatives=None,
                  internals=None,
                  project_transrot=False,
                  proj_coord=(3, 2, 1, 0),
                  try_reaction_path=True,
                  return_status=False
                  ) -> ((LocalizedModes, NormalModes,NormalModes), bool):
    if potential_derivatives is not None:
        me_ints = me_ints.modify(potential_derivatives=[np.array(v) for v in potential_derivatives])
    if nms is None:
        nms = me_ints.get_normal_modes(
            use_internals=False,
            project_transrot=project_transrot,
            # potential_derivatives=me_ints.potential_derivatives,
            zero_freq_cutoff=ZERO_FREQ_CUTOFF
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
        zero_freq_cutoff=ZERO_FREQ_CUTOFF
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
        tf = me_ints.get_internals_by_cartesians(1, strip_embedding=True)
        inv = me_ints.get_cartesians_by_internals(1, strip_embedding=True)
        nms_int = type(nms)(
            me_ints.internal_coordinates.system,
            inv[0] @ nms.remove_mass_weighting().modes_by_coords,
            freqs=nms.freqs,
            inverse=nms.remove_mass_weighting().coords_by_modes @ tf[0],
            g_matrix=me_ints.get_gmatrix(use_internals=True),
            masses=me_ints.atomic_masses,
            origin=me_ints.get_internals(strip_embedding=True),
            mass_weighted=False
        )
        # nms_int = me_ints.get_normal_modes(
        #     use_internals=True,
        #     project_transrot=False,
        #     zero_freq_cutoff=ZERO_FREQ_CUTOFF
        # )
        print(nms_int.freqs * UnitsData.hartrees_to_wavenumbers)
        print(nms_int.compute_freqs() * UnitsData.hartrees_to_wavenumbers)
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
        # print(locs.compute_freqs() * UnitsData.hartrees_to_wavenumbers)
        # with np.printoptions(linewidth=1e8):
        #     print(np.round(locs.localizing_transformation[0], 8)**2)
        #     print(locs.localizing_transformation[0] @ locs.localizing_transformation[1])
        # raise Exception(locs.freqs * UnitsData.hartrees_to_wavenumbers)

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
    from . import expansions
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

def get_oblique_rescaled_modes(cpmo3, modes, invs):
    g_red = nput.metric_tensor(modes, masses=cpmo3.atomic_masses)
    red_derivs = nput.tensor_reexpand(
        invs,
        cpmo3.potential_derivatives
    )
    f_red = red_derivs[1]
    f2, _, tf, inv = ObliqueModeGenerator(f_red, g_red).run()
    ob_inv = nput.tensor_reexpand([inv], invs)
    ob_modes = nput.tensor_reexpand(modes, [tf])
    return ob_inv, ob_modes


def get_oblique_diagonal_modes(cpmo3, modes, invs):
    g_red = nput.metric_tensor(modes, masses=cpmo3.atomic_masses)
    red_derivs = nput.tensor_reexpand(
        invs,
        cpmo3.potential_derivatives
    )
    f_red = red_derivs[1]
    f2, _, tf, inv = ObliqueModeGenerator(f_red, g_red).run()
    # print(f2 * UnitsData.hartrees_to_wavenumbers)
    a = np.diag(np.power(np.diag(f_red) / np.diag(g_red), 1 / 4))
    ai = np.diag(np.power(np.diag(g_red) / np.diag(f_red), 1 / 4))
    f_red2 = ai @ f_red @ ai
    b = np.diag(np.power(np.diag(f2) / np.diag(f_red2), 1/4))
    bi = np.diag(np.power(np.diag(f_red2) / np.diag(f2), 1/4))
    ob_inv = nput.tensor_reexpand([b @ ai], invs)
    ob_modes = nput.tensor_reexpand(modes, [a @ bi])
    return ob_modes, ob_inv

def get_local_rescaled_modes(cpmo3, modes, invs):
    g_red = nput.metric_tensor(modes, masses=cpmo3.atomic_masses)
    red_derivs = nput.tensor_reexpand(
        invs,
        cpmo3.potential_derivatives
    )
    f_red = red_derivs[1]
    a = np.diag(np.power(np.diag(f_red) / np.diag(g_red), 1 / 4))
    ai = np.diag(np.power(np.diag(g_red) / np.diag(f_red), 1 / 4))
    ob_inv = nput.tensor_reexpand([ai], invs)
    ob_modes = nput.tensor_reexpand(modes, [a])
    return ob_modes, ob_inv

def get_cpmo3_symmetry_coordinates(cpmo3:Molecule, rescale=True):
    from . import expansions
    cpmo3_int = cpmo3.modify(internals=expansions.cpmo3_zmatrix)
    pg, coeffs, intenrals, red_exps, base_exps = cpmo3_int.get_symmetrized_internals(
        atom_selection=[0, 1, 2, 3, 4, 10], tol=.8, permutation_tol=.9,
        return_expansions=True,
        return_point_group=True,
        reduce_redundant_coordinates={'untransformed_coordinates': [0, 1, 3, 8, 15]},
        extra_internals=[(0, 10), (0, 1)]
    )

    tf_data = (pg, coeffs, intenrals, red_exps[0], base_exps)

    expansion_data = red_exps[1]
    if dev.str_is(rescale, 'local'):
        expansion_data = get_local_rescaled_modes(cpmo3, *expansion_data)
    elif dev.str_is(rescale, 'oblique'):
        expansion_data = get_oblique_rescaled_modes(cpmo3, *expansion_data)
    elif rescale:
        expansion_data = get_oblique_diagonal_modes(cpmo3, *expansion_data)
    return expansion_data, tf_data

def get_contraction_modes(cpmo3, rescale=True):
    from . import expansions
    cpmo3_int = cpmo3.modify(internals=expansions.cpmo3_zmatrix)

    centroid = np.average(cpmo3.coords[:5], axis=0)
    axis = centroid - cpmo3.coords[10]
    g12 = cpmo3.get_gmatrix(use_internals=False, power=1/2)
    gi12 = cpmo3.get_gmatrix(use_internals=False, power=-1/2)
    extra_coord = np.concatenate([
        np.repeat(-axis[np.newaxis], 10, axis=0),
        np.repeat(axis[np.newaxis], 7, axis=0)
    ], axis=0).flatten() @ gi12

    tf = nput.maximum_similarity_transformation(
        gi12 @ cpmo3_int.get_cartesians_by_internals(1, strip_embedding=True)[0].T,
        extra_coord[:, np.newaxis],
        apply_transformation=False
    )
    full_basis = nput.find_basis(np.concatenate([tf, np.eye(tf.shape[0])], axis=1), method='qr')
    mode = nput.tensor_reexpand(cpmo3_int.get_internals_by_cartesians(1, strip_embedding=True), [full_basis])
    inv = nput.tensor_reexpand([full_basis.T], cpmo3_int.get_cartesians_by_internals(1, strip_embedding=True))

    expansion_data = (mode, inv)
    if dev.str_is(rescale, 'local'):
        expansion_data = get_local_rescaled_modes(cpmo3, *expansion_data)
    elif dev.str_is(rescale, 'oblique'):
        expansion_data = get_oblique_rescaled_modes(cpmo3, *expansion_data)
    elif rescale:
        expansion_data = get_oblique_diagonal_modes(cpmo3, *expansion_data)
    mode, inv = expansion_data

    return (mode, inv), (full_basis, coordops.extract_zmatrix_internals(expansions.cpmo3_zmatrix))

def get_cpmo3_modes(struct_id, projection_type='symmetrized', rescale=True, localize=False, return_structure=False):
    from . import expansions
    cpmo3 = expansions.load_cpmo3(struct_id)

    if projection_type == 'symmetrized':
        modes = get_cpmo3_symmetry_coordinates(cpmo3, rescale=rescale)
    else:
        modes = get_contraction_modes(cpmo3, rescale=rescale)

    if localize:
        modes = cpmo3.get_normal_modes().localize(target_modes=modes[0][0])

    if return_structure:
        return cpmo3, modes
    else:
        return modes



