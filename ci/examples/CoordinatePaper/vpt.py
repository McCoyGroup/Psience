import os
import enum
import numpy as np
from McUtils.Data import UnitsData
from Psience.Molecools import Molecule
from Psience.VPT2 import PerturbationTheoryHamiltonian

from . import paths
from .expansions import (
    get_aimnet_structure, get_aimnet_expansion,
    get_mace_structure, get_mace_expansion,
    methanol, methanol_zmatrix, methanol_gaussian_zmatrix
)
from .modes import prep_rp_modes, print_rpnm_hessian


class LevelsOfTheory(enum.Enum):
    B3LYP = 'b3lyp'
    wB97 = 'wb97'
    AIMNet2 = 'aimnet'
    AIMNet2Old = 'aimnet-old'
    MACE = 'mace'

def get_gaussian_logfile(lot, subkey, use_reaction_path=True, use_degeneracies=True, mode_selection=None, use_internals=False):
    ms = "".join(str(m) for m in
                 (mode_selection if mode_selection is not None else 'all')
                 )
    if not use_degeneracies:
        deg = "_nodegs"
    else:
        deg = ""

    if use_reaction_path:
        if use_internals:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_rp_ints_{subkey}_{ms}{deg}.out')
        else:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_rp_{subkey}_{ms}{deg}.out')
    else:
        if use_internals:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_norp_ints_{subkey}_{ms}{deg}.out')
        else:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_norp_{subkey}_{ms}{deg}.out')

    return log_file

def get_wb97_logfile(subkey, **opts):
    return get_gaussian_logfile('wb97', subkey, **opts)

def get_b3lyp_logfile(subkey, **opts):
    return get_gaussian_logfile('b3lyp', subkey, **opts)

def run_gaussian_vpt(
        lot,
        subkey,
        file_pattern='methanol_vpt_{key}.fchk',
        use_internals=False,
        use_reaction_path=True,
        mode_selection=None,
        use_degeneracies=True,
        log_file=None,
        return_runner=False,
        overwrite=False
):
    me_gaussian = Molecule.from_file(
        paths.torsion_scan_path(lot, file_pattern.format(key=subkey)),
        internals=methanol_gaussian_zmatrix if use_internals else None
    )

    # fp = file_pattern.format(i=i + 1,
    #                          modes="".join(str(m) for m in
    #                                        (mode_selection if mode_selection is not None else 'all')
    #                                        )
    #                          )

    if use_reaction_path:
        (locs, nms, rpnms), stat = prep_rp_modes(me_gaussian,
                                                 internals=methanol_gaussian_zmatrix,
                                                 proj_coord=(3, 2, 1, 0),
                                                 return_status=True
                                                 )
    else:
        stat = False
        locs = None
        rpnms = None
        nms = me_gaussian.get_normal_modes(project_transrot=False, use_internals=False).remove_mass_weighting().remove_frequency_scaling()

    os.makedirs(paths.torsion_scan_path(lot, 'results'), exist_ok=True)
    if log_file is None:
        log_file = get_gaussian_logfile(lot, subkey,
                                        use_reaction_path=(stat or use_reaction_path),
                                        mode_selection=mode_selection,
                                        use_degeneracies=use_degeneracies,
                                        use_internals=use_internals)

    if stat or use_reaction_path:
        tf = locs.localizing_transformation
        if mode_selection is not None:
            tf = locs[mode_selection].localizing_transformation

        print_rpnm_hessian(me_gaussian, nms, rpnms, locs)

        mode_selection = None
    else:
        tf = None

    if overwrite or not os.path.exists(log_file):
        try:
            os.remove(log_file)
        except:
            ...
        runner, _ = me_gaussian.setup_VPT(states=2,
                                          degeneracy_specs='auto' if use_degeneracies else None,
                                          cartesian_analytic_deriv_order=-1,
                                          cartesian_by_internal_derivative_method='fast',
                                          # gmatrix_tolerance=1*UnitsData.convert("Wavenumbers", "Hartrees"),
                                          modes=nms,
                                          mode_transformation=tf,
                                          mode_selection=mode_selection,
                                          logger=log_file
                                          )
        if return_runner:
            return runner
        try:
            runner.print_tables(logger=runner.hamiltonian.logger)
        except:
            os.remove(log_file)
            raise
    elif return_runner:
        runner, _ = me_gaussian.setup_VPT(states=2,
                                          degeneracy_specs='auto' if use_degeneracies else None,
                                          cartesian_analytic_deriv_order=-1,
                                          cartesian_by_internal_derivative_method='fast',
                                          modes=nms,
                                          mode_transformation=tf,
                                          mode_selection=mode_selection,
                                          logger=log_file
                                          )
        return runner

def run_b3lyp_vpt(
        subkey,
        file_pattern='methanol_vpt_{key}.fchk',
        **args
):
    return run_gaussian_vpt('b3lyp', subkey, file_pattern=file_pattern, **args)

def run_wb97_vpt(
        subkey,
        file_pattern='methanol_vpt_{key}_wb97.fchk',
        **args
):
    return run_gaussian_vpt('wb97', subkey, file_pattern=file_pattern, **args)

def get_mlip_logfile(lot, key, use_reaction_path=True, use_degeneracies=True,
                       mode_selection=None,
                       use_internals=False,
                       step_size=None,
                       analytic_derivative_order=None
                       ):
    # if use_reaction_path:
    #     if mode_selection is None:
    #         ms = "_all"
    #     else:
    #         ms = "_" + ("".join(str(int(s)) for s in mode_selection))
    #     if not use_degeneracies:
    #         deg = "_nodegs"
    #     else:
    #         deg = ""
    #     if use_internals:
    #         log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_rp_ints_{key}{ms}{deg}.out')
    #     else:
    #         log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_{key}{ms}{deg}.out')
    # else:
    #     if mode_selection is None:
    #         ms = "_all"
    #     else:
    #         ms = "_" + ("".join(str(int(s)) for s in mode_selection))
    #     if not use_degeneracies:
    #         deg = "_nodegs"
    #     else:
    #         deg = ""
    #     if use_internals:
    #         log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_norp_ints_{key}{ms}{deg}.out')
    #     else:
    #         log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_norp_{key}{ms}{deg}.out')

    if mode_selection is None:
        ms = "_all"
    else:
        ms = "_" + ("".join(str(int(s)) for s in mode_selection))
    if not use_degeneracies:
        deg = "_nodegs"
    else:
        deg = ""
    if use_internals:
        coord = "ints_"
        rp = "rp_"
    else:
        coord = ""
        rp = ""
    if not use_reaction_path:
        rp = "norp_"
    if step_size is not None:
        step = "_s" + str(round(1000*step_size))
    else:
        step = ""
    if analytic_derivative_order is not None:
        ad = f"_ad{analytic_derivative_order}"
    else:
        ad = ""
    log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_{rp}{coord}{key}{ms}{deg}{step}{ad}.out')
    return log_file

def get_aimnet_logfile(key, use_reaction_path=True, use_degeneracies=True,
                       mode_selection=None,
                       use_internals=False,
                       step_size=None,
                       analytic_derivative_order=None
                       ):
    return get_mlip_logfile('aimnet', key, use_reaction_path=use_reaction_path, use_degeneracies=use_degeneracies,
                       mode_selection=mode_selection,
                       use_internals=use_internals,
                       step_size=step_size,
                       analytic_derivative_order=analytic_derivative_order)

def get_aimnet_old_logfile(key, use_reaction_path=True, use_degeneracies=True,
                       mode_selection=None,
                       use_internals=False,
                       step_size=None,
                       analytic_derivative_order=None
                       ):
    return get_mlip_logfile('aimnet-old', key, use_reaction_path=use_reaction_path, use_degeneracies=use_degeneracies,
                       mode_selection=mode_selection,
                       use_internals=use_internals,
                       step_size=step_size,
                       analytic_derivative_order=analytic_derivative_order)

def get_mace_logfile(key, use_reaction_path=True, use_degeneracies=True,
                       mode_selection=None,
                       use_internals=False,
                       step_size=None,
                       analytic_derivative_order=None
                       ):
    return get_mlip_logfile('mace', key, use_reaction_path=use_reaction_path, use_degeneracies=use_degeneracies,
                       mode_selection=mode_selection,
                       use_internals=use_internals,
                       step_size=step_size,
                       analytic_derivative_order=analytic_derivative_order)


def run_mlip_vpt(lot, logfile_getter, structure_getter, expansion_getter, key,
                 use_internals=False,
                 use_reaction_path=True,
                 log_file=None,
                 overwrite=False,
                 mode_selection=None,
                 recompute_expansion=False,
                 step_size=None,
                 use_degeneracies=True,
                 return_runner=False,
                 analytic_derivative_order=None,
                 **opts
                 ):
    coords = structure_getter(key, 'optimized')
    coords = np.array(coords)
    (modes, status), potential_expansion = expansion_getter(key, 'optimized',
                                                            overwrite=recompute_expansion,
                                                            step_size=step_size,
                                                            analytic_derivative_order=analytic_derivative_order
                                                            )
    modes.matrix = np.array(modes.matrix)
    modes.inverse = np.array(modes.inverse)
    modes.masses = np.array(modes.masses)
    modes.freqs = np.array(modes.freqs)
    modes.origin = coords
    modes.g_matrix = (
        np.array(modes.g_matrix).reshape(18, 18)
            if modes.g_matrix is not None else
        None
    )
    me_aimnet = methanol.modify(
        coords=coords,
        internals=methanol_zmatrix if use_internals else None,
        potential_derivatives=[np.array(v) for v in potential_expansion[1:]]
    )

    if log_file is None:
        log_file = logfile_getter(key,
                                  use_degeneracies=use_degeneracies,
                                  use_reaction_path=(status or use_reaction_path),
                                  mode_selection=mode_selection, use_internals=use_internals,
                                  step_size=step_size,
                                  analytic_derivative_order=analytic_derivative_order
                                  )


    os.makedirs(paths.torsion_scan_path(lot, 'results'), exist_ok=True)
    if status or use_reaction_path:
        locs, nms, rpnms = prep_rp_modes(me_aimnet,
                                         # nms=modes,
                                         internals=methanol_zmatrix,
                                         proj_coord=(2, 0, 1, 5)
                                         )
        tf = locs.localizing_transformation
        if mode_selection is not None:
            tf = locs[mode_selection].localizing_transformation

        print_rpnm_hessian(me_aimnet, nms, rpnms, locs)
        # return nms.make_mass_weighted().apply_transformation(locs.localizing_transformation[1].T)

        mode_selection = None
        # print(tf[0].T @ tf[0] )
        # return
        # tf = np.eye(tf[0].shape[0])[:, 1:]
        # tf = (tf, tf.T)
    else:
        # nms = modes
        _, nms, _ = prep_rp_modes(me_aimnet,
                                         # nms=modes,
                                         internals=methanol_zmatrix,
                                         proj_coord=(2, 0, 1, 5)
                                         )
        tf = None

    # nms = modes
    opts = dict(
        dict(
            states=2,
            degeneracy_specs='auto' if use_degeneracies else None,
            cartesian_analytic_deriv_order=-1,
            cartesian_by_internal_derivative_method='fast',
            gmatrix_tolerance=1 * UnitsData.convert("Wavenumbers", "Hartrees"),
            modes=nms,
            mode_transformation=tf,
            mode_selection=mode_selection,
            logger=log_file
        ),
        **opts
    )
    if overwrite or not os.path.exists(log_file):
        try:
            os.remove(log_file)
        except:
            ...
        runner, _ = me_aimnet.setup_VPT(**opts)
        if return_runner:
            return runner

        try:
            runner.print_tables(logger=runner.hamiltonian.logger, print_intensities=False)
        except:
            os.remove(log_file)
            raise
    elif return_runner:
        runner, _ = me_aimnet.setup_VPT(**opts)
        return runner

def run_aimnet_vpt(key,
                 use_internals=False,
                 use_reaction_path=True,
                 log_file=None,
                 overwrite=False,
                 mode_selection=None,
                 recompute_expansion=False,
                 step_size=None,
                 use_degeneracies=True,
                 return_runner=False,
                 analytic_derivative_order=None,
                 **opts
                 ):
    return run_mlip_vpt(
        'aimnet', get_aimnet_logfile, get_aimnet_structure, get_aimnet_expansion,
        key,
        use_internals=use_internals,
        use_reaction_path=use_reaction_path,
        log_file=log_file,
        overwrite=overwrite,
        mode_selection=mode_selection,
        recompute_expansion=recompute_expansion,
        step_size=step_size,
        use_degeneracies=use_degeneracies,
        return_runner=return_runner,
        analytic_derivative_order=analytic_derivative_order,
        **opts
    )

def run_mace_vpt(key,
                 use_internals=False,
                 use_reaction_path=True,
                 log_file=None,
                 overwrite=False,
                 mode_selection=None,
                 recompute_expansion=False,
                 step_size=None,
                 use_degeneracies=True,
                 return_runner=False,
                 analytic_derivative_order=None,
                 **opts
                 ):
    return run_mlip_vpt(
        'mace', get_mace_logfile, get_mace_structure, get_mace_expansion,
        key,
        use_internals=use_internals,
        use_reaction_path=use_reaction_path,
        log_file=log_file,
        overwrite=overwrite,
        mode_selection=mode_selection,
        recompute_expansion=recompute_expansion,
        step_size=step_size,
        use_degeneracies=use_degeneracies,
        return_runner=return_runner,
        analytic_derivative_order=analytic_derivative_order,
        **opts
    )


def get_log_generator(lot_name:str):
    if isinstance(lot_name, str):
        lot = LevelsOfTheory(lot_name.lower())
    else:
        lot = lot_name
    if lot == LevelsOfTheory.B3LYP:
        return get_b3lyp_logfile
    elif lot == LevelsOfTheory.wB97:
        return get_wb97_logfile
    elif lot == LevelsOfTheory.AIMNet2:
        return get_aimnet_logfile
    elif lot == LevelsOfTheory.AIMNet2Old:
        return get_aimnet_old_logfile
    elif lot == LevelsOfTheory.MACE:
        return get_mace_logfile
    else:
        raise NotImplementedError(lot)


def get_xii_term(runner, i, mode_selection=None, include_coriolis=False):
    if include_coriolis: raise ValueError('dumb')

    freqs = np.diag(runner.hamiltonian.V_terms[0])
    ndim = len(freqs)
    zeta = np.zeros((ndim, ndim, ndim))
    Be = np.zeros((3,))

    v3 = runner.hamiltonian.V_terms[1]
    v4 = runner.hamiltonian.V_terms[2]
    if mode_selection is not None:
        mode_selection = np.arange(ndim)[mode_selection,]
        rem = np.setdiff1d(np.arange(ndim), mode_selection)
        v3 = v3.copy()
        v3[:, :, rem] = 0
        v3[:, rem, :] = 0
        v3[rem, :, :] = 0

        v4 = v4.copy()
        v4[:, :, :, rem] = 0
        v4[:, :, rem, :] = 0
        v4[:, rem, :, :] = 0
        v4[rem, :, :, :] = 0

    if i < 0:
        i = i + ndim

    return PerturbationTheoryHamiltonian._Nielsen_xss(
        i,
        freqs,
        v3,
        v4,
        zeta,
        Be,
        ndim,
        return_components=True
    )


def get_xij_term(runner, i, j, mode_selection=None, include_coriolis=False):
    if include_coriolis: raise ValueError('dumb')

    freqs = np.diag(runner.hamiltonian.V_terms[0])
    ndim = len(freqs)
    zeta = np.zeros((ndim, ndim, ndim))
    Be = np.zeros((3,))

    v3 = runner.hamiltonian.V_terms[1]
    v4 = runner.hamiltonian.V_terms[2]
    if mode_selection is not None:
        mode_selection = np.arange(ndim)[mode_selection,]
        rem = np.setdiff1d(np.arange(ndim), mode_selection)
        v3 = v3.copy()
        v3[:, :, rem] = 0
        v3[:, rem, :] = 0
        v3[rem, :, :] = 0

        v4 = v4.copy()
        v4[:, :, :, rem] = 0
        v4[:, :, rem, :] = 0
        v4[:, rem, :, :] = 0
        v4[rem, :, :, :] = 0

    if i < 0:
        i = i + ndim
    if j < 0:
        j = j + ndim

    return PerturbationTheoryHamiltonian._Nielsen_xst(
        i, j,
        freqs,
        v3,
        v4,
        zeta,
        Be,
        ndim,
        return_components=True
    )


def contrib_x_terms(terms):
    return [
        [y * UnitsData.hartrees_to_wavenumbers for y in x]
        for x in terms
    ]


def sum_x_terms(terms):
    return sum(sum(x) for x in terms) * UnitsData.hartrees_to_wavenumbers