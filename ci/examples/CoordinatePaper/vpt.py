import os

import numpy as np
from Psience.Molecools import Molecule

from . import paths
from .expansions import (
    get_aimnet_structure, get_aimnet_expansion,
    methanol, methanol_zmatrix, methanol_gaussian_zmatrix
)
from .modes import prep_rp_modes, print_rpnm_hessian

def get_gaussian_logfile(lot, subkey, use_reaction_path=True, mode_selection=None, use_internals=False):
    ms = "".join(str(m) for m in
                 (mode_selection if mode_selection is not None else 'all')
                 )

    if use_reaction_path:
        if use_internals:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_rp_ints_{subkey}_{ms}.out')
        else:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_rp_{subkey}_{ms}.out')
    else:
        if use_internals:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_norp_ints_{subkey}_{ms}.out')
        else:
            log_file = paths.torsion_scan_path(lot, 'results', f'methanol_vpt_norp_{subkey}_{ms}.out')

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
    log_file = get_gaussian_logfile(lot, subkey,
                                    use_reaction_path=(stat or use_reaction_path),
                                    mode_selection=mode_selection,
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
                                          degeneracy_specs='auto',
                                          cartesian_analytic_deriv_order=-1,
                                          internal_by_cartesian_derivative_method='fast',
                                          modes=nms,
                                          mode_transformation=tf,
                                          mode_selection=mode_selection,
                                          logger=log_file
                                          )
        runner.print_tables(logger=runner.hamiltonian.logger)

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

def get_aimnet_logfile(key, use_reaction_path=True, mode_selection=None, use_internals=False):
    if use_reaction_path:
        if mode_selection is None:
            ms = "_all"
        else:
            ms = "_" + ("".join(str(int(s)) for s in mode_selection))
        if use_internals:
            log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_rp_ints_{key}{ms}.out')
        else:
            log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_{key}{ms}.out')
    else:
        if mode_selection is None:
            ms = ""
        else:
            ms = "_" + ("".join(str(int(s)) for s in mode_selection))
        if use_internals:
            log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_norp_ints_{key}{ms}.out')
        else:
            log_file = paths.torsion_scan_path('aimnet', 'results', f'methanol_vpt_norp_{key}{ms}.out')
    return log_file

def run_aimnet_vpt(key,
                   use_internals=False,
                   use_reaction_path=True,
                   log_file=None,
                   overwrite=False,
                   mode_selection=None,
                   recompute_expansion=False,
                   **opts
                   ):
    coords = get_aimnet_structure(key, 'optimized')
    coords = np.array(coords)
    (modes, status), potential_expansion = get_aimnet_expansion(key, 'optimized', overwrite=recompute_expansion)
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

    os.makedirs(paths.torsion_scan_path('aimnet', 'results'), exist_ok=True)
    if status or use_reaction_path:
        locs, nms, rpnms = prep_rp_modes(me_aimnet,
                                         nms=modes,
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
        tf = None

    log_file = get_aimnet_logfile(key, use_reaction_path=status, mode_selection=mode_selection, use_internals=use_internals)

    if overwrite or not os.path.exists(log_file):
        try:
            os.remove(log_file)
        except:
            ...
        nms = modes
        opts = dict(
            dict(
                states=2,
                degeneracy_specs='auto',
                cartesian_analytic_deriv_order=-1,
                internal_by_cartesian_derivative_method='fast',
                modes=nms,
                mode_transformation=tf,
                mode_selection=mode_selection,
                logger=log_file
            ),
            **opts
        )

        runner, _ = me_aimnet.setup_VPT(**opts)
        runner.print_tables(logger=runner.hamiltonian.logger, print_intensities=False)