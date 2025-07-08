
import numpy as np
from Psience.Molecools import Molecule

from .expansions import (
    get_aimnet_structure, get_aimnet_expansion,
    methanol, methanol_zmatrix, methanol_gaussian_zmatrix
)
from .modes import prep_rp_modes, print_rpnm_hessian


def run_gaussian_vpt(target_dir,
                     which=None,
                     file_pattern='methanol_vpt_{i}.fchk',
                     use_internals=False,
                     use_reaction_path=True
                     ):
    if which is None:
        which = list(range(which))

    for i in which:
        print("=" * 100)

        fp = file_pattern.format(i=i + 1)

        me_gaussian = Molecule.from_file(
            f'{target_dir}/{fp}',
            internals=methanol_gaussian_zmatrix if use_internals else None
        )

        (locs, nms, rpnms), stat = prep_rp_modes(me_gaussian,
                                                 internals=methanol_gaussian_zmatrix,
                                                 proj_coord=(3, 2, 1, 0),
                                                 return_status=True
                                                 )

        if stat or use_reaction_path:
            tf = locs.localizing_transformation

            print_rpnm_hessian(me_gaussian, nms, rpnms, locs)
            if use_internals:
                log_file = f'{target_dir}/methanol_vpt_rp_ints_{i + 1}.out'
            else:
                log_file = f'{target_dir}/methanol_vpt_{i + 1}.out'
        else:
            tf = None

            if use_internals:
                log_file = f'{target_dir}/methanol_vpt_norp_ints_{i + 1}.out'
            else:
                log_file = f'{target_dir}/methanol_vpt_norp_{i + 1}.out'

        runner, _ = me_gaussian.setup_VPT(states=2,
                                          degeneracy_specs='auto',
                                          cartesian_analytic_deriv_order=-1,
                                          internal_by_cartesian_derivative_method='fast',
                                          modes=nms,
                                          mode_transformation=tf,
                                          logger=log_file
                                          )
        runner.print_tables(logger=runner.hamiltonian.logger)

def run_aimnet_vpt(key, target_dir="torsion_scan/aimnet",
                   use_internals=False,
                   use_reaction_path=True,
                   log_file=None,
                   overwrite=False,
                   mode_selection=None,
                   **opts
                   ):
    coords = get_aimnet_structure(key, 'optimized')
    coords = np.array(coords)
    (modes, status), potential_expansion = get_aimnet_expansion(key, 'optimized', overwrite=overwrite)
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
    if status or use_reaction_path:
        locs, nms, rpnms = prep_rp_modes(me_aimnet,
                                         nms=modes,
                                         internals=methanol_zmatrix,
                                         proj_coord=(2, 0, 1, 5)
                                         )
        tf = locs.localizing_transformation

        print_rpnm_hessian(me_aimnet, nms, rpnms, locs)
        # return nms.make_mass_weighted().apply_transformation(locs.localizing_transformation[1].T)

        if log_file is None:
            if mode_selection is None:
                ms = ""
            else:
                ms = "_ " +("-".join(str(int(s)) for s in mode_selection))
            if use_internals:
                log_file = f'{target_dir}/methanol_vpt_rp_ints{ms}_{key}.out'
            else:
                log_file = f'{target_dir}/methanol_vpt{ms}_{key}.out'

        # print(tf[0].T @ tf[0] )
        # return
        # tf = np.eye(tf[0].shape[0])[:, 1:]
        # tf = (tf, tf.T)
    else:
        tf = None


        if log_file is None:
            if mode_selection is None:
                ms = ""
            else:
                ms = "_ " +("-".join(str(int(s)) for s in mode_selection))
            if use_internals:
                log_file = f'{target_dir}/methanol_vpt_norp_ints{ms}_{key}{ms}.out'
            else:
                log_file = f'{target_dir}/methanol_vpt_norp{ms}_{key}.out'


    nms = modes

    opts = dict(
        dict(
            states=2,
            degeneracy_specs='auto',
            cartesian_analytic_deriv_order=-1,
            internal_by_cartesian_derivative_method='fast',
            modes=nms,
            mode_transformation=tf,
            logger=log_file,
            mode_selection=mode_selection
        ),
        **opts
    )

    runner, _ = me_aimnet.setup_VPT(**opts)
    runner.print_tables(logger=runner.hamiltonian.logger, print_intensities=False)