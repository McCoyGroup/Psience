
import numpy as np
from Psience.Molecools import Molecule
import McUtils.Coordinerds as coordops
import McUtils.Scaffolding as scaff
from McUtils.Data import UnitsData

from .paths import *

# methanol = Molecule.from_string(
#     'methanol',
#     energy_evaluator='aimnet2'
# ).optimize(max_iterations=500)
methanol = Molecule(
    ['C', 'O', 'H', 'H', 'H', 'H'],
    [[-0.71174571, 0.0161939, 0.02050266],
     [1.71884591, -1.07310118, -0.2778059],
     [-1.30426891, 0.02589585, 1.99632677],
     [-0.77962613, 1.94036941, -0.7197672],
     [-2.02413643, -1.14525287, -1.05166036],
     [2.91548382, -0.08353621, 0.65084457]],
    energy_evaluator='aimnet2'
)

methanol_mace = Molecule(
    ['C', 'O', 'H', 'H', 'H', 'H'],
    [[-0.7112156961155202, 0.016882573159026193, 0.021090468294894457],
     [1.719641249832828, -1.0741936736456394, -0.2784223473417176],
     [-1.3165664356893934, 0.029604401074940445, 1.9969457089839304],
     [-0.7918278489541754, 1.9444274221797408, -0.7196441670293944],
     [-2.0238536583678375, -1.1455436654754763, -1.0518106466388337],
     [2.9207668881810407, -0.08187322414482251, 0.653037418675572]],
    energy_evaluator='mace'
)

methanol_zmatrix = coordops.reindex_zmatrix(
    coordops.functionalized_zmatrix(
        3,
        {
            (2, 1, 0): [
                [0, -1, -2, -3],
                [1, -1, 0, -2],
                [2, -1, 0, 1],
            ]
        }
    ),
    [5, 1, 0, 2, 3, 4]
)

methanol_gaussian_zmatrix = coordops.functionalized_zmatrix(
    3,
    {
        (2, 1, 0): [
            [0, -1, -2, -3],
            [1, -1, 0, -2],
            [2, -1, 0, 1],
        ]
    }
)


def get_angle_key(key):
    if isinstance(key, (float, np.floating)):
        key = round(np.rad2deg(key))
    return str(int(key))


def get_aimnet_checkpointer(key):
    return scaff.Checkpointer.from_file(
        torsion_scan_path(f'checkpoints/optimized_torsion_potentials_{key}.json')
    )

def get_aimnet_methanol():
    base_chk = scaff.Checkpointer.from_file(
        torsion_scan_path('checkpoints/base_struct.json')
    )

    me_opt = base_chk.cached_eval(
        'opt_mol',
        lambda: methanol.optimize()
    )
    # me_opt = methanol.optimize()
    me_int = me_opt.modify(
        internals=methanol_zmatrix,
    )

    return me_opt, me_int

def get_mace_methanol():
    base_chk = scaff.Checkpointer.from_file(
        torsion_scan_path('checkpoints/base_struct_mace.json')
    )
    meth_mace = methanol.modify(energy_evaluator='mace')

    me_opt = base_chk.cached_eval(
        'opt_mol',
        lambda: meth_mace.optimize()
    )
    # me_opt = methanol.optimize()
    me_int = me_opt.modify(
        internals=methanol_zmatrix,
    )

    return me_opt, me_int

def get_aimnet_structure(key, subkey, coords=None, overwrite=False,
                         chk_pattern='checkpoints/optimized_torsion_potentials_{key}.json'):
    me_opt, me_int = get_aimnet_methanol()
    key = get_angle_key(key)
    chk = scaff.Checkpointer.from_file(
        torsion_scan_path(chk_pattern.format(key=key))
    )
    if coords is None:
        if overwrite:
            del chk[(key, 'unoptimized', 'coords')]
        coords = chk.cached_eval(
            (key, 'unoptimized', 'coords'),
            lambda: me_int.get_displaced_coordinates(
                [[np.deg2rad(float(key))]],
                which=coordops.zmatrix_indices(methanol_zmatrix, [(2, 0, 1, 5)]),
                # which=coordops.zmatrix_indices(methanol_zmatrix, [(4, 0, 2, 3)]),
                use_internals='reembed',
                strip_embedding=True
            )[0]
        )
    if subkey == 'optimized':
        if overwrite:
            del chk[(key, subkey, 'coords')]
        coords = chk.cached_eval(
            (key, subkey, 'coords'),
            lambda: me_int.modify(coords=coords).optimize(
                coordinate_constraints=[(2, 0, 1, 5)],
                tol=1e-8
            ).coords
        )
    return coords
def get_aimnet_expansion(key, subkey, coords=None, overwrite=False, step_size=None, analytic_derivative_order=None,
                         chk_pattern='checkpoints/optimized_torsion_potentials_{key}.json'
                         ):
    coords=get_aimnet_structure(key, subkey, coords=coords, chk_pattern=chk_pattern)
    key=get_angle_key(key)
    chk = scaff.Checkpointer.from_file(
        torsion_scan_path(chk_pattern.format(key=key))
    )
    with chk:
        mol = methanol.modify(coords=coords)
        # me_aimnet2 = methanol.get_energy_function('aimnet2')
        harmonic_expansion=chk.cached_eval(
            (key, subkey, 'harmonic_expansion'),
            lambda: mol.calculate_energy(order=2)
        )
        modes, status = chk.cached_eval(
            (key, subkey, 'modes'),
            lambda: mol.get_reaction_path_modes(
                potential_derivatives=harmonic_expansion[1:],
                return_status=True
            )
        )
        if step_size is None:
            if analytic_derivative_order is None:
                tag = (key, subkey, 'potential_derivatives')
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes
                    )
                )
            else:
                tag = (key, subkey, 'potential_derivatives_analytic_do', str(analytic_derivative_order))
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes,
                        analytic_derivative_order=analytic_derivative_order
                    )
                )
        else:
            tag = (key, subkey, 'potential_derivatives_stepped', str(round(step_size * 1000)))
            if analytic_derivative_order is None:
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes,
                        mesh_spacing=step_size
                    )
                )
            else:
                tag = (
                            key, subkey,
                            'potential_derivatives_stepped_analytic',
                            str(round(step_size * 1000)),
                            str(analytic_derivative_order)
                        )
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes,
                        mesh_spacing=step_size,
                        analytic_derivative_order=analytic_derivative_order
                    )
                )
    return (modes, status), potential_derivatives

def get_mace_structure(key, subkey, coords=None, overwrite=False):
    me_opt, me_int = get_mace_methanol()
    key = get_angle_key(key)
    chk = scaff.Checkpointer.from_file(
        torsion_scan_path(f'checkpoints/optimized_torsion_potentials_mace_{key}.json')
    )
    if coords is None:
        if overwrite:
            try:
                del chk[(key, 'unoptimized', 'coords')]
            except KeyError:
                pass
        coords = chk.cached_eval(
            (key, 'unoptimized', 'coords'),
            lambda: me_int.get_displaced_coordinates(
                [[np.deg2rad(float(key))]],
                which=coordops.zmatrix_indices(methanol_zmatrix, [(2, 0, 1, 5)]),
                # which=coordops.zmatrix_indices(methanol_zmatrix, [(4, 0, 2, 3)]),
                use_internals='reembed',
                strip_embedding=True
            )[0]
        )
    if subkey == 'optimized':
        if overwrite:
            try:
                del chk[(key, subkey, 'coords')]
            except KeyError:
                pass
        coords = chk.cached_eval(
            (key, subkey, 'coords'),
            lambda: me_int.modify(coords=coords).optimize(
                coordinate_constraints=[(2, 0, 1, 5)],
                max_iterations=150,
                method='scipy',
                logger=True
            ).coords
        )
    return coords

def get_mace_expansion(key, subkey, coords=None, overwrite=False, step_size=None, analytic_derivative_order=None,
                       return_harmonic=False
                       ):
    coords=get_mace_structure(key, subkey, coords=coords)
    key=get_angle_key(key)
    chk = scaff.Checkpointer.from_file(
        torsion_scan_path(f'checkpoints/optimized_torsion_potentials_mace_{key}.json')
    )
    with chk:
        mol = methanol_mace.modify(coords=coords)
        if overwrite:
            try:
                del chk[(key, subkey, 'harmonic_expansion')]
            except KeyError:
                pass
        harmonic_expansion=chk.cached_eval(
            (key, subkey, 'harmonic_expansion'),
            lambda: mol.calculate_energy(order=2)
        )
        if return_harmonic: return harmonic_expansion
        modes, status = chk.cached_eval(
            (key, subkey, 'modes'),
            lambda: mol.get_reaction_path_modes(
                potential_derivatives=harmonic_expansion[1:],
                return_status=True,
                zero_gradient_cutoff=25 / UnitsData.hartrees_to_wavenumbers
            )
        )
        if step_size is None:
            if analytic_derivative_order is None:
                tag = (key, subkey, 'potential_derivatives')
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes
                    )
                )
            else:
                tag = (key, subkey, 'potential_derivatives_analytic_do', str(analytic_derivative_order))
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes,
                        analytic_derivative_order=analytic_derivative_order
                    )
                )
        else:
            tag = (key, subkey, 'potential_derivatives_stepped', str(round(step_size * 1000)))
            if analytic_derivative_order is None:
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes,
                        mesh_spacing=step_size
                    )
                )
            else:
                tag = (
                            key, subkey,
                            'potential_derivatives_stepped_analytic',
                            str(round(step_size * 1000)),
                            str(analytic_derivative_order)
                        )
                if overwrite:
                    try:
                        del chk[tag]
                    except KeyError:
                        pass
                potential_derivatives = chk.cached_eval(
                    tag,
                    lambda: mol.partial_force_field(
                        order=4,
                        modes=modes,
                        mesh_spacing=step_size,
                        analytic_derivative_order=analytic_derivative_order
                    )
                )
    return (modes, status), potential_derivatives

cpmo3_zmatrix = coordops.reindex_zmatrix(
    coordops.functionalized_zmatrix(
        1,
        {
            (0, -1, -2): coordops.chain_zmatrix(2),
            (0, 1, 2): coordops.chain_zmatrix(2),
            (0, 3, 4): coordops.chain_zmatrix(2),
            (0, 5, 6): coordops.chain_zmatrix(2),
            (0, 7, 8): coordops.chain_zmatrix(2),
            (0, 9, 10): coordops.chain_zmatrix(2),
            (0, 11, 12): coordops.chain_zmatrix(2),
            (0, 13, 14): coordops.chain_zmatrix(2)
        }
    ),
    [10, 11, 12, 13, 14, 15, 16, 0, 9, 1, 8, 2, 7, 3, 6, 4, 5]
)

def load_cpmo3(struct_id):
    return Molecule.from_file(
        path('cpmo3_dist_scan', f'cpmo3m_opt_dist_{struct_id}.fchk')
    )
