import numpy as np

def symmetrize_potential(scan_spec, potential):
    # TODO: make symmetrize variable number of times

    symm_pot = np.concatenate([
        potential,
        np.flip(potential[:-1]),
        potential[1:]
    ])

    scan_spec = np.asanyarray(scan_spec)
    scan_min_diffs = scan_spec - np.min(scan_spec)
    scan_max_diffs = np.flip(np.max(scan_spec) - scan_spec)
    scan_max = np.max(scan_spec)
    scan_width = scan_max - np.min(scan_spec)

    symm_scan = np.concatenate([
        scan_spec,
        scan_max + scan_max_diffs[1:],
        scan_max + scan_width + scan_min_diffs[1:]
    ])

    return symm_scan, symm_pot