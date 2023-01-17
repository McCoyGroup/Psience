import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.expanduser("~/Documents/UW/Research/Development"))
from Psience.VPT2 import *

os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

def conv(r, t, f, **kwargs):
    cp1 = np.cos(f[..., 3])  # skip three for embedding
    ct1 = np.cos(t[..., 2])  # skip two for embedding
    ct2 = np.cos(t[..., 3])
    st1 = np.sin(t[..., 2])
    st2 = np.sin(t[..., 3])
    f[..., 3] = np.arccos(st1 * st2 * cp1 + ct1 * ct2)
    return np.array([r, t, f])

def inv(r, t, f, **kwargs):
    cp1 = np.cos(f[..., 3])
    ct1 = np.cos(t[..., 2])
    ct2 = np.cos(t[..., 3])
    st1 = np.sin(t[..., 2])
    st2 = np.sin(t[..., 3])
    f[..., 3] = np.arccos((cp1 - ct1 * ct2) / (st1 * st2))
    return np.array([r, t, f])

VPTRunner.helpers.run_anne_job(
    '.',
    2,
    order=2,
    # logger='nh3_r.out',
    degeneracy_specs='auto',
    # coordinate_transformation=[conv, inv],
    check_overlap=False,
    results_file='nh3.hdf5'
    # degeneracy_specs={
    #     'nT':[1, 1, 1, 1, 2, 2]
    # }
)