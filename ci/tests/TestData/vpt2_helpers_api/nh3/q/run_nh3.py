import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.expanduser("~/Documents/UW/Research/Development"))
from Psience.VPT2 import *

os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

VPTRunner.helpers.run_anne_job(
    '.',
    2,
    order=2,
    logger='nh3_r.out',
    degeneracy_specs='auto',
    # coordinate_transformation=[conv, inv]
    # degeneracy_specs={
    #     'nT':[1, 1, 1, 1, 2, 2]
    # }
)