import os, numpy as np
from McUtils.Extensions import DynamicFFILibrary, CLoader # we'll load this at runtime

__all__ = [
    "MBPol"
]

cur_dir = os.path.dirname(os.path.abspath(__file__))
lib_file = os.path.join(cur_dir, "libs", "libmbpol.so")
if not os.path.isfile(lib_file):
    CLoader(cur_dir).custom_make(
        True, os.path.join(cur_dir, "libs", 'mbpol')
    )

# lib_file = os.path.join(test_dir, "libmbpol.so")
MBPol = DynamicFFILibrary(
    lib_file,
    get_pot=dict(
        name='calcpot_',
        nwaters=(int,),
        energy=[float],
        coords=[float],
        return_type=float,
        prep_args=lambda kw:[
            kw.__setitem__('energy', np.zeros(len(kw['coords']) if kw['coords'].ndim > 2 else 1)),
            kw][-1],
        defaults={'energy': None},
        return_handler=lambda r, kw: kw['energy'] / 627.5094740631
    ),
    get_pot_grad = dict(
        name='calcpotg_',
        nwaters=(int,),
        energy=(float,),
        coords=[float],
        grad=[float],
        return_type=float,
        prep_args=lambda kw:[
            kw.__setitem__('grad', np.zeros(kw['coords'].shape)),
            kw.__setitem__('energy', np.zeros(len(kw['coords']) if kw['coords'].ndim > 2 else 1)),
            kw][-1],
        defaults={'grad':None, 'energy': None},
        return_handler=lambda r, kw: {
            'grad':kw['grad'] / 627.5094740631,
            'energy': kw['energy'] / 627.5094740631
        }
    ),
    compiler_options=dict(
        threaded=True
    )
)