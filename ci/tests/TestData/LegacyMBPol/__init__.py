import os, numpy as np
from McUtils.Extensions import DynamicFFILibrary, CLoader # we'll load this at runtime

__all__ = [
    "MBPol",
    "potential"
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

def potential(coords, deriv_order=0, nwaters=1, chunk_size=int(5e5)):
    from McUtils.Data import UnitsData
    from McUtils.Zachary import FiniteDifferenceDerivative

    nx = 9*nwaters
    b2a = UnitsData.convert("BohrRadius", "Angstroms")
    coords = coords.reshape(-1, nx)

    just_vals = deriv_order is None
    if just_vals:
        deriv_order = 0

    chunks = [[] for _ in range(deriv_order + 1)]
    num_chunks = int(len(coords) / chunk_size) + 1

    for coords in np.array_split(coords, num_chunks):
        # interp.logger.log_print("evaluating energies")
        energies = MBPol.get_pot(coords=coords.reshape(-1, 3*nwaters, 3) * b2a, nwaters=nwaters,
                                 threading_vars=['energy', 'coords'], threading_mode='omp')
        if deriv_order > 0:
            derivs = []
            grads = lambda c: MBPol.get_pot_grad(
                coords=c.reshape(-1, 3*nwaters, 3) * b2a,
                nwaters=nwaters,
                threading_vars=['energy', 'grad', 'coords'],
                threading_mode='omp'
            )['grad'].reshape(c.shape) * b2a
            # interp.logger.log_print("evaluating forces")
            derivs.append(grads(coords))
            if deriv_order > 1:
                # interp.logger.log_print("evaluating Hessians")
                hess_fun = FiniteDifferenceDerivative(
                    grads,
                    function_shape=(nx, nx)
                )

                # chunks = []
                # num_chunks = int(len(coords)/1000)
                # for a in np.array_split(coords, num_chunks):
                #     chunks.append(hess_fun.derivatives(a).compute_derivatives(1))
                # hess = np.concatenate(chunks, axis=0)
                new_derivs = hess_fun.derivatives(coords).derivative_tensor(list(range(1, deriv_order)))
                # print([d.shape for d in new_derivs])
                derivs.extend(
                    np.moveaxis(d, -2, 0)
                    for i, d in enumerate(new_derivs)
                )
            # interp.logger.log_print("done")
            for i, d in enumerate([energies] + derivs):
                chunks[i].append(d)
        else:
            chunks[0].append(energies)

    for i, c in enumerate(chunks):
        chunks[i] = np.concatenate(c, axis=0)

    if just_vals:
        chunks = chunks[0]

    return chunks
