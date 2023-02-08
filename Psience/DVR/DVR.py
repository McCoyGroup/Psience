
import numpy as np
from McUtils.Scaffolding import ParameterManager

from .ColbertMiller import PolarDVR, RingDVR, CartesianDVR
from .DirectProduct import DirectProductDVR
from .Extensions import SelfConsistentDVR, PotentialOptimizedDVR

__all__ = [
    "DVR"
]

class DVRConstructor:

    _domain_map = None
    @classmethod
    def load_domain_map(cls):

        return {
            (0, np.pi): PolarDVR,
            (0, 2*np.pi): RingDVR,
            None: CartesianDVR
        }
    @classmethod
    def infer_DVR_type(cls, domain):
        if cls._domain_map is None:
            cls._domain_map = cls.load_domain_map()
        for k,v in cls._domain_map.items():
            if k is not None:
                if np.allclose(k, domain):
                    return v
        else:
            return cls._domain_map[None]

    @classmethod
    def construct(cls,
                  domain=None,
                  divs=None,
                  potential_function=None,
                  g=None,
                  g_deriv=None,
                  mass=None,
                  po_divs=25,
                  classes=None,
                  scf=False,
                  potential_optimize=False,
                  logger=True,
                  **base_opts
                  ):

        # dispatches based on domain to construct the appropriate DVR
        if domain is None or divs is None:
            raise ValueError("can't have `None` for `domain` or `divs`")
        if isinstance(domain[0], (int, float, np.integer, np.floating)): # 1D
            domain = [domain]
            divs = [divs]
            if mass is not None:
                mass = [mass]
        if classes is None:
            classes = [None] * len(domain)
        if g is not None:
            if callable(g):
                subg = [g]
            else:
                subg = [g[i][i] for i in range(len(g))]
            mass = [None] * len(subg)
            if g_deriv is None:
                g_deriv = [None] * len(subg)
        else:
            subg = [None]*len(mass)
            g_deriv = [None]*len(mass)

        if isinstance(po_divs, int):
            po_divs = [po_divs]*len(mass)

        ndim = len(list(zip(domain, divs, classes, mass, subg, g_deriv)))
        if ndim == 1:
            classes = [
                cls.infer_DVR_type(r) if c is None else c
                for r, n, c, m, sg, gd, nwf in zip(domain, divs, classes, mass, subg, g_deriv, po_divs)
            ]

            dvr = classes[0](
                domain=domain[0],
                divs=divs[0],
                potential_function=potential_function,
                g=subg[0],
                mass=mass[0],
                g_deriv=g_deriv[0],
                logger=logger,
                **base_opts
            )
        else:
            dvrs_1D = [
                cls.infer_DVR_type(r)(domain=r, divs=n, mass=m, g=sg, g_deriv=gd, num_wfns=nwf) if c is None else c(domain=r, divs=n)
                for r, n, c, m, sg, gd, nwf in zip(domain, divs, classes, mass, subg, g_deriv, po_divs)
            ]
            dvr = DirectProductDVR(
                dvrs_1D,
                domain=domain,
                divs=divs,
                potential_function=potential_function,
                g=g,
                mass=mass,
                g_deriv=g_deriv,
                logger=logger if not potential_optimize or scf else None,
                **ParameterManager(base_opts).exclude((SelfConsistentDVR, PotentialOptimizedDVR))
            )
            
            if potential_optimize or scf:
                if potential_optimize and scf is False:
                    dvr = PotentialOptimizedDVR.from_minimum(dvr,
                                                         logger=logger,
                                                         **ParameterManager(base_opts).filter(PotentialOptimizedDVR)
                                                         )
                else:
                    dvr = SelfConsistentDVR(dvr, logger=logger if not potential_optimize else None, **ParameterManager(base_opts).filter(SelfConsistentDVR))
                    if potential_optimize:
                        dvr = PotentialOptimizedDVR.from_scf(dvr,
                                                             logger=logger,
                                                             **ParameterManager(base_opts).filter(PotentialOptimizedDVR)
                                                             )
        return dvr


def DVR(
        domain=None,
        divs=None,
        classes=None,
        potential_function=None,
        g=None,
        g_deriv=None,
        scf=False,
        potential_optimize=False,
        **base_opts
):
    """
    Constructs a DVR object

    :param domain:
    :type domain:
    :param divs:
    :type divs:
    :param classes:
    :type classes:
    :param potential_function:
    :type potential_function:
    :param g:
    :type g:
    :param g_deriv:
    :type g_deriv:
    :param base_opts:
    :type base_opts:
    :return:
    :rtype:
    """

    return DVRConstructor.construct(
        domain=domain,
        divs=divs,
        classes=classes,
        potential_function=potential_function,
        g=g,
        g_deriv=g_deriv,
        scf=scf,
        potential_optimize=potential_optimize,
        **base_opts
    )