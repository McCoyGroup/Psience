import enum
import numpy as np


__all__ = [
    "RotorTypes",
    "identify_rotor_type"
]

class RotorTypes(enum.Enum):
    Atom = "atom"
    Linear = "linear"
    Planar = "planar"
    Oblate = "oblate"
    Prolate = "prolate"
    Spherical = "spherical"
    Asymmetric = "asymmetric"

def identify_rotor_type(moms:np.ndarray, tol=1e-8):
    mom_deg, counts = np.unique(tol*np.round(moms/tol), return_counts=True)
    planar = np.abs((moms[0] + moms[1]) - moms[2]) < tol
    if np.min(mom_deg) < 1e-6:
        if len(np.where(moms < 1e-6)[0]) == 2:
            type = RotorTypes.Atom
        else:
            type = RotorTypes.Linear
    elif tuple(counts) == (3,):
        type = RotorTypes.Spherical
    elif tuple(counts) == (1, 2):
        type = RotorTypes.Prolate
    elif tuple(counts) == (2, 1):
        type = RotorTypes.Oblate
    else:
        type = RotorTypes.Asymmetric

    return type, planar