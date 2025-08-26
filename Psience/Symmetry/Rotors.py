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

def identify_rotor_type(moms:np.ndarray, tol=1e-8, zero_tol=1e-6):
    diffs = np.diff(moms)
    eq_pos = np.where(diffs < tol)
    if len(eq_pos) > 0:
        pos = eq_pos[0]
    else:
        pos = []
    planar = np.abs((moms[0] + moms[1]) - moms[2]) < tol
    if np.min(moms) < zero_tol:
        if len(np.where(moms < zero_tol)[0]) == 2:
            type = RotorTypes.Atom
        else:
            type = RotorTypes.Linear
    elif len(pos) == 2:
        type = RotorTypes.Spherical
    elif len(pos) == 1:
        if pos[0] == 0:
            type = RotorTypes.Oblate
        else:
            type = RotorTypes.Prolate
    else:
        type = RotorTypes.Asymmetric

    return type, planar