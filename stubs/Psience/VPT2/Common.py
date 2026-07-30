import numpy as np
__all__ = ['PerturbationTheoryException', 'Settings']

class PerturbationTheoryException(Exception):
    ...

class Settings:
    non_zero_cutoff = 1e-14

def _safe_dot(a, b):
    ...