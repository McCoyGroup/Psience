"""
Provides a class for handling a compiled set of atomic data
"""
import os, sys
from McUtils.Data import DataHandler
__all__ = ['KEData', 'KEDataHandler']

class KEDataHandler(DataHandler):
    """
    A DataHandler that's built for use with the G-matrix and V' terms
    from the 1999 Frederick and Woywood paper
    Usually used through the `KETermData` object.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))

    def __init__(self):
        ...
    equivalent_perms = {'r': [[1, 0]], 'a': [[2, 1, 0]], 't': [[3, 2, 1, 0]], 'y': []}

    def _get_remapping(self, i1, i2):
        ...

    def find_expressions(self, k, return_permutation=False):
        ...
KEData = KEDataHandler()
KEData.__doc__ = "An instance of KEDataHandler that can be used for looking up G-matrix and V' data"
KEData.__name__ = 'KEData'