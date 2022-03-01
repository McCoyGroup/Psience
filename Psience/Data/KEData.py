"""
Provides a class for handling a compiled set of atomic data
"""

import os, sys
from McUtils.Data import DataHandler

__all__ = [ "KEData", "KEDataHandler" ]

class KEDataHandler(DataHandler):
    """
    A DataHandler that's built for use with the G-matrix and V' terms
    from the 1999 Frederick and Woywood paper
    Usually used through the `KETermData` object.
    """

    base_dir = os.path.dirname(os.path.abspath(__file__))
    def __init__(self):
        super().__init__("KEData", data_dir=self.base_dir, data_pkg='PsiDatasets')#, alternate_keys=("Name", "Symbol", "CanonicalSymbol"))
KEData = KEDataHandler()
KEData.__doc__ = """An instance of KEDataHandler that can be used for looking up G-matrix and V' data"""
KEData.__name__ = "KEData"