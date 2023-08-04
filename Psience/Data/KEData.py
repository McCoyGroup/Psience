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
    equivalent_perms = {
        "r":[[1, 0]],
        "a":[[2, 1, 0]],
        "t":[[3, 2, 1, 0]],
        "y":[]
    }
    def _get_remapping(self, i1, i2):
        remapping = {}
        n = 1
        _ = []
        for i in i1:
            if i not in remapping:
                remapping[i] = n
                _.append(n)
                n += 1
            else:
                _.append(remapping[i])
        i1 = tuple(_)

        _ = []
        for i in i2:
            if i not in remapping:
                remapping[i] = n
                _.append(n)
                n += 1
            else:
                _.append(remapping[i])
        i2 = tuple(_)

        return i1, i2
    def find_expressions(self, k, return_permutation=False):
        try:
            val = self[(k,)]
        except KeyError:
            pass
        else:
            if return_permutation:
                return val, None
            else:
                return val

        try:
            ((type1, type2), inds1, inds2) = k
        except:
            ((type1, type2), (inds1, inds2)) = k

        i1, i2 = self._get_remapping(inds1, inds2)
        k = ((type1, type2), i1, i2)
        try:
            val = self[(k,)]
        except KeyError:
            pass
        else:
            if return_permutation:
                return val, None
            else:
                return val


        for p1 in self.equivalent_perms[type1]:
            i1, i2 = self._get_remapping(tuple(inds1[j] for j in p1), inds2)
            k = ((type1, type2), i1, i2)
            try:
                val = self[(k,)]
            except KeyError:
                pass
            else:
                if return_permutation:
                    return val, (p1, None)
                else:
                    return val

            for p2 in self.equivalent_perms[type2]:
                i1, i2 = self._get_remapping(inds1, tuple(inds2[j] for j in p2))
                k = ((type1, type2), i1, i2)
                try:
                    val = self[(k,)]
                except KeyError:
                    pass
                else:
                    if return_permutation:
                        return val, (None, p2)
                    else:
                        return val

                i1, i2 = self._get_remapping(tuple(inds1[j] for j in p1), tuple(inds2[j] for j in p2))
                k = ((type1, type2), i1, i2)
                try:
                    val = self[(k,)]
                except KeyError:
                    pass
                else:
                    if return_permutation:
                        return val, (p1, p2)
                    else:
                        return val

        if return_permutation:
            return [0, 0], None
        else:
            return [0, 0]

    # def load(self):
    #     # now update by max equivalent indices
    #     super().load()
    #     new_vals = {}
    #
    #     self._data.update(new_vals)
KEData = KEDataHandler()
KEData.__doc__ = """An instance of KEDataHandler that can be used for looking up G-matrix and V' data"""
KEData.__name__ = "KEData"