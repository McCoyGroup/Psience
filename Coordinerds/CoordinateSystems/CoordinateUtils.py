"""
Little utils that both CoordinateSet and CoordinateSystem needed
"""

import numpy as np

def is_multiconfig(coords):
    return len(coords.shape) > 2

def mc_safe_apply(self, fun, coords = None):
    """Applies fun to the coords in such a way that it will apply to an array of valid
    coordinates (as determined by dimension of the basis). This will either be a single configuration
    or multiple configurations

    :param fun:
    :type fun:
    :return:
    :rtype:
    """
    if coords is None:
        coords = self
    if is_multiconfig(coords):
        base_shape = coords.shape
        new_shape = (np.product(base_shape[:-2]),) + base_shape[-2:]
        coords = np.reshape(coords, new_shape)
        new_coords = fun(coords)
        revert_shape = tuple(base_shape[:-2]) + new_coords.shape[1:]
        new_coords = np.reshape(new_coords, revert_shape)
    else:
        new_coords = fun(coords)
    return new_coords