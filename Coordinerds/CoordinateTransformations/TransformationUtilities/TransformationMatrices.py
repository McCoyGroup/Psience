
from .VectorOps import vec_normalize
import math, numpy as np

#######################################################################################################################
#
#                                                 rotation_matrix
#

def rotation_matrix_basic(xyz, theta):
    """rotation matrix about x, y, or z axis

    :param xyz: x, y, or z axis
    :type xyz: str
    :param theta: counter clockwise angle in radians
    :type theta: float
    """

    axis = xyz.lower()
    if axis == "z": # most common case so it comes first
        mat = [
            [ math.cos(theta), math.sin(theta), 0.],
            [-math.sin(theta), math.cos(theta), 0.],
            [0.,               0.,              1.]
        ]
    elif axis == "y":
        mat = [
            [ math.cos(theta), 0., math.sin(theta)],
            [0.,               1.,              0.],
            [-math.sin(theta), 0., math.cos(theta)]
        ]
    elif axis == "x":
        mat = [
            [1.,               0.,              0.],
            [0.,  math.cos(theta), math.sin(theta)],
            [0., -math.sin(theta), math.cos(theta)]
        ]
    else:
        raise Exception("{}: axis '{}' invalid".format('rotation_matrix_basic', xyz))
    return np.array(mat)

def rotation_matrix_basic_vec(xyz, thetas):
    """rotation matrix about x, y, or z axis

    :param xyz: x, y, or z axis
    :type xyz: str
    :param thetas: counter clockwise angle in radians
    :type thetas: float
    """

    thetas = np.asarray(thetas)
    nmats = len(thetas)
    z = np.zeros((nmats,))
    o = np.ones((nmats,))
    c = np.cos(thetas)
    s = np.sin(thetas)
    axis = xyz.lower()
    if axis == "z": # most common case so it comes first
        mat = [
            [ c, s, z],
            [-s, c, z],
            [ z, z, o]
        ]
    elif axis == "y":
        mat = [
            [ c, z, s],
            [ z, o, z],
            [-s, z, c]
        ]
    elif axis == "x":
        mat = [
            [ o, z, z],
            [ z, c, s],
            [ z,-s, c]
        ]
    else:
        raise Exception("{}: axis '{}' invalid".format('rotation_matrix_basic', xyz))
    return np.array(mat).T

#thank you SE for the nice Euler-Rodrigues imp: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotation_matrix_ER(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad),     2 * (bd - ac)    ],
        [2 * (bc - ad),     aa + cc - bb - dd, 2 * (cd + ab)    ],
        [2 * (bd + ac),     2 * (cd - ab),     aa + dd - bb - cc]
    ])

def rotation_matrix_ER_vec(axes, thetas):
    """Vectorized version of baisc ER
    """

    axes = np.asarray(axes)
    thetas = np.asarray(thetas)
    if len(axes.shape) == 1:
        axes = axes/np.linalg.norm(axes)
        axes = np.broadcast_to(axes, (len(thetas), 3))
    else:
        axes = vec_normalize(axes)

    a = np.cos(thetas/2.0)
    b, c, d = ( -axes * np.reshape(np.sin(thetas / 2.0), (len(thetas), 1)) ).T
    # raise Exception(axes)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad),     2 * (bd - ac)    ],
        [2 * (bc - ad),     aa + cc - bb - dd, 2 * (cd + ab)    ],
        [2 * (bd + ac),    2 * (cd - ab),     aa + dd - bb - cc]
    ]).T

def rotation_matrix(axis, theta):
    try:
        flen = len(theta)
    except TypeError:
        flen = 0
    if type(axis) == str:
        if flen >0:
            mat_fun = rotation_matrix_basic_vec
        else:
            mat_fun = rotation_matrix_basic
    else:
        if flen > 0:
            mat_fun = rotation_matrix_ER_vec
        else:
            mat_fun = rotation_matrix_ER

    return mat_fun(axis, theta)

#######################################################################################################################
#
#                                                 translation_matrix
#

def translation_matrix(shift):
    share = np.asarray(shift)
    if len(share.shape) == 1:
        ss = share
        zs = 0.
        os = 1.
        mat = np.array(
            [
                [os, zs, zs, ss[0]],
                [zs, os, zs, ss[1]],
                [zs, zs, os, ss[2]],
                [zs, zs, zs, os   ]
            ]
        )
    else:
        zs = np.zeros((share.shape[0],))
        os = np.ones((share.shape[0],))
        ss = share.T
        mat = np.array(
            [
                [os, zs, zs, ss[0]],
                [zs, os, zs, ss[1]],
                [zs, zs, os, ss[2]],
                [zs, zs, zs, os   ]
            ]
        ).T
    return mat

#######################################################################################################################
#
#                                                 affine_matrix
#

def affine_matrix(tmat, shift):
    """Creates an affine transformation matrix from a 3x3 transformation matrix or set of matrices and a shift or set of vecs

    :param tmat: base transformation matrices
    :type tmat: np.ndarray
    :param shift:
    :type shift:
    :return:
    :rtype:
    """
    base_mat = np.asarray(tmat)
    if shift is None:
        return base_mat

    if len(base_mat.shape) > 2:
        shifts = np.asarray(shift)
        if len(shifts.shape) == 1:
            shifts = np.repeat(shifts, len(base_mat))
        ones = np.ones((len(base_mat), 1))
        shifts = np.concatenate((shifts, ones), axis=1)
        zeros = np.zeros((len(base_mat), 1, 3))
        mat = np.concatenate((base_mat, zeros), axis=1)
        shifts = np.reshape(shifts, shifts.shape + (1,))
        mat = np.concatenate((mat, shifts), axis=2)

    else:
        base_shift = np.append(np.array(shift), [1])
        np.reshape(base_shift, (4, 1))
        mat = np.append(
            np.append(base_mat, np.zeros((1, 3)), axis=0),
            base_shift.transpose(),
            axis=1
        )
    return mat