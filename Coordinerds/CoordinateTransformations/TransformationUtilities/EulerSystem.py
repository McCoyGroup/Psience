"""Implements calculation of Euler angles, Euler matrices, etc

"""

import numpy as np

# I pulled a bunch of Euler matrices for diferent conventions out of Mathematica with the
# code Resources/Mathematica/eulerSolver.wl
# Any new orientation can be easily added like that
# We could do some metaprogramming thing with eval(orientation+_mat)
# but I think it's better to just hard code the choices the dispatcher maps

#### ZYZ
def zyz_mat(c, s):
    return [
        [ c[0]*c[1]*c[2]+-s[0]*s[2], -c[2]*s[0]+-c[0]*c[1]*s[2], c[0]*s[1] ],
        [ c[1]*c[2]*s[0]+c[0]*s[2] , c[0]*c[2]+-c[1]*s[0]*s[2] , s[0]*s[1] ],
        [ -c[2]*s[1]               , s[1]*s[2]                 , c[1]      ]
    ]
def zyz_angles(basis):
    return [
        -np.arccos(-basis[2][0]*1/(np.sqrt(1+-(basis[2][2])**(2)))),
        np.arctan2(basis[2][2], -np.sqrt(1+-(basis[2][2])**(2))),
        -np.arccos(basis[0][2]*1/(np.sqrt(1+-(basis[2][2])**(2))))
    ]

##### XYZ
def xyz_mat(c, s):
    return [
        [ c[1]*c[2]                , -c[1]*s[2]               , s[1]       ],
        [ c[2]*s[0]*s[1]+c[0]*s[2] , c[0]*c[2]+-s[0]*s[1]*s[2], -c[1]*s[0] ],
        [ -c[0]*c[2]*s[1]+s[0]*s[2], c[2]*s[0]+c[0]*s[1]*s[2] , c[0]*c[1]  ]
    ]
def xyz_angles(basis):
    return [
        -np.arccos(-1/(np.sqrt((basis[0][0])**(2)+(basis[0][1])**(2)))*basis[2][2]),
        -np.arccos(-np.sqrt((basis[0][0])**(2)+(basis[0][1])**(2))),
        np.arctan2(-basis[0][0], basis[0][1])
    ]

##### XZY
def xzy_mat(c, s):
    return  [
        [ c[1]*c[2]                , -s[1]    , c[1]*s[2]                 ],
        [ c[0]*c[2]*s[1]+s[0]*s[2] , c[0]*c[1], -c[2]*s[0]+c[0]*s[1]*s[2] ],
        [ c[2]*s[0]*s[1]+-c[0]*s[2], c[1]*s[0], c[0]*c[2]+s[0]*s[1]*s[2]  ]
    ]
def xzy_angles(basis):
    return [
        -np.arccos(-1/(np.sqrt(1+-(basis[0][1])**(2)))*basis[1][1]),
        np.arctan2(-np.sqrt(1+-(basis[0][1])**(2)), -basis[0][1]),
        -np.arccos(-basis[0][0]*1/(np.sqrt(1+-(basis[0][1])**(2))))
    ]

##### ZYX
def zyx_mat(c, s):
    return [
        [ c[0]*c[1], -c[2]*s[0]+c[0]*s[1]*s[2], c[0]*c[2]*s[1]+s[0]*s[2]  ],
        [ c[1]*s[0], c[0]*c[2]+s[0]*s[1]*s[2] , c[2]*s[0]*s[1]+-c[0]*s[2] ],
        [ -s[1]    , c[1]*s[2]                , c[1]*c[2]                 ]
    ]
def zyx_angles(basis):
    return [
        np.arctan2(-basis[0][0], -basis[1][0]),
        -np.arccos(-np.sqrt((basis[0][0])**(2)+(basis[1][0])**(2))),
        -np.arccos(-1/(np.sqrt((basis[0][0])**(2)+(basis[1][0])**(2)))*basis[2][2])
    ]

### Define a dispatch map for these functions

euler_mat_map = {
    "zyz" : zyz_mat,
    "xyz" : xyz_mat,
    "zyx" : zyx_mat
}

# After we let Mathematica do the heavy lifting of actually generating the matrices
# we then simply take their sines and cosines

def euler_matrix(angles, ordering="xyz"):
    """Returns the Euler matrix for the specified angles

    :param angles:
    :type angles:
    :param ordering: the order in which the rotations should be performed
    :type ordering:
    """
    ordering = ordering.lower()
    if ordering in euler_mat_map:
        mat_gen = euler_mat_map[ordering]
    else:
        raise KeyError("Euler matrix for ordering '{}' not yet supported".format(ordering))
    c = np.cos(angles)
    s = np.sin(angles)
    return mat_gen(c, s)

def euler_angles(basis, ordering="xyz"):
    """Calculates the Euler angles for the basis

    :param basis: the basis to get the Euler angles for
    :type basis: np.ndarray
    """
    if ordering in euler_mat_map:
        ang_gen = euler_mat_map[ordering]
    else:
        raise KeyError("Euler matrix for orientation '{}' not yet supported".format(ordering))
    return ang_gen(basis)