import numpy as np

def merge_transformation_matrix(transf, other):
    """Merges two transformation matrices

    :param transf:
    :type transf: np.ndarray
    :param other:
    :type other: np.ndarray
    :return:
    :rtype:
    """

    other_shape = other.shape
    self_shape = transf.shape

    if other_shape[-1] > self_shape[-1]:
        transf = np.append(np.append(transf, np.zeros((1, 3)), axis=0), np.eye(4)[:, -1], axis=1)
    elif self_shape[-1] > other_shape[-1]:
        other = np.append(np.append(other, np.zeros((1, 3)), axis=0), np.eye(4)[:, -1], axis=1)

    return np.dot(transf, other)