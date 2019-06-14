import numpy as np
from .TransformationFunction import TransformationFunction
from .TransformationUtilities import affine_matrix
from .TransformationUtilities import merge_transformation_matrix

######################################################################################################
##
##                                   AffineTranform Class
##
######################################################################################################

class AffineTransform(TransformationFunction):
    """A simple AffineTranform implementation of the TransformationFunction abstract base class

    """

    def __init__(self, tmat, shift=None):
        """tmat must be a transformation matrix to work properly

        :param shift: the shift for the transformation
        :type shift: np.ndarray
        :param tmat: the matrix for the linear transformation
        :type tmat: np.ndarray
        """

        self.transf = affine_matrix(tmat, shift)
        super().__init__()

    @property
    def transform(self):
        return self.transf[:3, :3]

    @property
    def shift(self):
        transf = self.transf
        transf_shape = transf.shape
        if transf_shape[-1] == 4:
            vec = transf[:-1, 4] # hopefully this copies rather than just being a view...
        else:
            vec = None
        return vec

    def merge(self, other):
        """

        :param other:
        :type other: np.ndarray or AffineTransform
        """

        if isinstance(type(other), AffineTransform):
            other = other.transf
        transf = self.transf

        # I wanted to use type(self) but then I realized that'll fuck me over
        # if I want to merge like a ScalingTransform and a RotationTransform
        return AffineTransform(merge_transformation_matrix(transf, other))


    def reverse(self):
        """Inverts the matrix

        :return:
        :rtype:
        """
        
        inverse = np.linalg.inv(self.transf)
        return AffineTransform(inverse)


    def operate(self, coords):
        """

        :param coords: the array of coordinates passed in
        :type coords: np.ndarry
        """

        coords = np.array(coords)
        coord_shape = coords.shape
        if len(coord_shape) == 1:
            coords.reshape((1, coord_shape[1]))
        elif len(coord_shape) > 2:
            nels = np.product(coord_shape[:-1])
            coords.reshape((nels, 3))

        tmat = self.transf
        if tmat.shape[-1] == 4:
            adj_coord = np.append(coords, np.ones(coords.shape[0], 1), axis=0)
            adj_coord = np.dot(tmat, adj_coord)
            adj_coord = adj_coord[:, :3]
        else:
            adj_coord = np.dot(tmat, coords)

        adj_coord.reshape(coord_shape)

        return adj_coord

    def __repr__(self):
        ## we'll basically just leverage the ndarray repr:
        return "{}({})".format(type(self), str(self.transf))