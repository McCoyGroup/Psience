"""The coordinate mats class defines an architecture to mats coordinates

"""
import numpy as np

######################################################################################################
##
##                                   CoordinateTransform Class
##
######################################################################################################

class CoordinateTransform:
    """The CoordinateTransform class provides a simple, general way to represent a
    compound coordinate transformation

    :transforms:
    """

    def __init__(self, transforms):
        self.transform_list = [self.parse_transform(tf) for tf in transforms]
        self.condense_transforms()

    def condense_transforms(self):
        raise NotImplemented

    @staticmethod
    def parse_transform(tf):
        return NotImplemented

    @staticmethod
    def matrix_transform(mat):
        return np.array(mat)

