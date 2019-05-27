from .AffineTransform import AffineTransform
import numpy as np

######################################################################################################
##
##                                   ScalingTransform Class
##
######################################################################################################

class ScalingTransform(AffineTransform):
    """A simple ScalingTransform from the basic AffineTransformation class

    """

    def __init__(self, scalings):
        """creates an AffineTransform

        :param scalings:
        :type scalings:
        """
        super().__init__(np.diag(scalings), shift=None)
