from .AffineTransform import AffineTransform
import numpy as np

######################################################################################################
##
##                                   TranslationTransform Class
##
######################################################################################################

class TranslationTransform(AffineTransform):
    """A simple TranslationTransform from the basic AffineTransformation class

    """

    def __init__(self, shift):
        """
        :param shift: 3D vector for translation
        :type shift: np.ndarray
        """

        super().__init__(np.eye(3), shift=shift)
