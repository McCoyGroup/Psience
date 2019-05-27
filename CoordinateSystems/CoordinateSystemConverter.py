from collections import OrderedDict as odict
import os, abc, numpy as np
from ..Utilities.ExtensionLoader import ExtensionLoader

######################################################################################################
##
##                                   CoordinateSystemConverters Class
##
######################################################################################################

class CoordinateSystemConverters:
    """A coordinate converter class. It's a singleton so can't be instantiated.

    """

    converters = odict([])
    converters_dir = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "Resources",
        "Converters"
    )
    converters_package = ".".join(__name__.split(".")[:-1])

    def __init__(self):
        raise NotImplementedError("{} is a singleton".format(type(self)))

    @classmethod
    def get_coordinates(self, coordinate_set):
        """Extracts coordinates from a coordinate_set
        """
        pass

    @classmethod
    def _get_converter_file(self, file):
        if os.path.exists(file):
            abspath = file
        else:
            abspath = os.path.join(self.converters_dir, file)
        return abspath

    @classmethod
    def load_converter(self, converter):

        file = self._get_converter_file(converter)
        loader = ExtensionLoader(self.converters_dir, self.converters_package)
        env = loader.load(file)

        try:
            converters = env.__converters__
        except KeyError:
            raise KeyError("converter at {} missing field '{}'".format(file, "converters"))

        for conv in converters:
            type_pair = tuple(conv.types)
            self.converters[type_pair] = conv

    @classmethod
    def _preload_converters(self):
        for file in os.listdir(self.converters_dir):
            self.load_converter(file)

    @classmethod
    def get_converter(self, system1, system2):
        """Gets the appropriate converter for two CoordinateSystem objects

        :param system1:
        :type system1: CoordinateSystem
        :param system2:
        :type system2: CoordinateSystem
        :return:
        :rtype:
        """
        try:
            converter = self.converters[(system1, system2)]
        except KeyError:
            for key_pair, conv in reversed(self.converters.items()):
                if isinstance(system1, key_pair[0]) and isinstance(system2, key_pair[1]):
                    converter = conv
                    break
            else:
                raise KeyError(
                    "no rules for converting coordinate system {} to {}".format(system1, system2)
                )
        return converter

######################################################################################################
##
##                                   CoordinateSystemConverter Class
##
######################################################################################################
class CoordinateSystemConverter(metaclass=abc.ABCMeta):
    """A base class for type converters
    """

    @property
    @abc.abstractmethod
    def types(self):
        """The types property of a converter returns the types the converter converts

        """
        pass

    def convert_many(self, coords_list, **kwargs):
        """Converts many coordinates. Used in cases where a CoordinateSet has higher dimension
        than its basis dimension. Should be overridden by a converted to provide efficient conversions
        where necessary.

        :param coords_list: many sets of coords
        :type coords_list: np.ndarray
        :param kwargs:
        :type kwargs:
        """
        return np.array((self.convert(coords, **kwargs) for coords in coords_list))

    @abc.abstractmethod
    def convert(self, coords, **kwargs):
        """The main necessary implementation method for a converter class.
        Provides the actual function that converts the coords set

        :param coords:
        :type coords: np.ndarray
        :param kwargs:
        :type kwargs:
        """
        pass