
import abc

__all__ = [
]
class VibronicState:
    """
    Abstract base class for representing and electronic-vibrational product state
    """
    @abc.abstractmethod
    def get_overlaps(self):
        ...


class VibronicModel:
    """
    Provides an abstraction for calculating vibronic spectra
    """

    def __init__(self, *vibronic_states:'VibronicState', **opts):
        self.states = vibronic_states
        self.opts = opts

    def calculate_spectrum(self):
        ...