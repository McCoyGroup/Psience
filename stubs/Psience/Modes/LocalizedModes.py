"""
Provides a general class of localized modes using a potentially non-linear
transformation of normal modes
"""
import numpy as np, scipy
import McUtils.Numputils as nput
from .MixtureModes import *
from .NormalModes import *
__all__ = ['LocalizedModes']

class LocalizedModes(MixtureModes):

    def __init__(self, normal_modes: MixtureModes, transformation, inverse=None, origin=None, masses=None, freqs=None, mass_weighted=None, frequency_scaled=None, **etc):
        """
        **LLM Docstring**

        Build a `LocalizedModes` object by applying a (possibly non-square) linear transformation to an existing set of normal modes, storing the transformation and the base modes for later use (e.g. by `modify`/`apply_transformation`/`get_complement`).

        :param normal_modes: the base mode set to localize
        :type normal_modes: MixtureModes
        :param transformation: the localizing transformation to apply to the mode-by-coordinates matrix, either a single 2D array (with its transpose used as the inverse unless `inverse` is given) or a `(forward, inverse)` pair
        :type transformation: np.ndarray | tuple
        :param inverse: an explicit inverse for `transformation`, overriding the default
        :type inverse: np.ndarray | None
        :param origin: the reference geometry; defaults to `normal_modes.origin`
        :type origin: np.ndarray | None
        :param masses: the atomic masses; defaults to `normal_modes.masses`
        :type masses: np.ndarray | None
        :param freqs: explicit frequencies to associate with the localized modes; if `None`, computed lazily on first access
        :type freqs: np.ndarray | None
        :param mass_weighted: accepted for interface consistency but not used directly (mass-weighting is taken from `normal_modes.mass_weighted`)
        :type mass_weighted: bool | None
        :param frequency_scaled: whether the localized modes are frequency-scaled
        :type frequency_scaled: bool | None
        :param etc: extra options forwarded to the base `MixtureModes.__init__`
        :type etc: dict
        :return: None
        :rtype: None
        """
        ...

    def __getitem__(self, item):
        """
        Takes a slice of the modes
        :param item:
        :type item:
        :return:
        :rtype:
        """
        ...

    @property
    def freqs(self):
        """
        **LLM Docstring**

        Property getter/setter for the localized modes' frequencies. The getter lazily computes them via `compute_freqs` if not already cached; the setter is currently a no-op (frequencies can't be set directly).

        :param freqs: (setter only) ignored; setting `freqs` directly is not supported
        :type freqs: np.ndarray
        :return: (getter) the (cached or newly computed) frequencies
        :rtype: np.ndarray
        """
        ...

    @freqs.setter
    def freqs(self, freqs):
        """
        **LLM Docstring**

        Property getter/setter for the localized modes' frequencies. The getter lazily computes them via `compute_freqs` if not already cached; the setter is currently a no-op (frequencies can't be set directly).

        :param freqs: (setter only) ignored; setting `freqs` directly is not supported
        :type freqs: np.ndarray
        :return: (getter) the (cached or newly computed) frequencies
        :rtype: np.ndarray
        """
        ...

    @property
    def mass_weighted(self):
        """
        **LLM Docstring**

        Property getter/setter for whether the modes are mass-weighted, delegated to (and required to match) the base modes' own `mass_weighted` state, since a `LocalizedModes` can't independently change mass-weighting.

        :param new: (setter only) the requested mass-weighting state; must match `self.base_modes.mass_weighted`
        :type new: bool
        :return: (getter) `self.base_modes.mass_weighted`
        :rtype: bool
        :raises ValueError: if the setter is given a value that doesn't match the base modes' mass-weighting state
        """
        ...

    @mass_weighted.setter
    def mass_weighted(self, new):
        """
        **LLM Docstring**

        Property getter/setter for whether the modes are mass-weighted, delegated to (and required to match) the base modes' own `mass_weighted` state, since a `LocalizedModes` can't independently change mass-weighting.

        :param new: (setter only) the requested mass-weighting state; must match `self.base_modes.mass_weighted`
        :type new: bool
        :return: (getter) `self.base_modes.mass_weighted`
        :rtype: bool
        :raises ValueError: if the setter is given a value that doesn't match the base modes' mass-weighting state
        """
        ...

    @property
    def g_matrix(self):
        """
        **LLM Docstring**

        Property getter/setter for the G-matrix, delegated to `self.base_modes.g_matrix`. The setter is currently a no-op.

        :param g: (setter only) ignored
        :type g: np.ndarray
        :return: (getter) `self.base_modes.g_matrix`
        :rtype: np.ndarray | None
        """
        ...

    @g_matrix.setter
    def g_matrix(self, g):
        """
        **LLM Docstring**

        Property getter/setter for the G-matrix, delegated to `self.base_modes.g_matrix`. The setter is currently a no-op.

        :param g: (setter only) ignored
        :type g: np.ndarray
        :return: (getter) `self.base_modes.g_matrix`
        :rtype: np.ndarray | None
        """
        ...

    def modify(self, base_modes=None, *, transformation=None, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None):
        """
        **LLM Docstring**

        Build a new `LocalizedModes` with the given fields overridden, defaulting each unspecified field to this object's own value (base modes, localizing transformation/inverse, name), and re-normalizing a `(forward, inverse)`-tuple `transformation` argument the same way the constructor does.

        :param base_modes: replacement base mode set; defaults to `self.base_modes`
        :type base_modes: MixtureModes | None
        :param transformation: replacement localizing transformation; defaults to `self.localizing_transformation[0]`
        :type transformation: np.ndarray | tuple | None
        :param freqs: replacement frequencies
        :type freqs: np.ndarray | None
        :param origin: replacement reference geometry
        :type origin: np.ndarray | None
        :param masses: replacement masses
        :type masses: np.ndarray | None
        :param inverse: replacement inverse transformation; defaults to `self.localizing_transformation[1]`
        :type inverse: np.ndarray | None
        :param name: replacement name; defaults to `self.name`
        :type name: str | None
        :param mass_weighted: accepted but not forwarded to the new instance's constructor
        :type mass_weighted: bool | None
        :param frequency_scaled: replacement frequency-scaling flag
        :type frequency_scaled: bool | None
        :param g_matrix: replacement G-matrix
        :type g_matrix: np.ndarray | None
        :return: the new, modified `LocalizedModes`
        :rtype: LocalizedModes
        """
        ...

    def _frequency_scaling(self, freqs=None):
        """
        **LLM Docstring**

        Compute the per-mode scaling factor (`sign(freq) * sqrt(|freq|)`) used to convert between frequency-scaled and non-frequency-scaled localized-mode bases.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.local_freqs`
        :type freqs: np.ndarray | None
        :return: `(freqs, conv)` -- the frequencies used and their corresponding scaling factors
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def make_mass_weighted(self, **kwargs):
        """
        **LLM Docstring**

        Build a mass-weighted version of these localized modes by mass-weighting the underlying base modes and re-localizing on top of them.

        :param kwargs: extra options forwarded to `self.base_modes.make_mass_weighted`
        :type kwargs: dict
        :return: the mass-weighted localized modes
        :rtype: LocalizedModes
        """
        ...

    def remove_mass_weighting(self, **kwargs):
        """
        **LLM Docstring**

        Build a non-mass-weighted version of these localized modes by removing the mass-weighting from the underlying base modes and re-localizing on top of them.

        :param kwargs: extra options forwarded to `self.base_modes.remove_mass_weighting`
        :type kwargs: dict
        :return: the non-mass-weighted localized modes
        :rtype: LocalizedModes
        """
        ...

    def make_frequency_scaled(self, freqs=None, **kwargs):
        """
        **LLM Docstring**

        Build a frequency-scaled (dimensionless) version of these localized modes by rescaling the localizing transformation by the per-mode frequency-scaling factors; returns `self` unchanged if already frequency-scaled.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.local_freqs`
        :type freqs: np.ndarray | None
        :param kwargs: accepted but not used in this method's body
        :type kwargs: dict
        :return: the frequency-scaled localized modes
        :rtype: LocalizedModes
        """
        ...

    def remove_frequency_scaling(self, freqs=None, **kwargs):
        """
        **LLM Docstring**

        Build a non-frequency-scaled version of these localized modes by undoing the frequency-scaling factor in the localizing transformation; returns `self` unchanged if not already frequency-scaled.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.local_freqs`
        :type freqs: np.ndarray | None
        :param kwargs: accepted but not used in this method's body
        :type kwargs: dict
        :return: the non-frequency-scaled localized modes
        :rtype: LocalizedModes
        """
        ...

    def compute_hessian(self, system='modes'):
        """
        **LLM Docstring**

        Compute the Hessian in either the localized-mode basis (via the localizing transformation and the base modes' frequencies) or the Cartesian/coordinate basis (delegating to the base class).

        :param system: which basis to compute the Hessian in, `'modes'` or `'coords'`
        :type system: str
        :return: the Hessian matrix
        :rtype: np.ndarray
        :raises ValueError: if `system` is neither `'modes'` nor `'coords'`
        """
        ...

    def apply_transformation(self, transformation, inverse=None, **opts):
        """
        **LLM Docstring**

        Apply an additional linear transformation on top of the existing localizing transformation, returning a new `LocalizedModes` built from the same base modes with the combined transformation.

        :param transformation: the additional transformation to apply, either a single 2D array (transpose used as inverse unless `inverse` is given) or a `(forward, inverse)` pair
        :type transformation: np.ndarray | tuple
        :param inverse: an explicit inverse for `transformation`
        :type inverse: np.ndarray | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the new `LocalizedModes` with the combined transformation
        :rtype: LocalizedModes
        """
        ...

    def get_complement(self, concatenate=False):
        """
        **LLM Docstring**

        Build the mode set complementary to this localized set (i.e. the remaining degrees of freedom not spanned by the localizing transformation), either as a freshly re-diagonalized (locally harmonic) mode set, or -- if `concatenate` -- as a single combined mode set spanning both this localized subspace and its complement.

        :param concatenate: whether to return the complementary modes standalone, or concatenated together with this localized set into one combined mode set spanning the full space
        :type concatenate: bool
        :return: the complementary modes, or (if `concatenate`) the combined full-space mode set
        :rtype: MixtureModes
        """
        ...