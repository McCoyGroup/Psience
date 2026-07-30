import numpy as np, scipy.linalg as slag, collections
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
from .MixtureModes import MixtureModes
__all__ = ['NormalModes', 'ReactionPathModes']
__reload_hook__ = ['.MixtureModes']

class NormalModes(MixtureModes):
    name = 'NormalModes'

    def __init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=False, frequency_scaled=False, g_matrix=None):
        """
        **LLM Docstring**

        Build a `NormalModes` object (a `MixtureModes` specialized for harmonic normal modes) from a mode basis, coefficient matrix, and associated frequency/mass/origin data.

        :param basis: the coordinate system the modes are expressed in
        :type basis: CoordinateSystem
        :param coeffs: the mode coefficient (coordinates-by-modes) matrix
        :type coeffs: np.ndarray
        :param freqs: the vibrational frequencies associated with each mode
        :type freqs: np.ndarray | None
        :param origin: the reference geometry the modes are expanded about
        :type origin: np.ndarray | None
        :param masses: the atomic masses
        :type masses: np.ndarray | None
        :param inverse: the modes-by-coordinates (inverse) matrix
        :type inverse: np.ndarray | None
        :param name: a name for the mode set
        :type name: str | None
        :param mass_weighted: whether the modes are mass-weighted
        :type mass_weighted: bool
        :param frequency_scaled: whether the modes are frequency-scaled (dimensionless)
        :type frequency_scaled: bool
        :param g_matrix: the G-matrix associated with the modes
        :type g_matrix: np.ndarray | None
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_fg(cls, basis, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, zero_freq_cutoff=None, mass_weighted=None, origin=None, projector=None, **opts):
        """
        Generates normal modes from the specified F and G matrices

        :param basis:
        :param f_matrix: second derivatives of the potential
        :param mass_spec:
        :param mass_units:
        :param remove_transrot:
        :param opts:
        :return:
        """
        ...
    default_projected_zero_freq_cutoff = None

    @classmethod
    def from_molecule(cls, mol, dimensionless=False, use_internals=None, potential_derivatives=None, project_transrot=True, zero_freq_cutoff=None, masses=None, energy_evaluator=None, **opts):
        """
        **LLM Docstring**

        Build normal modes for a molecule, dispatching to `from_fg` with either the internal-coordinate Hessian/G-matrix (if `use_internals`) or the (optionally translation/rotation-projected) Cartesian Hessian/masses, computing missing potential derivatives via the molecule's energy evaluator or falling back to its existing normal modes where possible.

        :param mol: the molecule to build normal modes for
        :type mol: Molecule
        :param dimensionless: whether the resulting modes should be dimensionless (frequency-scaled)
        :type dimensionless: bool
        :param use_internals: whether to build the modes in internal coordinates; defaults to whether the molecule has internal coordinates defined
        :type use_internals: bool | None
        :param potential_derivatives: explicit potential derivative tensors to use instead of the molecule's own/computed ones
        :type potential_derivatives: list[np.ndarray] | None
        :param project_transrot: whether to project out translational/rotational degrees of freedom before diagonalizing (Cartesian case only)
        :type project_transrot: bool
        :param zero_freq_cutoff: the frequency cutoff below which modes are discarded as translation/rotation; defaults to `default_projected_zero_freq_cutoff` or 1 wavenumber when projecting
        :type zero_freq_cutoff: float | None
        :param masses: masses to use instead of the molecule's own
        :type masses: np.ndarray | None
        :param energy_evaluator: an explicit energy evaluator to use for computing missing potential derivatives
        :type energy_evaluator: object | None
        :param opts: extra options forwarded to `from_fg`
        :type opts: dict
        :return: the constructed normal modes
        :rtype: NormalModes
        :raises ValueError: if no potential derivatives are available and none can be computed or reused
        """
        ...

    def get_reaction_path_modes(self, grad: np.ndarray, zero_gradient_cutoff=None, return_status=False, **kwargs):
        """
        **LLM Docstring**

        Build reaction-path-following normal modes from this mode set and a gradient vector, via `ReactionPathModes.from_modes_and_grad`.

        :param grad: the gradient vector defining the reaction path direction
        :type grad: np.ndarray
        :param zero_gradient_cutoff: the gradient-magnitude cutoff below which the gradient direction is treated as zero (falling back to ordinary normal modes)
        :type zero_gradient_cutoff: float | None
        :param return_status: whether to also return whether the gradient direction was actually used
        :type return_status: bool
        :param kwargs: extra options forwarded to `ReactionPathModes.from_modes_and_grad`
        :type kwargs: dict
        :return: the reaction-path modes, or `(modes, status)` if `return_status` is set
        :rtype: ReactionPathModes | tuple
        """
        ...

class ReactionPathModes(NormalModes):
    zero_gradient_cutoff = 2.5e-05

    @classmethod
    def get_rp_modes(cls, gradient, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, mass_weighted=None, zero_freq_cutoff=None, return_gmatrix=False, projector=None, zero_gradient_cutoff=None, use_max_gradient_cutoff=True, gradient_check_transformation=None, return_indices=False):
        """
        **LLM Docstring**

        Compute reaction-path-following vibrational modes from a gradient and force-constant/mass specification: mass-weights and normalizes the gradient to form the leading "reaction coordinate" mode, projects the remaining Hessian (and any existing constraint projector) orthogonal to that direction to get the transverse vibrational modes via `get_normal_modes`, and (unless `remove_transrot`) drops the transverse mode most similar to the (near-zero-frequency) projected-out translation/rotation modes so the gradient mode takes its place. Points whose gradient magnitude falls below `zero_gradient_cutoff` are instead treated as ordinary stationary-point normal modes with no reaction coordinate.

        :param gradient: the gradient vector(s) defining the reaction-path direction(s)
        :type gradient: np.ndarray
        :param f_matrix: the force-constant (Hessian) matrix/matrices
        :type f_matrix: np.ndarray
        :param mass_spec: the atomic masses (or element symbols), or an already-built G-matrix/mass matrix
        :type mass_spec: np.ndarray | list[str]
        :param remove_transrot: whether to remove translation/rotation modes as usual (dropping the one most similar to the reaction-coordinate mode)
        :type remove_transrot: bool
        :param dimensionless: whether the resulting modes should be frequency-scaled (dimensionless)
        :type dimensionless: bool
        :param mass_weighted: whether the resulting modes should be mass-weighted
        :type mass_weighted: bool | None
        :param zero_freq_cutoff: frequency cutoff for discarding near-zero-frequency modes
        :type zero_freq_cutoff: float | None
        :param return_gmatrix: whether to also return the (broadcast) mass/G-matrix used
        :type return_gmatrix: bool
        :param projector: an existing constraint projector to combine with the gradient-orthogonal projector
        :type projector: np.ndarray | None
        :param zero_gradient_cutoff: the gradient-magnitude cutoff below which a point is treated as having no reaction-path direction; defaults to `cls.zero_gradient_cutoff`
        :type zero_gradient_cutoff: float | None
        :param use_max_gradient_cutoff: whether to test the cutoff against the maximum gradient component rather than its norm
        :type use_max_gradient_cutoff: bool
        :param gradient_check_transformation: an optional coordinate-transformation expansion used to re-express the gradient before testing it against the cutoff (without changing the gradient actually used to build the mode)
        :type gradient_check_transformation: list[np.ndarray] | None
        :param return_indices: whether to also return which batch positions used the reaction-path treatment versus the plain stationary-point treatment
        :type return_indices: bool
        :return: the mode data (`ModeData(freqs, modes, inv)`), or a tuple additionally including the mass/G-matrix and/or the `(regular_mode_pos, rem_pos)` index split depending on the `return_*` flags
        :rtype: object | tuple
        """
        ...

    @classmethod
    def from_grad_fg(cls, basis, gradient, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, zero_freq_cutoff=None, mass_weighted=None, origin=None, projector=None, zero_gradient_cutoff=None, gradient_check_transformation=None, return_status=False, **opts):
        """
        Generates normal modes from the specified F and G matrices

        :param basis:
        :param f_matrix: second derivatives of the potential
        :param mass_spec:
        :param mass_units:
        :param remove_transrot:
        :param opts:
        :return:
        """
        ...

    @classmethod
    def from_molecule(cls, mol, dimensionless=False, use_internals=None, potential_derivatives=None, project_transrot=True, zero_freq_cutoff=None, masses=None, zero_gradient_cutoff=None, return_status=False, gradient_check_internals=None, gradient_check_transformation=None, **opts):
        """
        **LLM Docstring**

        Build reaction-path-following normal modes for a molecule at a (generally non-stationary) point, using its gradient together with either its internal-coordinate Hessian/G-matrix (if `use_internals`) or its (optionally translation/rotation-projected) Cartesian Hessian/masses, dispatching to `from_grad_fg`.

        :param mol: the molecule to build reaction-path modes for
        :type mol: Molecule
        :param dimensionless: whether the resulting modes should be dimensionless
        :type dimensionless: bool
        :param use_internals: whether to build the modes in internal coordinates
        :type use_internals: bool | None
        :param potential_derivatives: explicit potential derivative tensors (including the gradient) to use instead of the molecule's own/computed ones
        :type potential_derivatives: list[np.ndarray] | None
        :param project_transrot: whether to project out translational/rotational degrees of freedom (Cartesian case only)
        :type project_transrot: bool
        :param zero_freq_cutoff: the frequency cutoff below which modes are discarded as translation/rotation
        :type zero_freq_cutoff: float | None
        :param masses: masses to use instead of the molecule's own
        :type masses: np.ndarray | None
        :param zero_gradient_cutoff: the gradient-magnitude cutoff below which the point is treated as stationary
        :type zero_gradient_cutoff: float | None
        :param return_status: whether to also return whether the reaction-path treatment was actually used
        :type return_status: bool
        :param gradient_check_internals: an alternate internal-coordinate specification used to build `gradient_check_transformation`, if that isn't given directly
        :type gradient_check_internals: object | None
        :param gradient_check_transformation: an optional coordinate-transformation expansion used only to test the gradient against the zero-gradient cutoff
        :type gradient_check_transformation: list[np.ndarray] | None
        :param opts: extra options forwarded to `from_grad_fg`
        :type opts: dict
        :return: the constructed reaction-path modes, or `(modes, status)` if `return_status` is set
        :rtype: ReactionPathModes | tuple
        """
        ...

    @classmethod
    def from_modes_and_grad(cls, modes: MixtureModes, grad: np.ndarray, zero_gradient_cutoff=None, use_max_gradient_cutoff=True, return_status=False, mass_weighted=None, **projection_opts):
        """
        **LLM Docstring**

        Given an existing mode set and a gradient vector, build reaction-path-following modes by localizing the modes against the (mass-weighted, normalized) gradient direction as an orthogonal projection constraint, unless the gradient falls below `zero_gradient_cutoff`, in which case the original modes are returned unchanged.

        :param modes: the mode set to re-localize around the gradient direction
        :type modes: MixtureModes
        :param grad: the gradient vector
        :type grad: np.ndarray
        :param zero_gradient_cutoff: the gradient-magnitude cutoff below which the modes are left unchanged
        :type zero_gradient_cutoff: float | None
        :param use_max_gradient_cutoff: whether to test the cutoff against the maximum gradient component rather than its norm
        :type use_max_gradient_cutoff: bool
        :param return_status: whether to also return whether the reaction-path treatment was actually applied
        :type return_status: bool
        :param mass_weighted: whether the resulting modes should be mass-weighted; defaults to `modes.mass_weighted`
        :type mass_weighted: bool | None
        :param projection_opts: extra options forwarded to `MixtureModes.localize`
        :type projection_opts: dict
        :return: the reaction-path modes, or `(modes, status)` if `return_status` is set
        :rtype: MixtureModes | tuple
        :raises ValueError: if localizing against the gradient direction doesn't actually reduce the mode-subspace size as expected
        """
        ...