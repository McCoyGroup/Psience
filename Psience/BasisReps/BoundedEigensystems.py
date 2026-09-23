"""Eigenvalue bounds assembled from spectra of additive Hermitian operators."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from McUtils.Plots import StairsPlot

__all__ = [
    "WeylEigenvalueBounds",
    "WeylBoundedEigenspace",
]


@dataclass(frozen=True)
class WeylEigenvalueBounds:
    """Lower and upper bounds for an increasingly ordered eigenspectrum."""

    lower: np.ndarray
    upper: np.ndarray

    @property
    def intervals(self) -> np.ndarray:
        """Return an ``(n, 2)`` array of lower/upper intervals."""

        return np.column_stack((self.lower, self.upper))

    @property
    def widths(self) -> np.ndarray:
        """Return the width of every bounded eigenvalue interval."""

        return self.upper - self.lower

    def __iter__(self):
        """Allow unpacking as ``lower, upper = bounds``."""

        yield self.lower
        yield self.upper


class WeylBoundedEigenspace:
    """Weyl-bounded spectrum of a sum of ``k`` Hermitian operators.

    Parameters
    ----------
    eigenspectra
        A collection of ``k`` equal-length one-dimensional eigenspectra.
        Input spectra need not already be sorted.
    eigenvectors
        Optional opaque eigensystem payload.  It is retained without coercion
        or interpretation, but is not used to compute eigenvalue bounds.  In
        particular, :meth:`from_dvrs` stores its input ``DVRWavefunctions``
        objects here for later approximate-eigenvector dispatch.
    labels
        Optional labels for the additive operator components.

    Notes
    -----
    Bounds for more than two spectra are evaluated by repeated max-plus and
    min-plus Weyl convolution.  For increasingly ordered spectra ``a`` and
    ``b``, the binary update is

    ``lower[m] = max_i(a[i] + b[m-i])``

    and

    ``upper[m] = min_i(a[i] + b[n-1+m-i])``.

    Repeating these updates gives valid generalized Weyl bounds for a sum of
    any number of equally sized Hermitian matrices without constructing them.
    """

    def __init__(
        self,
        eigenspectra: Iterable[Sequence[float]],
        *,
        eigenvectors: Iterable | None = None,
        labels: Sequence[str] | None = None,
    ):
        spectra = tuple(
            self._canonicalize_spectrum(spectrum, index)
            for index, spectrum in enumerate(eigenspectra)
        )
        if len(spectra) == 0:
            raise ValueError("at least one eigenspectrum is required")

        dimension = spectra[0].size
        if any(spectrum.size != dimension for spectrum in spectra[1:]):
            raise ValueError(
                "all eigenspectra must have the same length; got "
                f"{tuple(spectrum.size for spectrum in spectra)}"
            )

        if labels is None:
            labels = tuple(f"component {i}" for i in range(len(spectra)))
        else:
            labels = tuple(labels)
            if len(labels) != len(spectra):
                raise ValueError("labels must have one entry per eigenspectrum")

        if eigenvectors is not None:
            eigenvectors = tuple(eigenvectors)

        self._eigenspectra = spectra
        self._eigenvectors = eigenvectors
        self._labels = labels
        self._bounds = None
        self._component_grids = None
        self._composite_grid = None
        self._coupling_values = None
        self._coupling_potential = None

    @staticmethod
    def _wavefunction_grid(wavefunctions):
        """Build the coordinate-eigenvalue grid for one wavefunction basis."""

        grid = np.asarray(wavefunctions.grid)
        coefficients = np.asarray(wavefunctions.wavefunctions)
        if grid.ndim == 1:
            flat_grid = grid[:, np.newaxis]
        else:
            flat_grid = grid.reshape(-1, grid.shape[-1])
        if coefficients.ndim != 2 or coefficients.shape[0] != flat_grid.shape[0]:
            raise ValueError(
                "DVR wavefunction coefficients and grid points are inconsistent: "
                f"got {coefficients.shape} and {flat_grid.shape}"
            )

        coordinate_matrices = np.einsum(
            "pi,pd,pj->dij",
            coefficients,
            flat_grid,
            coefficients,
            optimize=True,
        )
        return np.linalg.eigvalsh(coordinate_matrices).T

    @classmethod
    def from_dvrs(cls, dvrs, coupling_potential):
        """Construct bounds from DVR eigensystems and a coupling potential.

        Each input ``DVRWavefunctions`` object supplies one component energy
        spectrum and wavefunction basis.  Coordinate operators are projected
        into those bases and diagonalized to form contracted component grids.
        Their Cartesian product defines the composite grid passed to
        ``coupling_potential`` as an ``(n_points, n_coordinates)`` array.

        The bounded sum contains two spectra: the direct-product sum of the
        component DVR energies and the sorted coupling-potential values.  The
        original wavefunction objects are stored directly in ``eigenvectors``.
        """

        wavefunctions = tuple(dvrs)
        if len(wavefunctions) == 0:
            raise ValueError("at least one DVRWavefunctions object is required")
        if not callable(coupling_potential):
            raise TypeError("coupling_potential must be callable")

        component_energies = []
        component_grids = []
        for index, wfns in enumerate(wavefunctions):
            if not hasattr(wfns, "energies") or not hasattr(wfns, "wavefunctions"):
                raise TypeError(
                    f"dvrs[{index}] does not provide DVR wavefunction data"
                )
            energies = cls._canonicalize_spectrum(wfns.energies, index)
            if energies.size != np.asarray(wfns.wavefunctions).shape[-1]:
                raise ValueError(
                    f"dvrs[{index}] has {energies.size} energies but "
                    f"{np.asarray(wfns.wavefunctions).shape[-1]} wavefunctions"
                )
            component_energies.append(energies)
            component_grids.append(cls._wavefunction_grid(wfns))

        product_shape = tuple(energies.size for energies in component_energies)
        product_indices = np.meshgrid(
            *[np.arange(size) for size in product_shape],
            indexing="ij",
        )
        uncoupled_energies = np.zeros(product_shape, dtype=float)
        composite_grid_parts = []
        for energies, grid, indices in zip(
            component_energies,
            component_grids,
            product_indices,
        ):
            uncoupled_energies += energies[indices]
            composite_grid_parts.append(grid[indices])
        composite_grid = np.concatenate(composite_grid_parts, axis=-1)
        flat_composite_grid = composite_grid.reshape(-1, composite_grid.shape[-1])

        coupling_values = np.asarray(coupling_potential(flat_composite_grid))
        if coupling_values.ndim == 0:
            coupling_values = np.full(flat_composite_grid.shape[0], coupling_values)
        coupling_values = np.reshape(coupling_values, -1)
        if coupling_values.size != flat_composite_grid.shape[0]:
            raise ValueError(
                "coupling_potential must return one value per composite grid point; "
                f"got {coupling_values.shape} for {flat_composite_grid.shape[0]} points"
            )
        if np.iscomplexobj(coupling_values) and not np.allclose(coupling_values.imag, 0):
            raise ValueError("coupling_potential returned complex values")
        coupling_values = np.asarray(coupling_values.real, dtype=float)
        if not np.all(np.isfinite(coupling_values)):
            raise ValueError("coupling_potential returned non-finite values")

        space = cls(
            [uncoupled_energies.ravel(), coupling_values],
            eigenvectors=wavefunctions,
            labels=("uncoupled DVR product", "coupling potential"),
        )
        space._component_grids = tuple(component_grids)
        space._composite_grid = composite_grid
        space._coupling_values = coupling_values.reshape(product_shape)
        space._coupling_potential = coupling_potential
        return space

    @staticmethod
    def _canonicalize_spectrum(values, index):
        spectrum = np.asarray(values)
        if spectrum.ndim != 1:
            raise ValueError(
                f"eigenspectra[{index}] must be one-dimensional; "
                f"got shape {spectrum.shape}"
            )
        if spectrum.size == 0:
            raise ValueError(f"eigenspectra[{index}] must not be empty")
        if np.iscomplexobj(spectrum):
            if not np.allclose(spectrum.imag, 0):
                raise ValueError(f"eigenspectra[{index}] contains complex eigenvalues")
            spectrum = spectrum.real
        spectrum = np.asarray(spectrum, dtype=float)
        if not np.all(np.isfinite(spectrum)):
            raise ValueError(f"eigenspectra[{index}] contains non-finite eigenvalues")
        spectrum = np.sort(spectrum).copy()
        spectrum.flags.writeable = False
        return spectrum

    @property
    def eigenspectra(self) -> tuple[np.ndarray, ...]:
        return self._eigenspectra

    @property
    def eigenvectors(self) -> tuple | None:
        """Optional component eigenvectors retained for future constructions."""

        return self._eigenvectors

    @property
    def component_grids(self) -> tuple[np.ndarray, ...] | None:
        """Contracted coordinate grids generated by :meth:`from_dvrs`."""

        return self._component_grids

    @property
    def composite_grid(self) -> np.ndarray | None:
        """Composite coordinate grid generated by :meth:`from_dvrs`."""

        return self._composite_grid

    @property
    def coupling_values(self) -> np.ndarray | None:
        """Coupling values on ``composite_grid``, when constructed from DVRs."""

        return self._coupling_values

    @property
    def labels(self) -> tuple[str, ...]:
        return self._labels

    @property
    def dimension(self) -> int:
        return self._eigenspectra[0].size

    @property
    def num_components(self) -> int:
        return len(self._eigenspectra)

    @staticmethod
    def _extend_bounds(
        bounds: WeylEigenvalueBounds,
        spectrum: np.ndarray,
    ) -> WeylEigenvalueBounds:
        """Add one spectrum by binary Weyl max/min convolution."""

        dimension = spectrum.size
        lower = np.empty(dimension, dtype=float)
        upper = np.empty(dimension, dtype=float)
        for level in range(dimension):
            lower[level] = np.max(
                bounds.lower[: level + 1] + spectrum[: level + 1][::-1]
            )
            upper[level] = np.min(
                bounds.upper[level:] + spectrum[level:][::-1]
            )
        return WeylEigenvalueBounds(lower=lower, upper=upper)

    def get_bounds(self) -> WeylEigenvalueBounds:
        """Compute generalized Weyl bounds for the summed eigensystem."""

        first = self._eigenspectra[0].copy()
        bounds = WeylEigenvalueBounds(lower=first, upper=first.copy())
        for spectrum in self._eigenspectra[1:]:
            bounds = self._extend_bounds(bounds, spectrum)

        midpoint = (bounds.lower + bounds.upper) / 2
        crossed = bounds.lower > bounds.upper
        close_crossings = crossed & np.isclose(
            bounds.lower,
            bounds.upper,
            rtol=1e-12,
            atol=1e-14,
        )
        bounds.lower[close_crossings] = midpoint[close_crossings]
        bounds.upper[close_crossings] = midpoint[close_crossings]
        if np.any(bounds.lower > bounds.upper):
            raise RuntimeError("internally inconsistent Weyl bounds")
        bounds.lower.flags.writeable = False
        bounds.upper.flags.writeable = False
        return bounds

    @property
    def bounds(self) -> WeylEigenvalueBounds:
        """Cached generalized Weyl bounds for the summed eigensystem."""

        if self._bounds is None:
            self._bounds = self.get_bounds()
        return self._bounds

    def plot(
        self,
        *,
        n_levels: int | None = None,
        energy_shift: float = 0,
        figure=None,
        lower_style: dict | None = None,
        upper_style: dict | None = None,
        plot_label: str = "Weyl-bounded eigenspectrum",
        **plot_options,
    ):
        """Plot lower and upper bounds as unfilled staircases."""

        if n_levels is None:
            n_levels = self.dimension
        if not 1 <= n_levels <= self.dimension:
            raise ValueError(f"n_levels must be between 1 and {self.dimension}")

        lower_style = dict(
            {"color": "#0077b6", "linewidth": 1.5, "label": "lower bound"},
            **({} if lower_style is None else lower_style),
        )
        upper_style = dict(
            {"color": "#d62828", "linewidth": 1.5, "label": "upper bound"},
            **({} if upper_style is None else upper_style),
        )
        edges = np.arange(n_levels + 1) - 0.5
        lower = self.bounds.lower[:n_levels] - energy_shift
        upper = self.bounds.upper[:n_levels] - energy_shift

        if figure is None:
            figure = StairsPlot(
                lower,
                edges,
                plot_legend=True,
                axes_labels=["Ordered eigenvalue index", "Energy"],
                plot_label=plot_label,
                **lower_style,
                **plot_options,
            )
        else:
            StairsPlot(lower, edges, figure=figure, **lower_style)
        StairsPlot(upper, edges, figure=figure, **upper_style)
        return figure

    def approximate_eigenvectors(self, *args, **kwargs):
        """Reserved for a later Ritz/subspace eigenvector construction."""

        if self._eigenvectors is None:
            raise ValueError(
                "component eigenvectors were not supplied when this eigenspace was built"
            )
        raise NotImplementedError(
            "approximate coupled eigenvectors require coupling-operator matrix "
            "elements and will be implemented separately"
        )

    def to_state(self, serializer=None):
        """Return the state needed to reconstruct this bounded eigenspace."""

        return {
            "eigenspectra": self.eigenspectra,
            "eigenvectors": self.eigenvectors,
            "labels": self.labels,
        }

    @classmethod
    def from_state(cls, state, serializer=None):
        """Reconstruct an eigenspace from :meth:`to_state` output."""

        deserialize = (lambda value: value) if serializer is None else serializer.deserialize
        return cls(
            deserialize(state["eigenspectra"]),
            eigenvectors=deserialize(state.get("eigenvectors")),
            labels=deserialize(state.get("labels")),
        )

    def __len__(self):
        return self.dimension

    def __repr__(self):
        return (
            f"{type(self).__name__}(components={self.num_components}, "
            f"dimension={self.dimension})"
        )
