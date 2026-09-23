"""Draft path-resolved density-matrix propagation for :mod:`Psience.Nonlinear`.

This module is intended to live beside ``NonlinearResponse.py``.  It reuses
the same transition-data preparation and ``prep_liouville_spaces`` pathway
construction as the response generator.  The complete density-matrix history
is propagated with unitary pulse kicks.  Each pathway also carries a masked
matrix history restricted to the density-matrix elements admitted by its
``path_paths`` transition maps.

Pathway matrices are *contributions*, not normalized density operators: the
ket/bra masks generally make them non-Hermitian.  Use ``density_matrices`` or
``populations`` for physical state populations.  No relaxation or orientation
average is applied.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from McUtils.Data import UnitsData
from .NonlinearResponse import (
    TransitionData, experiment_defaults, expand_transition_data,
    interpret_polarization, prep_liouville_spaces,
    prep_nonlinear_transition_data,
)

__all__ = [
    "DensityMatrixPathHistory",
    "PathResolvedDensityMatrixResult",
    "propagate_density_matrix_paths",
]


@dataclass
class DensityMatrixPathHistory:
    """Matrix-valued contribution selected by one ``path_paths`` entry."""

    sign: int
    liouville_path: tuple[int, ...]
    density_matrices: np.ndarray  # (initial + pulses, state, state)
    allowed_elements: np.ndarray  # Boolean mask with the same shape


@dataclass
class PathResolvedDensityMatrixResult:
    """Complete physical history plus pathway-restricted matrix histories."""

    total_space: BasisStateSpace
    energies: np.ndarray
    transition_moments: np.ndarray
    density_matrices: np.ndarray
    path_histories: dict[str | None, list[DensityMatrixPathHistory]]

    @property
    def populations(self) -> np.ndarray:
        """Nonnegative state populations after each driving pulse."""
        return np.real(np.diagonal(self.density_matrices, axis1=-2, axis2=-1))


def _initial_density(total_space, initial_states, initial_density):
    size = len(total_space)
    if initial_density is not None:
        rho = np.asarray(initial_density, dtype=complex)
        if rho.shape != (size, size):
            raise ValueError(f"initial_density must have shape ({size}, {size})")
        if not np.allclose(rho, rho.conj().T, atol=1e-12):
            raise ValueError("initial_density must be Hermitian")
        if not np.isclose(np.trace(rho), 1., atol=1e-12):
            raise ValueError("initial_density must have unit trace")
        if np.linalg.eigvalsh(rho).min() < -1e-10:
            raise ValueError("initial_density must be positive semidefinite")
        return rho.copy()

    indices = np.asarray(total_space.find(initial_states), dtype=int).ravel()
    if not len(indices) or np.any(indices < 0):
        raise ValueError("initial_states must be present in total_space")
    rho = np.zeros((size, size), dtype=complex)
    rho[indices, indices] = 1 / len(indices)
    return rho


def _advance_path(index_set, liouville_path, pulse_count, size):
    """Build reachable matrix-element and transition-edge masks per stage."""
    masks = np.zeros((pulse_count + 1, size, size), dtype=bool)
    edges = np.zeros((pulse_count, size, size), dtype=bool)
    left, right = index_set[0]
    current = {(int(i), int(j)) for i, j in zip(left, right)}
    for i, j in current:
        masks[0, i, j] = True

    for stage in range(1, pulse_count + 1):
        # The final detection interaction has no fixed side in path_paths.
        sides = ((int(liouville_path[stage - 1]),)
                 if stage <= len(liouville_path) else (0, 1))
        new = set()
        for side in sides:
            transitions = index_set[stage][side]
            if not isinstance(transitions, dict):
                raise ValueError("path_paths contains no interaction map at this stage")
            for ket, bra in current:
                old = ket if side == 0 else bra
                for target in transitions.get(old, ()):
                    target = int(target)
                    pair = (target, bra) if side == 0 else (ket, target)
                    new.add(pair)
                    edges[stage - 1, old, target] = True
                    edges[stage - 1, target, old] = True
        current = new
        for i, j in current:
            masks[stage, i, j] = True
    return masks, edges


def _interaction_operator(field, transition_moments):
    """Contract a field vector/tensor with Psience's transition dipoles."""
    field = np.asarray(field, dtype=complex)
    if field.shape == (3,):
        operator = np.einsum("ijc,c->ij", transition_moments, field)
    elif field.shape == transition_moments.shape:
        operator = np.einsum("ijc,ijc->ij", transition_moments, field)
    else:
        raise ValueError(
            "interaction_tensor must return a 3-vector or an "
            f"array shaped {transition_moments.shape}"
        )
    if not np.all(np.isfinite(operator)) or not np.allclose(operator, operator.conj().T, atol=1e-12):
        raise ValueError("the dipole-contracted interaction must be Hermitian and finite")
    return operator


def propagate_density_matrix_paths(
        transition_data,
        *,
        initial_states=0,
        initial_density=None,
        experiment_type="2dir",
        paths=None,
        phases=None,
        selection_rules=None,
        num_interactions=None,
        num_pulses=3,
        pulse_areas=(1., 1., 1.),
        polarizations="XXX",
        delay_times=(0., 0.),
        interaction_tensor=None,
        orientational_averaging=None,
        state_filter_opts=None,
        prepared_spaces=None,
        frequency_unit="Wavenumbers",
        time_unit="PicoSeconds",
) -> PathResolvedDensityMatrixResult:
    """Propagate density matrices on the basis selected by Liouville paths.

    ``interaction_tensor(total_space, delay_times, pulse_index)`` may return
    either a Cartesian 3-vector of pulse field components or an
    ``(nstates, nstates, 3)`` field tensor.  Both are contracted with the
    transition dipoles from ``prep_nonlinear_transition_data`` to form the
    pulse generator.  Pulse indices are zero-based.  When no callable is
    supplied, ``pulse_areas[k] * polarizations[k]`` is used instead.

    The complete history uses the union of edges permitted by all paths at
    each step.  Each path history uses only its own transition edges and
    projects onto its reachable ket/bra matrix elements after the kick.
    Thus path histories can be complex, non-Hermitian, and non-normalized.
    They should not be interpreted as independent physical populations.

    ``prepared_spaces`` may supply the ``(total_space, path_paths)`` pair
    returned by ``prep_liouville_spaces`` to reuse an existing enumeration.
    """
    if orientational_averaging is not None:
        raise NotImplementedError("orientational averaging of pulse propagation is not implemented")
    if experiment_type not in experiment_defaults:
        raise ValueError(f"unknown experiment_type {experiment_type!r}")
    defaults = experiment_defaults[experiment_type]
    if paths is None:
        paths = defaults["paths"]
    if phases is None and selection_rules is None:
        phases = defaults.get("phases")
    if num_interactions is None:
        num_interactions = defaults["num_interactions"]
    if not 1 <= num_pulses <= num_interactions:
        raise ValueError("num_pulses must lie between 1 and num_interactions")
    delays = np.asarray(delay_times, dtype=float)
    if delays.shape != (num_pulses - 1,) or not np.all(np.isfinite(delays)) or np.any(delays < 0):
        raise ValueError(f"delay_times must contain {num_pulses - 1} finite nonnegative values")

    if isinstance(transition_data, TransitionData):
        td = transition_data
    else:
        td = prep_nonlinear_transition_data(transition_data)
    if not hasattr(td.states, "as_excitations"):
        raise ValueError("transition_data must contain state vectors")
    if not hasattr(initial_states, "as_excitations"):
        if isinstance(initial_states, (int, np.integer)):
            initial_states = BasisStateSpace(
                HarmonicOscillatorProductBasis(td.states.ndim),
                [int(initial_states)], mode="indices",
            )
        else:
            initial_states = BasisStateSpace(
                HarmonicOscillatorProductBasis(td.states.ndim), initial_states,
            )

    if prepared_spaces is None:
        total_space, path_paths = prep_liouville_spaces(
            initial_states, paths, num_interactions=num_interactions,
            phases=phases, selection_rules=selection_rules,
            **({} if state_filter_opts is None else state_filter_opts),
        )
    else:
        total_space, path_paths = prepared_spaces
    expanded = expand_transition_data(td, total_space)
    transition_moments = np.asarray(expanded.transition_moments, dtype=float)
    if not np.allclose(transition_moments, transition_moments.swapaxes(0, 1), atol=1e-12):
        raise ValueError("transition moments must be Hermitian")

    initial = _initial_density(total_space, initial_states, initial_density)
    reference = int(np.flatnonzero(np.real(np.diag(initial)) > 0)[0])
    energies = np.asarray(expanded.frequencies[reference], dtype=float)
    if not np.all(np.isfinite(energies)):
        raise ValueError("state energies must be finite")
    frequency_time_factor = (
        2 * np.pi * UnitsData.convert(frequency_unit, "Hertz")
        * UnitsData.convert(time_unit, "Seconds")
    )

    if interaction_tensor is None:
        areas = np.asarray(pulse_areas, dtype=float)
        pol = interpret_polarization(polarizations)
        if areas.shape != (num_pulses,) or pol.shape != (num_pulses, 3):
            raise ValueError("pulse_areas and polarizations must match num_pulses")
        norms = np.linalg.norm(pol, axis=1)
        if not np.all(np.isfinite(areas)) or not np.all(np.isfinite(pol)) or np.any(norms <= 0):
            raise ValueError("pulse areas and polarizations must be finite, with nonzero directions")
        fields = areas[:, None] * pol / norms[:, None]
        interaction_tensor = lambda _space, _delays, k: fields[k]
    elif not callable(interaction_tensor):
        raise TypeError("interaction_tensor must be callable")

    path_map = path_paths if isinstance(path_paths, dict) else {None: path_paths}
    masks_by_signal = {}
    edges_by_signal = {}
    all_edges = np.zeros((num_pulses, len(total_space), len(total_space)), dtype=bool)
    for signal, entries in path_map.items():
        masks_by_signal[signal] = []
        edges_by_signal[signal] = []
        for sign, path, index_set in entries:
            masks, edges = _advance_path(index_set, path, num_pulses, len(total_space))
            masks_by_signal[signal].append(masks)
            edges_by_signal[signal].append(edges)
            all_edges |= edges

    density = initial.copy()
    density_history = [density.copy()]
    path_density = {
        signal: [initial * masks[0] for masks in masks_by_signal[signal]]
        for signal in path_map
    }
    path_history = {
        signal: [[rho.copy()] for rho in path_density[signal]]
        for signal in path_map
    }
    for pulse_index in range(num_pulses):
        field = interaction_tensor(total_space, tuple(delays), pulse_index)
        generator = _interaction_operator(field, transition_moments)
        if pulse_index:
            phase = np.exp(-1j * frequency_time_factor * energies * delays[pulse_index - 1])
            density = phase[:, None] * density * phase.conj()[None, :]
        full_kick = generator * all_edges[pulse_index]
        unitary = expm(1j * full_kick)
        density = unitary @ density @ unitary.conj().T
        density_history.append(density.copy())

        for signal, entries in path_map.items():
            for path_index, _entry in enumerate(entries):
                rho = path_density[signal][path_index]
                if pulse_index:
                    rho = phase[:, None] * rho * phase.conj()[None, :]
                path_kick = generator * edges_by_signal[signal][path_index][pulse_index]
                path_unitary = expm(1j * path_kick)
                rho = path_unitary @ rho @ path_unitary.conj().T
                rho = np.where(masks_by_signal[signal][path_index][pulse_index + 1], rho, 0)
                path_density[signal][path_index] = rho
                path_history[signal][path_index].append(rho.copy())

    histories = {
        signal: [
            DensityMatrixPathHistory(
                sign=int(sign),
                liouville_path=tuple(int(x) for x in path),
                density_matrices=np.asarray(path_history[signal][index]),
                allowed_elements=masks_by_signal[signal][index],
            )
            for index, (sign, path, _index_set) in enumerate(entries)
        ]
        for signal, entries in path_map.items()
    }
    return PathResolvedDensityMatrixResult(
        total_space=total_space,
        energies=energies,
        transition_moments=transition_moments,
        density_matrices=np.asarray(density_history),
        path_histories=histories,
    )
