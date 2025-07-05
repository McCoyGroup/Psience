import numpy as np
import scipy.linalg

import McUtils.Coordinerds as coordops
from .StateSpaces import AbstractStateSpace, BasisStateSpace, BraKetSpace
# from .Representations import ExpansionRepresentation
from .HarmonicOscillator import HarmonicOscillatorProductBasis
from . import Util as util

__all__ = [
    "LocalHarmonicModel"
]

class LocalHarmonicModel:
    def __init__(self, f, g,
                 internals=None,
                 anharmonic_scalings=None,
                 anharmonic_constants=None,
                 freqs=None,
                 **operator_settings):
        self.f = np.asanyarray(f)
        self.g = np.asanyarray(g)
        if freqs is None:
            freqs = scipy.linalg.eigvalsh(self.f, self.g, type=3)
            freqs = np.sign(freqs) * np.sqrt(np.abs(freqs))
        self.freqs = np.asanyarray(freqs)
        self.basis = HarmonicOscillatorProductBasis(len(f))
        self.rep = 1 / 2 * (
                self.basis.representation(
                    "p", "p",
                    coeffs=self.g,
                    name="T",
                    **operator_settings
                ) +
                self.basis.representation(
                    "x", "x",
                    coeffs=self.f,
                    axes=[[0, 1], [0, 1]],
                    name="V",
                    **operator_settings
                )
        )
        self.internals = internals
        self.scalings = anharmonic_scalings
        self.constants = anharmonic_constants

    def get_hamiltonian(self, states, coupled_space=None, remove_zpe=True, **rep_opts):
        if not isinstance(states, AbstractStateSpace):
            states = BasisStateSpace(self.basis, states)
        if coupled_space is None:
            if isinstance(states, BraKetSpace):
                coupled_space = states
                states = states.bras.take_unique()
            else:
                coupled_space = states.get_representation_brakets()
        ham = self.rep.get_representation_matrix(
            coupled_space,
            states,
            zero_element_warning=False,
            **rep_opts
        ).asarray()

        if remove_zpe:
            zpe = np.sum(np.abs(self.freqs)) * 1/2
            ham = ham - np.eye(len(ham)) * zpe

        if self.scalings is not None or self.constants is not None:
            if self.internals is not None:
                ham = util.modify_internal_hamiltonian(
                    ham,
                    self.internals,
                    states=states.excitations,
                    scaling_types=self.scalings,
                    coupling_types=self.constants,
                )
            else:
                ham = util.modify_hamiltonian(
                    ham,
                    [tuple(s) for s in states.excitations],
                    scaling_types=self.scalings,
                    coupling_types=self.constants,
                )

        return ham

    @classmethod
    def from_modes(cls, nms, internals, oblique=True, include_complement=False, **opts):
        from ..Modes import NormalModes
        nms:NormalModes

        nms = nms.remove_frequency_scaling().remove_mass_weighting()
        loc_modes = nms.localize(internals=internals)
        if oblique:
            loc_modes = loc_modes.make_oblique()
        if include_complement:
            loc_modes = loc_modes.get_complement(concatenate=True)

        loc_modes = loc_modes.make_mass_weighted()
        f = loc_modes.compute_hessian()
        g = loc_modes.compute_gmatrix()

        return cls(
            f, g,
            internals=internals,
            **opts
        )

    @classmethod
    def from_molecule(cls, mol, internals=None,
                      coordinate_filter=None,
                      include_stretches=True,
                      include_bends=True,
                      include_dihedrals=False,
                      **opts):
        from ..Molecools import Molecule
        mol: Molecule

        if internals is None:
            st,bo,di = coordops.get_stretch_coordinate_system(
                [tuple(b[:2]) for b in mol.bonds]
            )
            bits = []
            if include_stretches:
                bits.append(st)
            if include_bends:
                bits.append(bo)
            if include_dihedrals:
                bits.append(di)
            internals = bits[0]
            for b in bits[1:]:
                internals = internals + b

            labels = mol.edge_graph.get_label_types()
            internals = {
                c:coordops.get_coordinate_label(
                    c,
                    labels
                )
                for c in internals
            }

            if coordinate_filter is not None:
                internals = coordinate_filter(internals)

        return cls.from_modes(
            mol.get_normal_modes(),
            internals=internals,
            **opts
        )