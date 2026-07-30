"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""
import numpy as np
from McUtils.Zachary import Mesh
import McUtils.Plots as plt
from Psience.Wavefun import Wavefunction, Wavefunctions
from .Gaussians import DGBGaussians
from .Coordinates import DGBCartesians, DGBWatsonModes
__all__ = ['DGBWavefunctions', 'DGBWavefunction']
__reload_hook__ = ['.Components']

class DGBWavefunction(Wavefunction):

    def __init__(self, energy, data, gaussians: DGBGaussians=None, **opts):
        ...

    def get_dimension(self):
        ...

    def plot(self, figure=None, domain=None, plot_centers=False, return_values=False, **opts):
        ...

    def to_cartesian_wavefunction(self):
        """
        Projects the wavefunction back to Cartesians
        :return:
        """
        ...

    def plot_cartesians(self, xyz_sel=None, *, atom_sel=None, figure=None, plot_centers=False, atom_styles=None, adjust_levels=False, projection_plot_density_cutoff=None, return_values=False, **plot_styles):
        ...

    def evaluate(self, points):
        ...

    def marginalize_out(self, dofs, rescale=True, return_scaling=False) -> 'DGBWavefunction':
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom, returning a projected wave function object

        :return:
        :rtype: Wavefunction
        """
        ...

class DGBWavefunctions(Wavefunctions):
    wavefunction_class = DGBWavefunction

    def __init__(self, energies=None, wavefunctions=None, hamiltonian=None, **opts):
        ...

    def as_cartesian_wavefunction(self):
        """
        Projects the wavefunction back to Cartesians
        :return:
        """
        ...

    @property
    def hamiltonian(self):
        ...

    def operator_representation(self, op, embed=True, expansion_degree=None, quadrature_degree=None, expansion_type=None):
        ...

    def expectation(self, op, expansion_degree=None, quadrature_degree=None, expansion_type=None, embed=True, other=None):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        ...

    def localize(self, criterion, which=None):
        """
        Find a transformation that maximally localizes the wavefunctions in the Boys' sense
        by minimizing <r^2> - <r>^2 over unitary transformations

        :param criterion:
        :param which:
        :return:
        """
        ...

    def __repr__(self):
        ...