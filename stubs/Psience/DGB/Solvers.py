import numpy as np, scipy as sp, collections

class DGBEigensolver:

    @classmethod
    def get_orthogonal_transform(self, S, min_singular_value=None, subspace_size=None):
        ...

    @classmethod
    def classic_eigensolver(self, H, S, hamiltonian, min_singular_value=None, subspace_size=None, nodeless_ground_state=False):
        ...
    eigensimilarity_cutoff = None
    eigensimilarity_chunk_size = 3
    similar_determinant_cutoff = 0.05
    similar_det_diff_cutoff = 0.05
    similarity_cutoff_rounding = 0.05

    @classmethod
    def get_eigensimilarity_subspace_size(cls, H, S, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, similar_det_diff_cutoff=None, similarity_cutoff_rounding=None):
        ...

    @classmethod
    def similarity_mapped_solver(self, H, S, hamiltonian, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None):
        ...
    default_shift = 2

    @classmethod
    def shift_similarity_solver(self, H, S, hamiltonian, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, similarity_shift=None, similar_det_diff_cutoff=None, similarity_cutoff_rounding=None):
        ...

    @classmethod
    def low_rank_solver(cls, H, S, hamiltonian, low_rank_energy_cutoff=None, low_rank_overlap_cutoff=None, low_rank_shift=None):
        ...

    @classmethod
    def cholesky_solver(cls, H, S, hamiltonian):
        ...

    @classmethod
    def fix_heiberger(self, H, S, hamiltonian, eps=1e-05):
        ...