import numpy as np

import McUtils.Numputils as nput
import McUtils.Combinatorics as comb

def nested_commutator_expansion(k, side='left'):
    if side == 'left':
        comm = 0
        for i in range(1, k):
            comm = [comm, i]
    else:
        comm = k - 1
        for i in range(k - 2, -1, -1):
            comm = [i, comm]
    return comm

def closed_system_propagator(interaction_generator, time_delays):
    iteractions = [interaction_generator(t) for t in time_delays]
    interaction_n = len(iteractions)
    comms = nput.commutator_terms(nested_commutator_expansion(interaction_n))
    #TODO: allow for pathway filtering on the commutator terms
    return nput.commutator_evaluate(
        comms, iteractions,
        direct=False,
        normalized=True
    )

def liouville_pathways(k):
    final_states = [[i, k-i] for i in range(k+1)]
    return np.concatenate([ # just lattice paths, but people don't call them that
        comb.UniquePermutations([0]*l + [1]*r)
        for l,r in final_states
    ])

def evaluate_pathway(pathway, initial_state, transition_graph):
    ...

def liouville_pathway_propagator(pathways):
    signs = (-1)**np.sum(pathways, axis=1)
    signal = ...
    for s,p in pathways:
        ...

def kronecker_delta_product_evaluate():
    ...

def isotropic_tensor_basis(k):
    ...

def rotational_averaging_data(k) -> (dict, np.ndarray, dict):
    ...

def rotational_average(isotropic_basis):
    ...