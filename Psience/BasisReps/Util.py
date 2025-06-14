
import numpy as np
import McUtils.Coordinerds as coordops

__all__ = [
    "StateMaker",
    "modify_hamiltonian",
    "modify_internal_hamiltonian"
]

class StateMaker:
    """
    A tiny but useful class to make states based on their quanta
    of excitation
    """

    def __init__(self, ndim, mode='low-high'):
        self.ndim = ndim
        self.mode = mode

    @classmethod
    def parse_state(cls, state, mode='low-high'):
        nzp = np.nonzero(state)
        if len(nzp) > 0: nzp = nzp[0]
        if len(nzp) == 0: return "()"
        if mode == 'low-high':
            pos = len(state) - nzp
            state_quants = reversed(list(zip(pos, nzp)))
        elif mode == 'high-low':
            pos = nzp + 1
            state_quants = list(zip(pos, nzp))
        else:
            pos = nzp
            state_quants = list(zip(pos, nzp))
        return "".join("{}({})".format(p, state[i]) for p,i in state_quants)

    def make_state(self, *specs, mode=None):

        if mode is None:
            mode = self.mode

        state = [0] * self.ndim
        for s in specs:
            if isinstance(s, (int, np.integer)):
                i = s
                q = 1
            else:
                i,q = s
            if mode == 'low-high':
                state[-i] = q
            elif mode == 'high-low':
                state[i-1] = q
            elif mode == 'normal':
                state[i] = q
            else:
                raise ValueError("don't know what to do with filling mode '{}'".format(mode))
        return state

    def __call__(self, *specs, mode=None):
        return self.make_state(*specs, mode=mode)


def modify_hamiltonian(ham,
                       index_specifiers,
                       index_key_function=None,
                       coupling_key_function=None,
                       scaling_types=None,
                       coupling_types=None
                       ):

    ham = np.array(ham)
    if scaling_types is None:
        scaling_types = {}
    if coupling_types is None:
        coupling_types = {}
    if index_key_function is None:
        index_key_function = lambda i,t:(t1, (t1,))
    if coupling_key_function is None:
        coupling_key_function = lambda ij,t1,t2:((t1,t2), (t2, t1))

    for i,t1 in enumerate(index_specifiers):
        keys = index_key_function(i,t1)
        for t in keys:
            if t in coupling_types:
                ham[i, i] = coupling_types[t]
                break
            elif t in scaling_types:
                ham[i, i] *= scaling_types[t]
                break

        for j,t2 in enumerate(index_specifiers[i+1:]):
            k = j + i + 1
            keys = coupling_key_function((i,k), t1, t2)
            for t in keys:
                if t in coupling_types:
                    ham[i, k] = coupling_types[t]
                    ham[k, i] = coupling_types[t]
                    break
                elif t in scaling_types:
                    ham[i, k] *= scaling_types[t]
                    ham[k, i] *= scaling_types[t]
                    break

    return ham

def _internal_index_key_function(internals, tags, atoms):
    def index_key_function(i, internal):
        tag = tags[i]
        if tag is None and atoms is not None:
            tag = "".join(atoms[j] for j in internal)
        return (internal, (internal,),
                (internal, internal),
                tag,
                (tag,),
                (len(internal), tag, tag)
                )
    return index_key_function

def _internal_coupling_key_function(internals, tags, atoms):
    def coupling_key_function(ij, internal_1, internal_2):
        i, j = ij
        tag1 = tags[i]
        tag2 = tags[j]
        print(tag1, tag2, internal_1, internal_2)
        if tag1 is None and atoms is not None:
            tag1 = "".join(atoms[k] for k in internal_1)
        if tag2 is None and atoms is not None:
            tag2 = "".join(atoms[k] for k in internal_2)
        num_shared = sum(
            1 if i1 in internal_2 else 0
            for i1 in internal_1
        )

        if num_shared > 0:
            keys = (
                (internal_1, internal_2), (internal_2, internal_1),
                (num_shared, tag1, tag2),
                (num_shared, tag2, tag1),
                (tag1, tag2),
                (tag2, tag1),
            )
        else:
            keys = ((internal_1, internal_2), (internal_2, internal_1))
        return keys
    return coupling_key_function

def modify_internal_hamiltonian(ham,
                                internals,
                                atoms=None,
                                tags=None,
                                index_key_function=None,
                                coupling_key_function=None,
                                scaling_types=None,
                                coupling_types=None
                                ):
    if hasattr(internals, 'items'):
        internals = {
            coordops.canonicalize_internal(k):t
            for k,t in internals.items()
        }
        if tags is None:
            tags = list(internals.values())
        internals = list(internals.keys())
    else:
        if tags is None:
            tags = [
                list(i.values())[0]
                    if hasattr(i, 'values') else
                None
                for i in internals
            ]
        internals = [
            i
                if hasattr(i, 'values') else
            coordops.canonicalize_internal(i)
            for i in internals
        ]

    if index_key_function is None:
        index_key_function = _internal_index_key_function(internals, tags, atoms)
    if coupling_key_function is None:
        coupling_key_function = _internal_coupling_key_function(internals, tags, atoms)

    return modify_hamiltonian(ham, internals,
                              index_key_function=index_key_function,
                              coupling_key_function=coupling_key_function,
                              scaling_types=scaling_types,
                              coupling_types=coupling_types
                              )


