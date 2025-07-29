import itertools

import numpy as np
from collections import namedtuple
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
import McUtils.Numputils as nput

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

def _get_tag_match(tag, mod_types):
    return mod_types.get(tag)

def modify_hamiltonian(ham,
                       index_specifiers,
                       index_key_function=None,
                       coupling_key_function=None,
                       scaling_types=None,
                       coupling_types=None,
                       shift_types=None,
                       tag_match=None
                       ):

    ham = np.array(ham)
    if scaling_types is None:
        scaling_types = {}
    if coupling_types is None:
        coupling_types = {}
    if shift_types is None:
        shift_types = {}
    if index_key_function is None:
        index_key_function = lambda i,t:(t1, (t1,))
    if coupling_key_function is None:
        coupling_key_function = lambda ij,t1,t2:((t1,t2), (t2, t1))

    if tag_match is None:
        tag_match = _get_tag_match

    for i,t1 in enumerate(index_specifiers):
        keys = index_key_function(i,t1)
        for t in keys:
            s = tag_match(t, coupling_types)
            if s is not None:
                ham[i, i] = s
                # break
            s = tag_match(t, scaling_types)
            if s is not None:
                ham[i, i] *= s
                # break
            s = tag_match(t, shift_types)
            if s is not None:
                ham[i, i] += s
                # break

        for j,t2 in enumerate(index_specifiers[i+1:]):
            k = j + i + 1
            keys = coupling_key_function((i,k), t1, t2)
            for t in keys:
                s = tag_match(t, coupling_types)
                if s is not None:
                    ham[i, k] = s
                    ham[k, i] = s
                    # break
                s = tag_match(t, scaling_types)
                if s is not None:
                    ham[i, k] *= s
                    ham[k, i] *= s
                    # break
                s = tag_match(t, shift_types)
                if s is not None:
                    ham[i, k] += s
                    ham[k, i] += s
                    # break

    return ham


product_state = namedtuple("product_state", ["indices", "quanta", "tags"])
def _canonicalize_internal_state_key(k):
    if isinstance(k, str):
        return product_state(None, None, (k,))
    elif all(isinstance(s, str) for s in k):
        return product_state(None, None, (tuple(k),))
    elif all(nput.is_numeric(x) for x in k):
        return product_state(
            tuple(n for n,s in enumerate(k) if s != 0),
            tuple(s for n,s in enumerate(k) if s != 0),
            None
        )
    else:
        if nput.is_numeric(k[0]):
            k = [k]
        return product_state(
            None,
            tuple(n for n,s in k),
            tuple(s for n,s in k)
        )

def _tag_match(tags1, tags2):
    if tags2 is None or tags1 is None: return False
    ml = min([len(tags1), len(tags2)])
    return tags1[-ml:] == tags2[-ml:]


def _product_state_match(state1:product_state, state2:product_state):
    if (
            state1.indices is not None
            and state2.indices is not None
    ):
        if state1.indices != state2.indices: return False
        if state1.quanta is None: return True
        if state2.quanta is None: return True
        return state1.quanta == state2.quanta
    else:
        if (
                state2.tags is None
                or state1.tags is None
        ):
            return (
                    state1.quanta is None
                    or state2.quanta is None
                    or list(sorted(state1.quanta)) == list(sorted(state2.quanta))
            )
        else:
            if len(state1.tags) != len(state2.tags): return False
            if state1.quanta is None or state2.quanta is None:
                # any match is valid of these can match
                for p in itertools.permutations(range(len(state1.tags)), len(state1.tags)):
                    if all(
                        _tag_match(s1, s2)
                        for s1, s2 in zip(
                            [state1.tags[i] for i in p],
                            state2.tags
                        )
                    ):
                        return p
            else:
                for p in itertools.permutations(range(len(state1.tags)), len(state1.tags)):
                    if all(
                            (
                                    n1 is None or n2 is None
                                    or n1 == n2
                            ) and _tag_match(s1, s2)
                            for n1,s1, n2,s2 in zip(
                                [state1.quanta[i] for i in p],
                                [state1.tags[i] for i in p],
                                state2.quanta,
                                state2.tags
                            )
                    ):
                        return p

def _canonicalize_scaling_key(key):
    if isinstance(key, product_state):
        return key
    elif (
            len(key) == 3
            and isinstance(key[1], product_state)
            and isinstance(key[2], product_state)
    ):
        shats = key[0]
        if nput.is_numeric(shats):
            shats = [[shats]]
        else:
            shats = [[s] if nput.is_numeric(s) else s for s in shats]
        shats = tuple(tuple(int(i) for i in s) for s in shats)
        return (
            shats,
            key[1],
            key[2]
        )
    elif (
            len(key) == 2
            and isinstance(key[0], product_state)
            and isinstance(key[1], product_state)
    ):
        return (
            key[0],
            key[1]
        )
    if (
            isinstance(key, str)
            or all(isinstance(k, str) for k in key)
            or all(nput.is_numeric(x) for x in key)
    ):
        return _canonicalize_internal_state_key(key)
    if len(key) == 3:
        shats = key[0]
        if nput.is_numeric(shats):
            shats = [[shats]]
        else:
            shats = [[s] if nput.is_numeric(s) else s for s in shats]
        shats = tuple(tuple(int(i) for i in s) for s in shats)
        return (
            shats,
            _canonicalize_internal_state_key(key[1]),
            _canonicalize_internal_state_key(key[2])
        )
    else:
        return (
            _canonicalize_internal_state_key(key[0]),
            _canonicalize_internal_state_key(key[1])
        )

def _find_internal_coupling_match(spec, couplings):
    base = couplings.get(spec)
    if base is not None: return base
    if nput.is_numeric(spec) or not (
            isinstance(spec, product_state)
            or any(isinstance(s, product_state) for s in spec)
    ): return None


    if isinstance(spec, product_state):
        for s,v in couplings.items():
            if isinstance(s, product_state) and _product_state_match(s, spec):
                return v
    else:
        if not isinstance(spec[0], product_state):
            shared_ats = spec[0]
            spec = spec[1:]
        else:
            shared_ats = None

        spec1, spec2 = spec
        for s,v in couplings.items():
            if not isinstance(s, product_state):
                if not isinstance(s[0], product_state):
                    shared_ats_test = s[0]
                    s = s[1:]
                else:
                    shared_ats_test = None

                good_test = (
                    shared_ats_test is None and
                    shared_ats is None
                ) or (
                        shared_ats_test is not None and
                        shared_ats is not None and
                        list(sorted(len(t) for t in shared_ats_test))
                            == list(sorted(len(t) for t in shared_ats))
                )
                if not good_test: continue

                s1,s2 = s
                perm_1 = _product_state_match(s1, spec1)
                perm_2 = _product_state_match(s2, spec2)
                if perm_1 is None or perm_2 is None:
                    perm_1 = _product_state_match(s1, spec2)
                    perm_2 = _product_state_match(s2, spec1)

                if not perm_1: continue
                if not perm_2: continue

                if perm_1 is not None and perm_2 is not None:
                    if shared_ats is not None:
                        if perm_1 is True or perm_2 is True:
                            return list(sorted(shared_ats)) == list(sorted(shared_ats_test))
                        elif all(
                                shared_ats[p1] == shared_ats_test[p2]
                                for p1,p2 in zip(perm_1, perm_2)
                        ):
                            return v
                    else:
                        return v

def _internal_index_key_function(internals, tags, atoms):
    def index_key_function(i, state):
        return (i,
                (state,),
                (state, state),
                product_state(
                    tuple(n for n, s in enumerate(state) if s != 0),
                    tuple(s for n, s in enumerate(state) if s != 0),
                    tuple(tags[n] for n, s in enumerate(state) if s != 0)
                )
                )
    return index_key_function

def _internal_coupling_key_function(internals, tags, atoms):
    def coupling_key_function(ij, state1, state2):
        i, j = ij
        matches = []
        for n1,s1 in enumerate(state1):
            if s1 == 0: continue
            sublist = []
            for n2,s2 in enumerate(state2):
                if s2 == 0: continue
                internal_1 = internals[n1]
                internal_2 = internals[n2]
                if not (nput.is_atomic(internal_1) or  nput.is_atomic(internal_2)):
                    num_shared = sum(
                        1 if i1 in internal_2 else 0
                        for i1 in internal_1
                    )
                else:
                    num_shared = 0
                sublist.append(num_shared)
            matches.append(tuple(sublist))

        return (
            (i, j),
            (state1, state2),
            (state2, state1),

            (
                tuple(matches),
                product_state(
                    tuple(n for n, s in enumerate(state1) if s != 0),
                    tuple(s for n, s in enumerate(state1) if s != 0),
                    tuple(tags[n] for n, s in enumerate(state1) if s != 0)
                ),
                product_state(
                    tuple(n for n, s in enumerate(state2) if s != 0),
                    tuple(s for n, s in enumerate(state2) if s != 0),
                    tuple(tags[n] for n, s in enumerate(state2) if s != 0)
                )
            ),
            (
                product_state(
                    tuple(n for n, s in enumerate(state1) if s != 0),
                    tuple(s for n, s in enumerate(state1) if s != 0),
                    tuple(tags[n] for n, s in enumerate(state1) if s != 0)
                ),
                product_state(
                    tuple(n for n, s in enumerate(state2) if s != 0),
                    tuple(s for n, s in enumerate(state2) if s != 0),
                    tuple(tags[n] for n, s in enumerate(state2) if s != 0)
                )
            )

        )
    return coupling_key_function


def modify_internal_hamiltonian(ham,
                                internals,
                                states=None,
                                atoms=None,
                                tags=None,
                                index_key_function=None,
                                coupling_key_function=None,
                                scaling_types=None,
                                coupling_types=None,
                                shift_types=None,
                                tag_match=None
                                ):
    if hasattr(internals, 'items'):
        internals = {
            (
                coordops.canonicalize_internal(k)
                    if not dev.is_atomic(k) and all(nput.is_int(cc) for cc in k) else
                k
            ):t
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
                if hasattr(i, 'values') or dev.is_atomic(i) or not all(nput.is_int(cc) for cc in i) else
            coordops.canonicalize_internal(i)
            for i in internals
        ]

    if states is None:
        states = np.eye(ham.shape[0], dtype=int)

    if index_key_function is None:
        index_key_function = _internal_index_key_function(internals, tags, atoms)
    if coupling_key_function is None:
        coupling_key_function = _internal_coupling_key_function(internals, tags, atoms)

    if scaling_types is not None:
        scaling_types = {
            _canonicalize_scaling_key(k):v
            for k,v in scaling_types.items()
        }
    if coupling_types is not None:
        coupling_types = {
            _canonicalize_scaling_key(k):v
            for k,v in coupling_types.items()
        }
    if shift_types is not None:
        shift_types = {
            _canonicalize_scaling_key(k):v
            for k,v in shift_types.items()
        }

    if tag_match is None:
        tag_match = _find_internal_coupling_match

    return modify_hamiltonian(ham,
                              [tuple(np.asanyarray(s).astype(int)) for s in states],
                              index_key_function=index_key_function,
                              coupling_key_function=coupling_key_function,
                              scaling_types=scaling_types,
                              coupling_types=coupling_types,
                              shift_types=shift_types,
                              tag_match=tag_match
                              )
