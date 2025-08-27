import numpy as np
import McUtils.Numputils as nput
from .Elements import *
from .PointGroups import PointGroup

__all__ = [
    "apply_symmetries",
    "symmetrize_coordinates"
]

def _symmetry_reduce(coords, op:np.ndarray, labels=None):
    perm = nput.symmetry_permutation(coords, op)
    cycles = nput.permutation_cycles(perm, return_groups=True)
    coords = np.array([
        coords[p[0]]
        for p in cycles
    ])
    if labels is not None:
        labels = [
            labels[p[0]]
            for p in cycles
        ]
        return coords, labels
    else:
        return coords

def prep_symmetry_operations(symmetry_elements: 'PointGroup|list[SymmetryElement|np.ndarray]'):
    if isinstance(symmetry_elements, PointGroup):
        symmetry_elements = symmetry_elements.elements
    return [
        e.get_transformation()
            if isinstance(e, SymmetryElement) else
        np.asanyarray(e)
        for e in symmetry_elements
    ]

def apply_symmetries(coords, symmetry_elements: 'PointGroup|list[SymmetryElement|np.ndarray]', labels=None, tol=1e-1):
    coords = np.asanyarray(coords)
    symmetry_elements = prep_symmetry_operations(symmetry_elements)

    for e in symmetry_elements:
        new_coords = coords @ e
        coord_diffs = np.linalg.norm(coords[:, np.newaxis, :] - new_coords[np.newaxis, :, :], axis=-1)
        dupe_pos = np.where(coord_diffs < tol)
        new_labs = None
        if len(dupe_pos) > 0 and len(dupe_pos[0]) > 0:
            rem = np.setdiff1d(np.arange(len(new_coords)), dupe_pos[0])
            if labels is not None:
                new_labs = [labels[l] for l in rem]
            new_coords = new_coords[rem,]
        elif labels is not None:
            new_labs = labels
        if labels is not None:
            labels = labels + new_labs
        coords = np.concatenate([coords, new_coords], axis=0)

    if labels is not None:
        return coords, labels
    else:
        return coords

def symmetrize_coordinates(coords,
                           symmetry_elements: 'PointGroup|list[SymmetryElement|np.ndarray]',
                           labels=None,
                           masses=None,
                           groups=None,
                           tol=1e-1,
                           mass_tol=1,
                           expand=True
                           ):
    coords = np.asanyarray(coords)
    if groups is None:
        if labels is not None:
            label_map = {}
            g_vec = []
            for l in labels:
                label_map[l] = label_map.get(l, len(label_map))
                g_vec.append(l)
            (_, groups), _ = nput.group_by(np.arange(len(g_vec)), g_vec)
        elif masses is not None:
            (_, groups), _ = nput.group_by(np.arange(len(masses)), np.round(masses/mass_tol))
        else:
            groups = [np.arange(len(coords))]
    symmetry_elements = prep_symmetry_operations(symmetry_elements)

    all_coords = []
    all_labels = [] if labels is not None else None
    for g in groups:
        subcoords = coords[g,]
        if labels is not None:
            sublabels = [labels[gg] for gg in g]
        else:
            sublabels = None
        for e in symmetry_elements:
            subcoords = _symmetry_reduce(subcoords, e, labels=sublabels)
            if sublabels is not None:
                subcoords, sublabels = subcoords
        if expand:
            subcoords = apply_symmetries(subcoords, symmetry_elements, labels=sublabels, tol=tol)
            if sublabels is not None:
                subcoords, sublabels = subcoords
        if all_labels is not None:
            all_labels.append(sublabels)
        all_coords.append(subcoords)

    all_coords = np.concatenate(all_coords)
    all_labels = sum(all_labels, [])
    if labels is not None:
        return all_coords, all_labels
    else:
        return all_coords

def get_mode_symmetries():
    ...