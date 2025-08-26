
import abc
import numpy as np
import enum
import itertools
import McUtils.Numputils as nput
import McUtils.Combinatorics as comb

from .Elements import *
from .PointGroups import *
from .Rotors import RotorTypes, identify_rotor_type

__all__ = [
    "PointGroupIdentifier"
]


class SymmetryEquivalentAtomData:
    __slots__ = ['coords', 'moms', 'com', 'axes', 'rotor_type', 'planar']
    coords: np.ndarray
    def __init__(self, coords, moms, com, axes, rotor_type, planar):
        self.coords = coords
        self.moms = moms
        self.com = com
        self.axes = axes
        self.rotor_type = rotor_type
        self.planar = planar

class PointGroupIdentifier:
    def __init__(self, coords, masses=None, groups=None,
                 tol=1e-8, mass_tol=1, mom_tol=1,
                 grouping_tol=1e-2,
                 verbose=False
                 ):
        self.mass_tol = mass_tol
        self.mom_tol = mom_tol
        self.tol = tol
        self.grouping_tol = grouping_tol
        self.verbose = verbose
        if masses is None:
            masses = np.ones(len(coords))
        masses = np.asanyarray(masses)
        if groups is None:
            groups = self.get_groups_from_masses(masses)
        self.groups = self.get_groups(coords, groups)
        self.group_orders = self.get_group_orders()
        self.coord_data = self.prep_coords(coords, masses=masses)
        self.group_coords = [
            self.prep_coords(self.coord_data.coords[g,], masses=masses[g,])
            for g in self.groups
        ]

    def get_groups(self, coords, base_groups):
        full_groups = []
        for g in base_groups:
            if len(g) < 2: continue
            subcoords = coords[g,]
            dists = nput.distance_matrix(subcoords)
            dists = np.sort(dists, axis=1)

            subgroups = []
            for i,d in enumerate(dists):
                new_sg = False
                for subgroup in subgroups:
                    if i in subgroup: break
                else:
                    new_sg = True
                    subgroup = {i}
                for j,d2 in enumerate(dists[i+1:]):
                    if np.max(np.abs(d - d2)) < self.grouping_tol:
                        subgroup.add(i+j+1)
                if new_sg and len(subgroup) > 1:
                    subgroups.append(subgroup)

            merged_groups = []
            for n,gg in enumerate(subgroups):
                for g2 in merged_groups:
                    if len(gg.intersection(g2)) > 0:
                        g2.update(gg)
                        break
                else:
                    merged_groups.append(gg)
            # print(dists)
            # dist_dm = nput.distance_matrix(dists, return_triu=True)
            # equiv_pos = np.where(dist_dm < self.tol)
            # if len(equiv_pos) == 0 or len(equiv_pos[0]) == 0:
            #     continue
            # else:
            #     subgroups = []
            #     row, col = np.triu_indices(len(dists), k=1)
            #     for x in equiv_pos[0]:
            #         e = row[x]
            #         p = col[x]
            #         for gg in subgroups:
            #             if e in gg:
            #                 gg.add(p)
            #                 break
            #             elif p in gg:
            #                 gg.add(e)
            #                 break
            #         else:
            #             subgroups.append({e,p})
            full_groups.extend(
                list(sorted(g[ss] for ss in s))
                for s in merged_groups
            )
        if self.verbose:
            print("Using groups:", full_groups)
        return full_groups

    def get_rotor_type(self, moms):
        return identify_rotor_type(moms, self.mom_tol)

    def get_group_orders(self):
        group_factors = None
        for g in self.groups:
            if len(g) < 2: continue
            if len(g) == 2:
                if group_factors is None:
                    group_factors = {2}
                continue
            subfactors = set()
            primes, ords = comb.prime_factorize(len(g))
            ords = np.array(ords)
            ord_pos = np.where(ords > 0)
            primes = primes[ord_pos]
            for exp in itertools.islice(
                    itertools.product(*[
                        range(p+1)
                        for p in ords[ord_pos]
                    ]),
                    1, None # skip 000
            ):
                exp = np.power(primes, np.array(exp, dtype=int))
                val = int(np.prod(exp, dtype=int))
                subfactors.add(val)
            if group_factors is None:
                group_factors = subfactors
            else:
                if len(group_factors) == 1 and list(group_factors)[0] == 2:
                    group_factors = group_factors.union(subfactors)
                else:
                    group_factors = group_factors.intersection(subfactors)
        return list(sorted(group_factors, key=lambda n:-n))

    def check_element(self, elem:SymmetryElement, verbose=False):
        tf = elem.get_transformation()
        for g in self.groups:
            coords = self.coord_data.coords[g,]
            new_coords = coords @ tf.T
            dm = np.linalg.norm(coords[:, np.newaxis] - new_coords[np.newaxis, :], axis=-1)
            # if verbose:
            #     print(dm)
            min_dists = np.argmin(dm, axis=1)
            max_min_dist = np.max(dm[np.arange(len(coords)), min_dists])
            if verbose:
                print("Max min-dist:", max_min_dist)
                print("Tol.:", len(coords) * self.tol)
            if max_min_dist > len(coords) * self.tol:
                return False
            diff_comp = coords - new_coords[min_dists,]
            max_comp = np.max(np.abs(diff_comp))
            if verbose:
                print("Max abs. displacement:", max_min_dist)
                print("Tol.:", self.tol)
            if max_comp > self.tol:
                return False
        return True

    def get_groups_from_masses(self, masses):
        (_, groups), _ = nput.group_by(np.arange(len(masses)), np.round(masses, self.mass_tol))
        return groups

    def prep_coords(self, coords, masses=None):
        coords = np.asanyarray(coords)
        if len(coords) == 1:
            return np.zeros(3), np.array([1, 0, 0]), np.eye(3)
        com = nput.center_of_mass(coords, masses)
        coords = coords - com[np.newaxis]
        inerts = nput.inertia_tensors(coords, masses)
        moms, axes = np.linalg.eigh(inerts)
        coords = coords @ axes
        rotor_type, planar = self.get_rotor_type(moms)
        return SymmetryEquivalentAtomData(coords, moms, com, axes, rotor_type, planar)

    def rotation_axis_iterator(self):
        for e in np.eye(3):
            yield e

        for g in self.groups:
            for n, i in enumerate(g):
                ax = self.coord_data.coords[i]
                if np.linalg.norm(ax) > 1e-6:
                    yield ax
                for j in g[n + 1:]:
                    ax = (self.coord_data.coords[i] + self.coord_data.coords[j]) / 2
                    if np.linalg.norm(ax) > 1e-6:
                        yield ax

    def reflection_plane_iterator(self):
        for e in np.eye(3):
            yield e
        for g in self.groups:
            for i,j in itertools.combinations(g, 2):
                ax = self.coord_data.coords[i] - self.coord_data.coords[j]
                yield ax

    def identify_point_group(self, realign=True):
        elements = []
        elem = InversionElement()
        has_inversion = self.check_element(elem)
        if has_inversion:
            elements.append(elem)

        if self.verbose:
            print("Rotor type:", self.coord_data.rotor_type)
        if self.coord_data.rotor_type == RotorTypes.Linear:
            if has_inversion:
                pg = PointGroup.from_name("D", -1)
            else:
                pg = PointGroup.from_name("C", -1)
        elif self.coord_data.rotor_type == RotorTypes.Spherical:
            raise NotImplementedError("cubic point groups require face identification")
        else:
            # search for proper rotation axes
            primary_axis = None
            group_order = None

            if self.verbose:
                print("Checking for primary rotation axis from orders:", self.group_orders)
            for o in self.group_orders:
                if all(d.planar for d in self.group_coords):
                    possible_axes = np.eye(3)
                else:
                    possible_axes = self.rotation_axis_iterator()
                for ax in possible_axes:
                    elem = RotationElement(o, ax)
                    if self.check_element(elem, verbose=False):
                        primary_axis = ax
                        group_order = o
                        elements.append(elem)
                        break
                if primary_axis is not None:
                    break

            if self.verbose:
                if primary_axis is not None:
                    print("found axis:", elements[-1])

            if group_order is not None:
                # check for any perpendicular C2 axis
                c2_axis = None
                if self.verbose:
                    print("Checking for perpendicular C2")
                for ax in self.rotation_axis_iterator():
                    if np.abs(np.dot(ax, primary_axis)) < self.tol: # has to be parallel
                        elem = RotationElement(2, ax)
                        if self.check_element(elem, verbose=False):
                            ax = np.cross(np.cross(ax, primary_axis), primary_axis) # remove any tolerances
                            c2_axis = ax
                            elements.append(elem)
                            break

                if self.verbose:
                    if c2_axis is not None:
                        print("found axis:", RotationElement(2, c2_axis))

                # check for s_h plane
                elem = ReflectionElement(axis=primary_axis)
                has_sh = self.check_element(elem)
                if has_sh:
                    elements.append(elem)
                    if c2_axis is None:
                        axes = nput.view_matrix(
                            primary_axis,
                            output_order=(0, 2, 1)
                        )
                        pg = PointGroup.from_name("Ch", group_order, axes=axes)
                    else:
                        axes = nput.view_matrix(
                            primary_axis,
                            view_vector=c2_axis,
                            output_order=(0, 2, 1)
                        )
                        pg = PointGroup.from_name("Dh", group_order, axes=axes)
                else:
                    if c2_axis is not None:
                        axes = nput.view_matrix(
                            primary_axis,
                            view_vector=c2_axis,
                            output_order=(0, 2, 1)
                        )
                        # check for s_d plane
                        elem = ReflectionElement(c2_axis)
                        has_sd = self.check_element(elem)
                        if not has_sd:
                            ax = axes[:, 2]
                            elem = ReflectionElement(ax)
                            has_sd = self.check_element(elem)
                        if has_sd:
                            elements.append(elem)
                            pg = PointGroup.from_name("Dd", group_order, axes=axes)
                        else:
                            pg = PointGroup.from_name("D", group_order, axes=axes)
                    else:
                        #TODO: throw in a special case for C2v
                        axes = nput.view_matrix(
                            primary_axis,
                            output_order=(0, 2, 1)
                        )
                        sv_axis = None
                        for ax in self.reflection_plane_iterator():
                            if np.abs(np.dot(ax, primary_axis)) < 1 - self.tol:  # not parallel
                                ax = np.cross(primary_axis, ax)
                            else:
                                ax = np.cross(np.cross(ax, primary_axis), primary_axis) # remove any tolerances
                            elem = ReflectionElement(ax)
                            if self.check_element(elem):
                                sv_axis = ax
                                elements.append(elem)
                                break
                        if sv_axis is not None:
                            pg = PointGroup.from_name("Cv", group_order, axes=axes)
                        else:
                            elem = ImproperRotationElement(2*group_order, primary_axis)
                            if self.check_element(elem):
                                pg = PointGroup.from_name("S", group_order, axes=axes)
                            else:
                                pg = PointGroup.from_name("C", group_order, axes=axes)
            else:
                if has_inversion:
                    pg = PointGroup.from_name("Ci")
                else:
                    primary_axis = np.array([0, 0, 1])
                    sd_axis = None
                    for ax in self.reflection_plane_iterator():
                        if np.abs(np.dot(ax, primary_axis)) < 1 - self.tol:  # not parallel
                            ax = np.cross(primary_axis, ax)
                        else:
                            ax = np.cross(np.cross(ax, primary_axis), primary_axis) # remove any tolerances
                        elem = ReflectionElement(ax)
                        if self.check_element(elem):
                            sd_axis = ax
                            elements.append(elem)
                            break
                    if sd_axis is not None:
                        axes = nput.view_matrix(
                            sd_axis,
                            output_order=(0, 2, 1)
                        )
                        pg = PointGroup.from_name("Cs", axes=axes)
                    else:
                        pg = PointGroup.from_name("C", 1)

        if realign:
            pg = pg.transform(self.coord_data.axes)
        return elements, pg