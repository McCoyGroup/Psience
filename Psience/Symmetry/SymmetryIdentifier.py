
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
    "PointGroupIdentifier",
    "identify_symmetry_equivalent_atoms",
    "identify_point_group"
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

def identify_symmetry_equivalent_atoms(coords,
                                       masses=None,
                                       base_groups=None,
                                       mass_tol=1,
                                       tol=1e-2
                                       ):
    if base_groups is None:
        if masses is not None:
            (_, base_groups), _ = nput.group_by(np.arange(len(masses)), np.round(masses/mass_tol))
        else:
            base_groups = [np.arange(len(coords))]
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
                if np.max(np.abs(d - d2)) < tol:
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
    return full_groups

class PointGroupIdentifier:
    def __init__(self, coords, masses=None, groups=None,
                 tol=1e-2, mass_tol=1, mom_tol=1,
                 grouping_tol=1e-2,
                 verbose=False
                 ):
        self.mass_tol = mass_tol
        self.mom_tol = mom_tol
        self.tol = tol
        self.grouping_tol = grouping_tol
        self.verbose = verbose
        self.groups = identify_symmetry_equivalent_atoms(coords,
                                                         masses=masses,
                                                         tol=grouping_tol,
                                                         mass_tol=mass_tol)
        if masses is None:
            masses = np.ones(len(coords))
        masses = np.asanyarray(masses)
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

        seen_axes = []
        for g in self.groups:
            for n, i in enumerate(g):
                ax = self.coord_data.coords[i]
                norm = np.linalg.norm(ax)
                if norm > 1e-2:
                    ax = ax / norm
                    if all(abs(np.dot(ax, s)) < 1 - 1e-2 for s in seen_axes):
                        seen_axes.append(ax)
                        yield ax
                for j in g[n + 1:]:
                    ax = (self.coord_data.coords[i] + self.coord_data.coords[j]) / 2
                    norm = np.linalg.norm(ax)
                    if norm > 1e-2:
                        ax = ax / norm
                        if all(abs(np.dot(ax, s)) < 1 - 1e-2 for s in seen_axes):
                            seen_axes.append(ax)
                            yield ax

    def reflection_plane_iterator(self):
        for e in np.eye(3):
            yield e
        for g in self.groups:
            for i,j in itertools.combinations(g, 2):
                ax = self.coord_data.coords[i] - self.coord_data.coords[j]
                yield ax

    def rotation_face_iterator(self):
        for g in self.groups:
            if len(g) < 3: continue
            for i,j,k in itertools.combinations(g, 3):
                v1 = self.coord_data.coords[i] - self.coord_data.coords[j]
                v2 = self.coord_data.coords[i] - self.coord_data.coords[k]
                v3 = self.coord_data.coords[j] - self.coord_data.coords[k]
                norm1 = np.linalg.norm(v1)
                norm2 = np.linalg.norm(v2)
                norm3 = np.linalg.norm(v3)
                if abs(norm1 - norm2) < 1e-2 and abs(norm1 - norm3) < 1e-2:
                    yield np.cross(v1, v2)
                    
    def embed_point_group(self, point_group:'PointGroup|list[SymmetryElement]'):
        axes = self.find_point_group_alignment_axes(point_group)
        if isinstance(point_group, PointGroup):
            return point_group.align(axes)
        else:
            base_axes = PointGroup.get_axes_from_symmetry_elements(point_group)
            tf = axes @ base_axes.T
            return [e.transform(tf) for e in point_group]

    def find_point_group_alignment_axes(self, point_group:'PointGroup|list[SymmetryElement]'):
        pg_elements = None
        # pg_axes = None
        if isinstance(point_group, PointGroup):
            pg_elements = point_group.elements
            # pg_axes = point_group.axes
        else:
            pg_elements = point_group
            # axes = PointGroup.get_axes_from_symmetry_elements(point_group)

        #TODO: cut down on duplicaiton by having `PointGroup.get_axes_from_symmetry_elements` return
        #      the elements used

        pg_primary_axis, counts = PointGroup.get_symmetry_element_primary_rotation(pg_elements)
        if pg_primary_axis is not None:
            primary_axis = None
            for ax in self.rotation_axis_iterator():
                test_elem = RotationElement(pg_primary_axis.order, ax)
                if self.check_element(test_elem):
                    primary_axis = test_elem
                    break
            else:
                raise ValueError(f"couldn't find {pg_primary_axis} type element")

            pg_perp_c2 = None
            for e in pg_elements:
                if (
                        isinstance(e, RotationElement)
                        and e.order == 2
                        and np.abs(np.dot(e.axis, pg_primary_axis.axis)) < 1e-2  # perp
                ):
                    pg_perp_c2 = e
                    break

            if pg_perp_c2 is not None:
                perp_c2 = None
                for ax in self.rotation_axis_iterator():
                    elem = RotationElement(2, ax)
                    if (
                            np.dot(elem.axis, primary_axis.axis) < 1e-2
                            and self.check_element(elem)
                    ):
                        perp_c2 = elem
                        break
                else:
                    raise ValueError("couldn't find perpendicular C2 axis")

                axes = nput.view_matrix(
                    primary_axis.axis,
                    view_vector=perp_c2.axis,
                    output_order=(0, 2, 1)
                )
            else:
                pg_sv_plane = None
                for e in pg_elements:
                    if (
                        isinstance(e, ReflectionElement)
                        and abs(np.dot(pg_primary_axis.axis, e.axis)) < 1e-2
                    ):
                        pg_sv_plane = e
                        break

                if pg_sv_plane is not None:
                    sv_plane = None
                    for ax in self.reflection_plane_iterator():
                        elem = ReflectionElement(ax)
                        if (
                                abs(np.dot(primary_axis.axis, elem.axis)) < 1e-2
                            and self.check_element(elem)
                        ):
                            sv_plane = elem
                            break
                    else:
                        raise ValueError("couldn't find vertical reflection plane")
                    axes = nput.view_matrix(
                        primary_axis.axis,
                        view_vector=sv_plane.axis,
                        output_order=(0, 2, 1)
                    )
                else:
                    # look for any c2
                    pg_c2 = None
                    for e in pg_elements:
                        if (
                                isinstance(e, RotationElement)
                            and e.order == 2
                            and abs(np.dot(pg_primary_axis.axis, e.axis)) < 1-1e-2
                        ):
                            pg_c2 = e
                            break

                    if pg_c2 is not None:
                        c2_axis = None
                        for ax in self.rotation_axis_iterator():
                            elem = RotationElement(2, ax)
                            if (
                                    abs(np.dot(primary_axis.axis, elem.axis)) < 1 - 1e-2
                                and self.check_element(elem)
                            ):
                                c2_axis = elem
                                break
                        else:
                            raise ValueError("couldn't find non-colinear C2 axis")

                        axes = nput.view_matrix(
                            primary_axis.axis,
                            view_vector=c2_axis.axis,
                            output_order=(0, 2, 1)
                        )
                    else:
                        axes = nput.view_matrix(
                            primary_axis.axis,
                            output_order=(0, 2, 1)
                        )
        else:
            # search for reflection plane in elements
            pg_plane = None
            for e in pg_elements:
                if isinstance(e, ReflectionElement):
                    pg_plane = e
                    break
            if pg_plane is None:
                # no orientation to speak of
                axes = np.eye(3)
            else:
                plane = None
                for ax in self.reflection_plane_iterator():
                    test_elem = ReflectionElement(ax)
                    if self.check_element(test_elem):
                        plane = test_elem
                        break
                else:
                    raise ValueError(f"couldn't find {pg_plane} type element")

                axes = nput.view_matrix(
                    plane.axis,
                    output_order=(0, 2, 1)
                )

        return self.coord_data.axes @ axes


    def identify_point_group(self, realign=True):
        elements = []
        elem = InversionElement()
        has_inversion = self.check_element(elem)
        if has_inversion:
            elements.append(elem)

        if self.verbose:
            print("Using groups:", self.groups)
            print("Rotor type:", self.coord_data.rotor_type)
        if self.coord_data.rotor_type == RotorTypes.Linear:
            if has_inversion:
                pg = PointGroup.from_name("D", -1)
            else:
                pg = PointGroup.from_name("C", -1)
        elif self.coord_data.rotor_type == RotorTypes.Spherical:
            n_c2 = 0
            for ax in self.rotation_axis_iterator():
                elem = RotationElement(2, ax)
                if self.check_element(elem, verbose=False):
                    n_c2 += 1
                    elements.append(elem)
                    if n_c2 >= 15: break
            if n_c2 >= 15:
                c5_axis = None
                for ax in self.rotation_face_iterator():
                    elem = RotationElement(5, ax)
                    if self.check_element(elem):
                        c5_axis = elem.axis
                        elements.append(elem)
                        break
                else:
                    raise ValueError("couldn't identify C4 axis for O-type structure")
                axes = nput.view_matrix(
                    c5_axis,
                    view_vector=elements[-2].axis,  # arbitrary C2
                    output_order=(0, 2, 1)
                )
                if has_inversion:
                    return PointGroup.from_name("Ih", axes=axes)
                else:
                    return PointGroup.from_name("I", axes=axes)
            elif n_c2 >= 9:
                c4_axis = None
                for ax in self.rotation_face_iterator():
                    elem = RotationElement(4, ax)
                    if self.check_element(elem):
                        c4_axis = elem.axis
                        elements.append(elem)
                        break
                else:
                    raise ValueError("couldn't identify C4 axis for O-type structure")
                axes = nput.view_matrix(
                    c4_axis,
                    view_vector=elements[-2].axis, # arbitrary C2
                    output_order=(0, 2, 1)
                )
                if has_inversion:
                    pg = PointGroup.from_name("Oh", axes=axes)
                else:
                    pg = PointGroup.from_name("O", axes=axes)
            elif n_c2 >= 3:
                c3_axis = None
                n_c3 = 0
                for ax in self.rotation_face_iterator():
                    elem = RotationElement(3, ax)
                    if self.check_element(elem):
                        n_c3 += 1
                        elements.append(elem)
                        c3_axis = elem.axis
                        if has_inversion or n_c3 > 4:
                            break
                if c3_axis is None:
                    raise ValueError("couldn't identify C3 axis for T-type structure")
                # print(elements[-(n_c3 + 1)])
                # print(elements)
                if has_inversion:
                    axes = nput.view_matrix(
                        c3_axis,
                        view_vector=elements[-(n_c3 + 2)].axis,  # arbitrary C2
                        output_order=(0, 2, 1)
                    )
                    pg = PointGroup.from_name("Th", axes=axes)
                elif n_c3 >= 4:
                    # search for Td plane by noting it C3 and at least one C2
                    sd_plane = None
                    for e in elements:
                        if isinstance(e, RotationElement) and e.order == 2:
                            ax = np.cross(c3_axis, e.axis)
                            if np.linalg.norm(ax) > 1e-3:
                                elem = ReflectionElement(ax)
                                if self.check_element(elem):
                                    sd_plane = elem
                                    elements.append(sd_plane)
                                    break
                    if sd_plane is None:
                        raise ValueError("no vertical reflection plane found?")
                        axes = nput.view_matrix(
                            c3_axis,
                            view_vector=elements[-(n_c3 + 2)].axis,  # arbitrary C2
                            output_order=(0, 2, 1)
                        )
                    else:
                        axes = nput.view_matrix(
                            c3_axis,
                            view_vector=sd_plane.axis,
                            output_order=(0, 2, 1)
                        )
                    pg = PointGroup.from_name("Td", axes=axes)
                else:
                    axes = nput.view_matrix(
                        c3_axis,
                        view_vector=elements[-(n_c3 + 2)].axis,  # arbitrary C2
                        output_order=(0, 2, 1)
                    )
                    pg = PointGroup.from_name("T", axes=axes)
            else:
                raise ValueError(f"bad number of C2 axes ({n_c2})")
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

def identify_point_group(
        coords, masses=None, groups=None,
        tol=1e-8, mass_tol=1, mom_tol=1,
        grouping_tol=1e-2,
        realign=True,
        verbose=False
):
    return PointGroupIdentifier(coords,
                                masses=masses, groups=groups,
                                tol=tol, mass_tol=mass_tol, mom_tol=mom_tol,
                                grouping_tol=grouping_tol,
                                verbose=verbose
                                ).identify_point_group(realign=realign)
