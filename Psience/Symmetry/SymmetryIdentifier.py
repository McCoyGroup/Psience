
import abc
import numpy as np
import enum
import itertools
import McUtils.Numputils as nput
import McUtils.Combinatorics as comb

__all__ = [
    "PointGroupIdentifier"
]

class SymmetryElement(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_transformation(self):
        ...

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}()"

class InversionElement(SymmetryElement):
    def get_transformation(self):
        return np.diag(-np.ones(3))

class ProperRotationElement(SymmetryElement):
    def __init__(self, order, axis):
        self.order = order
        self.axis = axis

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.order})"

    def get_transformation(self):
        return nput.rotation_matrix(self.axis, 2*np.pi/self.order)

class ReflectionElement(SymmetryElement):
    def __init__(self, axis):
        self.axis = axis

    def get_transformation(self):
        tf = nput.view_matrix(self.axis)
        return tf @ np.diag([1, -1, 1]) @ tf.T

class ImproperRotationElement(SymmetryElement):
    def __init__(self, order, axis):
        self.order = order
        self.axis = axis
    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.order})"
    def get_transformation(self):
        rot = nput.rotation_matrix(self.axis, 2*np.pi/self.order)
        tf = nput.view_matrix(self.axis)
        reflect = tf @ np.diag([1, -1, 1]) @ tf.T
        return reflect @ rot

class RotorTypes:
    Atom = "atom"
    Linear = "linear"
    Planar = "planar"
    Oblate = "oblate"
    Prolate = "prolate"
    Spherical = "spherical"
    Asymmetric = "asymmetric"

class PointGroup(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_name(self):
        ...
    @property
    def name(self):
        return self.get_name()
    def __repr__(self):
        return "PointGroup<{}>".format(self.name)

class NamedPointGroup(PointGroup):
    base_name: str
    def get_name(self):
        return self.base_name

class PointGroupIh(NamedPointGroup): base_name = "Ih"
class PointGroupI(NamedPointGroup): base_name = "I"
class PointGroupOh(NamedPointGroup): base_name = "Oh"
class PointGroupO(NamedPointGroup): base_name = "O"
class PointGroupTh(NamedPointGroup): base_name = "Th"
class PointGroupTd(NamedPointGroup): base_name = "Td"
class PointGroupT(NamedPointGroup): base_name = "T"
class PointGroupCs(NamedPointGroup): base_name = "Cs"
class PointGroupCi(NamedPointGroup): base_name = "Ci"

class ParametrizedPointGroup(PointGroup):
    base_name: str
    n: int
    def __init__(self, n):
        self.n = n
    def get_name(self):
        return self.base_name + "_" + str(self.n)

class PointGroupC(ParametrizedPointGroup): base_name = "C"
class PointGroupCh(ParametrizedPointGroup): base_name = "Ch"
class PointGroupCv(ParametrizedPointGroup): base_name = "Cv"
class PointGroupD(ParametrizedPointGroup): base_name = "D"
class PointGroupDh(ParametrizedPointGroup): base_name = "Dh"
class PointGroupDd(ParametrizedPointGroup): base_name = "Dd"
class PointGroupS(ParametrizedPointGroup): base_name = "S"

class NamedPointGroups(enum.Enum):
    Ih = "Ih"
    I = "I"
    Oh = "Oh"
    O = "O"
    Th = "Th"
    Td = "Td"
    T = "T"
    Cs = "Cs"
class ParametrizedPointGroup(enum.Enum):
    C = "C"
    Ch = "Ch"
    Cv = "Cv"
    D = "D"
    Dh = "Dh"
    Dd = "Dd"
    S = "S"

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
    def __init__(self, coords, masses=None, groups=None, tol=1e-8, mass_tol=1, mom_tol=1):
        self.mass_tol = mass_tol
        self.mom_tol = mom_tol
        self.tol = tol
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
                    if np.max(np.abs(d - d2)) < self.tol:
                        subgroup.add(i+j+1)
                if new_sg and len(subgroup) > 1:
                    subgroups.append(subgroup)
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
                for s in subgroups
            )
        return full_groups

    def get_rotor_type(self, moms):
        mom_deg, counts = np.unique(np.round(moms, self.mom_tol), return_counts=True)
        planar = np.abs((moms[0] + moms[1]) - moms[2]) < self.mom_tol
        if np.min(mom_deg) < 1e-6:
            if len(np.where(moms < 1e-6)[0]) == 2:
                type = RotorTypes.Atom
            else:
                type = RotorTypes.Linear
        elif tuple(counts) == (3,):
            type = RotorTypes.Spherical
        elif tuple(counts) == (1,2):
            type = RotorTypes.Prolate
        elif tuple(counts) == (2,1):
            type = RotorTypes.Oblate
        else:
            type = RotorTypes.Asymmetric

        return type, planar

    def get_group_orders(self):
        group_factors = None
        for g in self.groups:
            if len(g) < 2: continue
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
                group_factors = group_factors.intersection(subfactors)
        return list(sorted(group_factors, key=lambda n:-n))

    def check_element(self, elem:SymmetryElement):
        tf = elem.get_transformation()
        for g in self.groups:
            coords = self.coord_data.coords[g,]
            new_coords = coords @ tf.T
            dm = np.linalg.norm(coords[:, np.newaxis] - new_coords[np.newaxis, :], axis=-1)
            min_dists = np.argmin(dm, axis=1)
            max_min_dist = np.max(dm[np.arange(len(coords)), min_dists])
            if max_min_dist > len(coords) * self.tol:
                return False
            diff_comp = coords - new_coords[min_dists,]
            if np.max(np.abs(diff_comp)) > self.tol:
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
        for g in self.groups:
            for n, i in enumerate(g):
                ax = self.coord_data.coords[i]
                yield ax
                for j in g[n + 1:]:
                    ax = (self.coord_data.coords[i] + self.coord_data.coords[j]) / 2
                    yield ax

    def reflection_plane_iterator(self):
        for e in np.eye(3):
            yield e
        for g in self.groups:
            for i,j in itertools.combinations(g, 2):
                ax = self.coord_data.coords[i] - self.coord_data.coords[j]
                yield ax

    def identify_point_group(self):
        elements = []
        elem = InversionElement()
        has_inversion = self.check_element(elem)
        if has_inversion:
            elements.append(elem)

        if self.coord_data.rotor_type == RotorTypes.Linear:
            if has_inversion:
                point_group = PointGroupD(-1)
            else:
                point_group = PointGroupC(-1)
        elif self.coord_data.rotor_type == RotorTypes.Spherical:
            raise NotImplementedError("cubic point groups require face identification")
        else:
            # search for proper rotation axes
            primary_axis = None
            group_order = None
            for o in self.group_orders:
                if all(d.planar for d in self.group_coords):
                    possible_axes = np.eye(3)
                else:
                    possible_axes = self.rotation_axis_iterator()
                for ax in possible_axes:
                    elem = ProperRotationElement(o, ax)
                    if self.check_element(elem):
                        primary_axis = ax
                        group_order = o
                        elements.append(elem)
                        break
                if primary_axis is not None:
                    break

            if group_order is not None:
                # check for any perpendicular C2 axis
                c2_axis = None
                for ax in self.rotation_axis_iterator():
                    if np.abs(np.dot(ax, primary_axis)) < self.tol: # has to be parallel
                        elem = ProperRotationElement(2, ax)
                        if self.check_element(elem):
                            c2_axis = ax
                            elements.append(elem)
                            break

                # check for s_h plane
                elem = ReflectionElement(axis=primary_axis)
                has_sh = self.check_element(elem)
                if has_sh:
                    elements.append(elem)
                    if c2_axis is None:
                        point_group = PointGroupCh(group_order)
                    else:
                        point_group = PointGroupDh(group_order)
                else:
                    if c2_axis is not None:
                        # check for s_d plane
                        ax = np.cross(primary_axis, c2_axis)
                        elem = ReflectionElement(ax)
                        has_sd = self.check_element(elem)
                        if has_sd is not None:
                            point_group = PointGroupDh(group_order)
                        else:
                            point_group = PointGroupD(group_order)
                    else:
                        sv_axis = None
                        for ax in self.reflection_plane_iterator():
                            if np.abs(np.dot(ax, primary_axis)) < 1 - self.tol:  # not parallel
                                ax = np.cross(primary_axis, ax)
                            elem = ReflectionElement(ax)
                            if self.check_element(elem):
                                sv_axis = ax
                                elements.append(elem)
                                break
                        if sv_axis is not None:
                            point_group = PointGroupCv(group_order)
                        else:
                            elem = ImproperRotationElement(2*group_order, primary_axis)
                            if self.check_element(elem):
                                point_group = PointGroupS(group_order)
                            else:
                                point_group = PointGroupC(group_order)
            else:
                if has_inversion:
                    point_group = PointGroupCi()
                else:
                    primary_axis = np.array([0, 0, 1])
                    sd_axis = None
                    for ax in self.reflection_plane_iterator():
                        if np.abs(np.dot(ax, primary_axis)) < 1 - self.tol:  # not parallel
                            ax = np.cross(primary_axis, ax)
                        elem = ReflectionElement(ax)
                        if self.check_element(elem):
                            sd_axis = ax
                            elements.append(elem)
                            break
                    if sd_axis is not None:
                        point_group = PointGroupCs()
                    else:
                        point_group = PointGroupC(1)

        return elements, point_group