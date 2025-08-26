import abc
import enum
import numpy as np
import McUtils.Devutils as dev
import McUtils.Numputils as nput
import McUtils.Combinatorics as comb

from .Elements import *

__all__ = [
    "NamedPointGroups",
    "ParametrizedPointGroups",
    "PointGroup"
]

class NamedPointGroups(enum.Enum):
    Ih = "Ih"
    I = "I"
    Oh = "Oh"
    O = "O"
    Th = "Th"
    Td = "Td"
    T = "T"
    Cs = "Cs"
class ParametrizedPointGroups(enum.Enum):
    C = "C"
    Ch = "Ch"
    Cv = "Cv"
    D = "D"
    Dh = "Dh"
    Dd = "Dd"
    S = "S"


class PointGroup(metaclass=abc.ABCMeta):
    group: 'NamedPointGroups|ParametrizedPointGroups'

    def __init__(self, character_table=None, elements=None, axes=None):
        self._elements = elements
        self._character_table = character_table
        self._axes = axes

    def get_modification_kwargs(
            self,
            character_table=dev.default,
            elements=dev.default,
            axes=dev.default
    ):
        return dict(
            character_table=self.character_table if dev.is_default(character_table) else character_table,
            elements=elements if dev.is_default(elements) else elements,
            axes=axes if dev.is_default(axes) else axes,
        )
    def modify(self,
               character_table=dev.default,
               elements=dev.default,
               axes=dev.default
               ):
        return type(self)(**self.get_modification_kwargs(
            character_table=character_table,
            elements=elements,
            axes=axes
        ))

    @abc.abstractmethod
    def get_name(self):
        ...
    @property
    def name(self):
        return self.get_name()
    def __repr__(self):
        return "PointGroup<{}>".format(self.name)

    @classmethod
    def from_name(cls, key, n=None, **etc):
        if isinstance(key, NamedPointGroups):
            return NamedPointGroup(key, **etc)
        elif isinstance(key, ParametrizedPointGroups):
            return ParametrizedPointGroup(key, n, **etc)
        else:
            if any(x.isdigit() for x in key):
                n = int("".join(x for x in key if x.isdigit()))
                key = "".join("".join(x for x in key if not x.isdigit()))
                return cls.from_name(key, n=n, **etc)
            try:
                key = NamedPointGroups(key)
            except ValueError:
                key = ParametrizedPointGroups(key)
            return cls.from_name(key, n=n, **etc)

    @classmethod
    def get_symmetry_element_primary_rotation(cls, elements: 'Iterable[SymmetryElement]'):
        rotors = [e for e in elements if isinstance(e, RotationElement)]
        if len(rotors) == 0:
            return None, None
        else:
            orders, loc, counts = np.unique([e.order for e in rotors], return_counts=True, return_index=True)
            return rotors[loc[-1]], counts[-1]

    @classmethod
    def from_symmetry_elements(cls, elements: list[SymmetryElement], tol=1e-2):
        primary_axis, counts = cls.get_symmetry_element_primary_rotation(elements)
        if primary_axis is None:
            if any(isinstance(e, InversionElement) for e in elements):
                return cls.from_name("Ci")
            elif any(isinstance(e, ReflectionElement) for e in elements):
                reflections = [e for e in elements if isinstance(e, ReflectionElement)]
                axes = nput.view_matrix(reflections[0].axis)[:, (0, 2, 1)]
                return cls.from_name("Cs", axes=axes)
            else:
                return cls.from_name("C", 1)
        elif counts > 1:  # high symmetry
            raise NotImplementedError(...)
        else:
            group_order = primary_axis.order
            primary_axis = primary_axis.axis
            c2_axis = None
            for e in elements:
                if (
                        isinstance(e, RotationElement)
                        and e.order == 2
                        and np.abs(np.dot(e.axis, primary_axis)) < tol  # perp
                ):
                    c2_axis = e
                    break

            sh_plane = None
            for e in elements:
                if (
                        isinstance(e, ReflectionElement)
                        and np.abs(np.dot(e.axis, primary_axis)) > 1 - tol  # parallel
                ):
                    sh_plane = e
                    break

            if c2_axis is not None:
                if sh_plane is not None:
                    return cls.from_name("Dh", group_order)
                else:
                    sd_plane = None
                    axes = nput.view_matrix(
                        primary_axis.axis,
                        view_vector=c2_axis.axis
                    )[:, (0, 2, 1)]
                    for e in elements:
                        if (
                                isinstance(e, ReflectionElement)
                                and np.abs(np.dot(e.axis, primary_axis)) < tol  # perp
                        ):
                            sd_plane = e
                    if sd_plane is not None:
                        return cls.from_name("Dd", group_order, axes=axes)
                    else:
                        return cls.from_name("D", group_order, axes=axes)
            else:
                if sh_plane is not None:
                    return cls.from_name("Ch", group_order)
                else:
                    sv_plane = None
                    for e in elements:
                        if (
                                isinstance(e, ReflectionElement)
                                and np.abs(np.dot(e.axis, primary_axis)) < tol  # perp
                        ):
                            sv_plane = e
                    if sv_plane is not None:
                        return cls.from_name("Cv", group_order)
                    else:
                        if any(isinstance(e, ImproperRotationElement) for e in elements):
                            return cls.from_name("S", group_order)
                        else:
                            return cls.from_name("C", group_order)

    def get_symmetry_elements(self):
        ct = self.get_character_table()
        mats = ct.space_representation
        if mats is None:
            raise NotImplementedError(f"no implementation for operator matrices for point group {self}")
        base_elems = tuple(SymmetryElement.from_transformation_matrix(x) for x in mats)
        if self._axes is not None:
            base_ax = self.get_axes(base_elems)
            tf = self._axes @ base_ax.T
            if abs(np.linalg.det(tf)) < 1e-6:
                raise ValueError("bad transformation")
            base_elems = [e.transform(tf) for e in base_elems]
        return base_elems

    @property
    def character_table(self) -> comb.CharacterTable:
        if self._character_table is None:
            self._character_table = self.get_character_table()
        return self._character_table

    @property
    def elements(self) -> 'tuple[SymmetryElement]':
        if self._elements is None:
            self._elements = self.get_symmetry_elements()
        return self._elements

    @abc.abstractmethod
    def get_character_table(self):
        ...

    def get_axes(self, elements=None):
        # could use point group identity for short circuiting, but don't want to

        if elements is None:
            elements = self.elements

        primary_axis = None
        secondary_axis = None

        primary, count = self.get_symmetry_element_primary_rotation(elements)
        if primary is None:
            # search for reflection plane
            planes = [e for e in elements if isinstance(e, ReflectionElement)]
            if len(planes) == 0:
                return np.eye(3)
            primary_axis = planes[0].axis
            if len(planes) > 1: # I don't think this is possible
                secondary_axis = planes[1].axis
        else:
            primary_axis = primary.axis
            c2_axes = [e for e in elements if isinstance(e, RotationElement) and e.order == 2]
            perp_axes = [
                e for e in c2_axes
                if abs(np.dot(primary_axis, e.axis)) < 1e-2
            ]
            if len(perp_axes) == 0:
                # search for vertical reflection planes
                planes = [e for e in elements if isinstance(e, ReflectionElement)]
                v_planes = [
                    e for e in planes
                    if abs(np.dot(primary_axis, e.axis)) < 1e-2
                ]
                if len(v_planes) > 0:
                    secondary_axis = v_planes[0].axis
                else:
                    non_par = [
                        e for e in c2_axes
                        if abs(np.dot(primary_axis, e.axis)) < .9
                    ]
                    if len(non_par) > 0:
                        secondary_axis = non_par[0].axis
            else:
                secondary_axis = perp_axes[0].axis

        if primary_axis is None:
            raise ValueError("no primary axis found") # how??

        if secondary_axis is None:
            secondary_axis = [1, 0, 0]
            if abs(np.dot(secondary_axis, primary_axis)) > 1 - 1e-2:
                secondary_axis = [0, 0, 1]

        return nput.view_matrix(primary_axis, secondary_axis, output_order=(0, 2, 1))
        # tertiary_axis = nput.vec_crosses(primary_axis, secondary_axis, normalize=True)
        # secondary_axis = nput.vec_crosses(tertiary_axis, primary_axis, normalize=True)
        #
        # return np.array([secondary_axis, tertiary_axis, primary_axis]).T

    @property
    def axes(self):
        if self._axes is None:
            self._axes = self.get_axes()
        return self._axes

    def align(self, axes):
        tf = axes @ self.axes.T
        return self.modify(
            axes=axes,
            elements=[e.transform(tf) for e in self.elements]
        )
    def transform(self, tf):
        return self.modify(
            axes=tf @ self.axes,
            elements=[e.transform(tf) for e in self.elements]
        )

    def plot(self,
             figure=None,
             elements=None,
             origin=None,
             inversion_styles=None,
             rotation_styles=None,
             reflection_styles=None,
             improper_rotation_styles=None,
             **opts):
        from McUtils.Plots import Graphics3D

        if figure is None:
            figure = Graphics3D(backend='x3d')

        if elements is None:
            elements = self.elements
        if inversion_styles is None:
            inversion_styles = {}
        if rotation_styles is None:
            rotation_styles = {}
        if reflection_styles is None:
            reflection_styles = {}
        if improper_rotation_styles is None:
            improper_rotation_styles = rotation_styles
        for e in elements:
            base_styles = (
                inversion_styles
                    if isinstance(e, InversionElement) else
                reflection_styles
                    if isinstance(e, ReflectionElement) else
                rotation_styles
                    if isinstance(e, RotationElement) else
                improper_rotation_styles
                    if isinstance(e, ImproperRotationElement) else
                {}
            )
            if base_styles is False: continue
            styles = dict(
                opts,
                **base_styles
            )
            e.plot(figure, origin=origin, **styles)

        return figure
class NamedPointGroup(PointGroup):
    def __init__(self, name:'str|NamedPointGroups', character_table=None, elements=None, axes=None):
        super().__init__(character_table=character_table, elements=elements, axes=axes)
        self.group = NamedPointGroups(name)
    def modify(self,
               character_table=dev.default,
               elements=dev.default,
               axes=dev.default
               ):
        return type(self)(
            self.group,
            **self.get_modification_kwargs(
            character_table=character_table,
            elements=elements,
            axes=axes
        ))
    def get_name(self):
        return self.group.value

    def get_character_table(self):
        return comb.CharacterTable.fixed_size_point_group(self.group.value)

class ParametrizedPointGroup(PointGroup):
    group: ParametrizedPointGroups
    n: int
    def __init__(self, name:'str|ParametrizedPointGroups', n, character_table=None, elements=None, axes=None):
        super().__init__(character_table=character_table, elements=elements, axes=axes)
        self.group = ParametrizedPointGroups(name)
        self.n = n
    def modify(self,
               character_table=dev.default,
               elements=dev.default,
               axes=dev.default
               ):
        return type(self)(
            self.group,
            self.n,
            **self.get_modification_kwargs(
                character_table=character_table,
                elements=elements,
                axes=axes
            ))
    def get_name(self):
        return self.group.value[:1] + str(self.n) + self.group.value[1:]

    def get_character_table(self):
        return comb.CharacterTable.point_group(self.group.value, self.n)


