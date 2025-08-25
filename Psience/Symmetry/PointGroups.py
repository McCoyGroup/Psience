import abc
import enum

import numpy as np

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

    @abc.abstractmethod
    def get_name(self):
        ...
    @property
    def name(self):
        return self.get_name()
    def __repr__(self):
        return "PointGroup<{}>".format(self.name)

    @classmethod
    def from_name(cls, key, n=None):
        if isinstance(key, NamedPointGroups):
            return NamedPointGroup(key)
        elif isinstance(key, ParametrizedPointGroups):
            return ParametrizedPointGroup(key, n)
        else:
            if any(x.isdigit() for x in key):
                n = int("".join(x for x in key if x.isdigit()))
                key = "".join("".join(x for x in key if not x.isdigit()))
                return cls.from_name(key, n=n)
            try:
                key = NamedPointGroups(key)
            except ValueError:
                key = ParametrizedPointGroups(key)
            return cls.from_name(key, n=n)

    @classmethod
    def get_symmetry_element_primary_rotation(cls, elements: list[SymmetryElement]):
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
                return cls.from_name("Cs")
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
                    for e in elements:
                        if (
                                isinstance(e, ReflectionElement)
                                and np.abs(np.dot(e.axis, primary_axis)) < tol  # perp
                        ):
                            sd_plane = e
                    if sd_plane is not None:
                        return cls.from_name("Dd", group_order)
                    else:
                        return cls.from_name("D", group_order)
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

    def _get_cn_elements(self, n, axis=None):
        if axis is None:
            axis = [0, 0, 1]
        return [
            IdentityElement()
                if i == 0 else
            RotationElement(n, axis, root=i)
            for i in range(n)
        ]

    def _complete_dn_elements(self, elem):
        ...

class NamedPointGroup(PointGroup):
    def __init__(self, name:str):
        self.group = NamedPointGroups(name)
    def get_name(self):
        return self.group.value

class ParametrizedPointGroup(PointGroup):
    group: ParametrizedPointGroups
    n: int
    def __init__(self, name:str, n):
        self.group = ParametrizedPointGroups(name)
        self.n = n
    def get_name(self):
        return self.group.value[:1] + str(self.n) + self.group.value[1:]



