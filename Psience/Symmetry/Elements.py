import abc
import numpy as np
import McUtils.Numputils as nput

__all__ = [
    "IdentityElement",
    "SymmetryElement",
    "InversionElement",
    "RotationElement",
    "ReflectionElement",
    "ImproperRotationElement"
]

class SymmetryElement(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_transformation(self):
        ...
    @abc.abstractmethod
    def inverse(self):
        ...
    def __eq__(self, other):
        return np.allclose(self.get_transformation(), other.get_transformation())
    def compose(self, other):
        return ComposedSymmetryElement(self, other)

    def __matmul__(self, other):
        return self.compose(other)

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}()"

class ComposedSymmetryElement(SymmetryElement):
    def __init__(self, *bits:SymmetryElement):
        self.bits = bits
    def get_transformation(self):
        mat = self.bits[0].get_transformation()
        for b in self.bits[1:]:
            mat = mat @ b.get_transformation()
        return mat
    def inverse(self):
        return ComposedSymmetryElement(*(b.inverse() for b in reversed(self.bits)))

class IdentityElement(SymmetryElement):
    def get_transformation(self):
        return np.eye(3)
    def compose(self, other):
        return other
    def inverse(self):
        return self

class InversionElement(SymmetryElement):
    def get_transformation(self):
        return np.diag(-np.ones(3))
    def inverse(self):
        return self

    def compose(self, other):
        if isinstance(other, InversionElement):
            return IdentityElement()
        elif isinstance(other, RotationElement):
            if other.order == 2:
                return ReflectionElement(other.axis)
            else:
                return (
                    RotationElement(2, other.axis) @
                        ImproperRotationElement(other.order, other.axis, root=other.root)
                )
        elif isinstance(other, ImproperRotationElement):
            new_root = (other.order - 2*other.root) % (2*other.order)
            return RotationElement(2*other.order, -other.axis, root=new_root)
        # elif isinstance(other, (RotationElement, ImproperRotationElement)):
        #     return RotationElement(other.order, other.axis, root=other.order-other.root)
        # elif isinstance(other, ReflectionElement):
        #     return ReflectionElement(-other.axis)
        return super().compose(other)

class RotationElement(SymmetryElement):
    def __init__(self, order, axis, root=1):
        if root > 1 and order % root == 0:
            order = order // root
            root = 1
        self.root = root
        self.order = order
        self.axis = nput.vec_normalize(axis)
    def inverse(self):
        return type(self)(self.order, self.axis, self.order-self.root)

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.order})"

    def compose(self, other):
        if isinstance(other, (RotationElement, ImproperRotationElement)):
            if (
                    isinstance(other, ImproperRotationElement)
                    and other.order == 2
            ):
                return self.compose(InversionElement())
            rot_type = type(other)
            check = np.dot(other.axis, self.axis)
            if abs(check) > 1-1e-6:
                gcd = np.gcd(self.order, other.order)
                order = self.order*other.order//gcd
                root = np.sign(check)*other.root*order//other.order + self.root*order//self.order
                return rot_type(order, self.axis, root=root)
            elif self.order == 2 and check < 1e-6:
                x = self.get_transformation() @ other.get_transformation() # symmetric by construction
                vals, axes = np.linalg.eigh(x)
                if isinstance(other, ImproperRotationElement):
                    pos = np.where(vals < -.9)[0][0]
                    return InversionElement() @ RotationElement(self.order, axes[:, pos], self.root)
                else:
                    pos = np.where(vals > .9)[0][0]
                    return RotationElement(self.order, axes[:, pos], self.root)
            elif other.order == 2 and check < 1e-6:
                x = self.get_transformation() @ other.get_transformation() # symmetric by construction
                vals, axes = np.linalg.eigh(x)
                pos = np.where(vals > .9)[0][0]
                return rot_type(other.order, axes[:, pos], other.root)
            else:
                return super().compose(other)

        return other @ self

    def get_transformation(self):
        return nput.rotation_matrix(self.axis, 2*np.pi*self.root/self.order)

class ReflectionElement(SymmetryElement):
    def __init__(self, axis):
        self.axis = nput.vec_normalize(axis)

    def get_transformation(self):
        tf = nput.view_matrix(self.axis)
        return tf @ np.diag([1, -1, 1]) @ tf.T
    def inverse(self):
        return self

    def compose(self, other):
        if isinstance(other, ReflectionElement):
            check = np.dot(other.axis, self.axis)
            if abs(check) > 1-1e-6:
                return IdentityElement()
            elif check < 1e-6:
                return RotationElement(2, np.cross(self.axis, other.axis))
        elif isinstance(other, RotationElement):
            check = np.dot(other.axis, self.axis)
            if abs(check) > 1-1e-6:
                return ImproperRotationElement(other.order, other.axis)
        elif isinstance(other, ImproperRotationElement):
            if other.order == 2:
                return self.compose(InversionElement())
            check = np.dot(other.axis, self.axis)
            if abs(check) > 1 - 1e-6:
                return RotationElement(other.order, other.axis)
        elif isinstance(other, InversionElement):
            return RotationElement(2, self.axis)

        return super().compose(other)

class ImproperRotationElement(SymmetryElement):
    def __init__(self, order, axis, root=1):
        if root > 1 and order % root == 0:
            order = order // root
            root = 1
        self.root = root
        self.order = order
        self.axis = nput.vec_normalize(axis)
    def inverse(self):
        return type(self)(self.order, self.axis, self.order-self.root)

    def compose(self, other):
        if self.order == 2:
            return InversionElement().compose(other)

        if isinstance(other, RotationElement):
            return other.compose(self).inverse()
            # need to get inverse
        elif isinstance(other, ImproperRotationElement):
            if other.order == 2:
                return self.compose(InversionElement())
            check = np.dot(other.axis, self.axis)
            if abs(check) > 1 - 1e-6:
                gcd = np.gcd(self.order, other.order)
                order = self.order * other.order // gcd
                root = np.sign(check) * other.root * order // other.order + self.root * order // self.order
                return ImproperRotationElement(order, self.axis, root=root)
            else:
                return super().compose(other)
        else:
            return other.compose(self) # commute

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.order})"
    def get_transformation(self):
        rot = nput.rotation_matrix(self.axis, 2*np.pi*self.root/self.order)
        tf = nput.view_matrix(self.axis)
        reflect = tf @ np.diag([1, -1, 1]) @ tf.T
        return reflect @ rot