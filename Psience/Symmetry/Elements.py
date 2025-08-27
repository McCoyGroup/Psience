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
        raise NotImplementedError("abstract")
    @abc.abstractmethod
    def inverse(self):
        raise NotImplementedError("abstract")
    def __eq__(self, other):
        return np.allclose(self.get_transformation(), other.get_transformation())
    def compose(self, other):
        return ComposedSymmetryElement(self, other)

    def __matmul__(self, other):
        return self.compose(other)

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}()"

    @classmethod
    def from_transformation_matrix(cls, x, max_rotation_order=60):
        # check unitary
        x = np.asanyarray(x)
        y = x @ x.T
        if np.trace(y) < 3 - 1e-2:
            raise ValueError("transformation matrix not unitary")

        t = np.trace(x)
        if t < -3 + 1e-2:
            return InversionElement()
        elif t > 3 - 1e-2:
            return IdentityElement()

        d = np.linalg.det(x)
        if d > 1-1e-2:
            # pure rotation
            ang, ax = nput.extract_rotation_angle_axis(x)
            rational = ang / (2*np.pi)
            if max_rotation_order is None:
                return RotationElement(1, ax, rational)

            for o in range(1, max_rotation_order+1):
                root = rational * o
                if abs(root - np.round(root)) < 1e-6:
                    return RotationElement(o, ax, int(np.round(root)))
            else:
                raise ValueError(
                    f"angle ratio {ang} doesn't correspond to a rational number up to order {max_rotation_order} rotations"
                )
        else:
            if t > -3 + 1e-2:
                # reflection or improper rotation
                # check symmetric
                if np.max(np.abs(x - x.T)) < 1e-8:
                    # reflection
                    eigs, ax = np.linalg.eigh(x)
                    pos = np.where(eigs < -1+1e-2)[0][0]
                    ax = ax[:, pos]
                    return ReflectionElement(ax)
                else:
                    eigs, ax = np.linalg.eig(x)
                    pos = np.where(np.real(eigs) < -1+1e-2)[0][0]
                    ax = np.real(ax[:, pos])
                    x2 = ReflectionElement(ax).get_transformation() @ x
                    ang, ax2 = nput.extract_rotation_angle_axis(x2)
                    #TODO: just double check ax == ax2
                    rational = ang / (2 * np.pi)
                    if max_rotation_order is None:
                        return ImproperRotationElement(1, ax, rational)

                    for o in range(1, max_rotation_order + 1):
                        root = rational * o
                        if abs(root - np.round(root)) < 1e-6:
                            return ImproperRotationElement(o, ax, int(np.round(root)))
                    else:
                        raise ValueError(
                            f"angle ratio {ang} doesn't correspond to a rational number up to order {max_rotation_order} improper rotations"
                        )
            else:
                return InversionElement()
    @abc.abstractmethod
    def transform(self, tf):
        ...

    @abc.abstractmethod
    def plot(self, figure, **graphics_options):
        ...

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
    def transform(self, tf):
        return type(self)(*(b.transform(tf) for b in reversed(self.bits)))
    def plot(self, figure, **graphics_options):
        symm = SymmetryElement.from_transformation_matrix(self.get_transformation(), max_rotation_order=None)
        return symm.plot(figure=figure, **graphics_options)

class IdentityElement(SymmetryElement):
    def get_transformation(self):
        return np.eye(3)
    def compose(self, other):
        return other
    def inverse(self):
        return self
    def transform(self, tf):
        return self

    def plot(self, figure, **graphics_options):
        ...

class InversionElement(SymmetryElement):
    def get_transformation(self):
        return np.diag(-np.ones(3))
    def inverse(self):
        return self
    def transform(self, tf):
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

    def plot(self, figure, *, origin=None, point_type=None, radius=.1, color='red', **graphics_options):
        from McUtils.Plots import Sphere
        if point_type is None:
            point_type = Sphere
        if origin is None:
            origin = (0, 0, 0)
        return point_type(
            origin,
            radius=radius,
            color=color,
            **graphics_options
        ).plot(figure)

class RotationElement(SymmetryElement):
    def __init__(self, order, axis, root=1):
        ax, norm = nput.vec_normalize(axis, return_norms=True)
        if norm < 1e-6: raise ValueError("can't have rotation element with no axis")
        if root > 1 and order % root == 0:
            order = order // root
            root = 1
        self.root = root
        self.order = order
        self.axis = ax
    def inverse(self):
        return type(self)(self.order, self.axis, self.order-self.root)
    def transform(self, tf):
        return type(self)(self.order, tf @ self.axis, root=self.root)

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.root}/{self.order}, {self.axis})"

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
                return rot_type(order, self.axis, root=int(np.round(root)))
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

    def plot(self, figure, *,
             origin=None,
             line_type=None,
             disk_type=None,
             color='black',
             radius=2,
             spoke_radius=.3,
             disk_color='black',
             disk_transparency=.8,
             **graphics_options
             ):
        from McUtils.Plots import Line, Disk
        if line_type is None:
            line_type = Line
        if origin is None:
            origin = (0, 0, 0)
        origin = np.asarray(origin)
        points1 = origin - radius * self.axis
        points2 = origin + radius * self.axis
        objects = []
        axis_line = line_type(
                [points1, points2],
                color=color,
                **graphics_options
            )
        objects.append(axis_line)
        if spoke_radius is not None:
            if disk_color is not None:
                if disk_type is None:
                    disk_type = Disk
                disk = disk_type(
                    origin,
                    radius=spoke_radius,
                    normal=self.axis,
                    color=None,
                    line_color=disk_color
                    # transparency=disk_transparency
                )
                objects.append(disk)

            tf = self.get_transformation()
            baseline = nput.view_matrix(self.axis)[:, 0] * spoke_radius
            for r in range(self.order):
                new_line = line_type(
                    [origin, origin+baseline],
                    color=color,
                    **graphics_options
                )
                objects.append(new_line)
                baseline = tf @ baseline

        return [
            o.plot(figure)
            for o in objects
        ]

class ReflectionElement(SymmetryElement):
    def __init__(self, axis):
        ax, norm = nput.vec_normalize(axis, return_norms=True)
        if norm < 1e-6: raise ValueError("can't have reflection element with no axis")
        self.axis = ax
    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.axis})"

    def get_transformation(self):
        tf = nput.view_matrix(self.axis)
        return tf @ np.diag([1, -1, 1]) @ tf.T
    def inverse(self):
        return self
    def transform(self, tf):
        return type(self)(tf@self.axis)

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

    def plot(self, figure, *,
             origin=None,
             disk_type=None,
             color='black',
             radius=2,
             disk_transparency=.8,
             **graphics_options
             ):
        from McUtils.Plots import Disk
        if origin is None:
            origin = (0, 0, 0)
        origin = np.asarray(origin)
        if disk_type is None:
            disk_type = Disk
        disk = disk_type(
            origin,
            radius=radius,
            normal=self.axis,
            color=color,
            transparency=disk_transparency
        )
        return disk.plot(figure)

class ImproperRotationElement(SymmetryElement):
    def __init__(self, order, axis, root=1):
        ax, norm = nput.vec_normalize(axis, return_norms=True)
        if norm < 1e-6: raise ValueError("can't have improper rotation element with no axis")
        if root > 1 and order % root == 0:
            order = order // root
            root = 1
        self.root = root
        self.order = order
        self.axis = ax
    def inverse(self):
        return type(self)(self.order, self.axis, self.order-self.root)
    def transform(self, tf):
        return type(self)(self.order, tf @ self.axis, root=self.root)

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
                return ImproperRotationElement(order, self.axis, root=int(np.round(root)))
            else:
                return super().compose(other)
        else:
            return other.compose(self) # commute

    def __repr__(self):
        cls = type(self)
        return f"{cls.__name__}({self.root}/{self.order}, {self.axis})"
    def get_transformation(self):
        rot = nput.rotation_matrix(self.axis, 2*np.pi*self.root/self.order)
        tf = nput.view_matrix(self.axis)
        reflect = tf @ np.diag([1, -1, 1]) @ tf.T
        return reflect @ rot

    def plot(self, figure, *,
             origin=None,
             line_type=None,
             disk_type=None,
             color='black',
             line_style='dashed',
             size=2,
             spoke_radius=.3,
             disk_color='black',
             disk_transparency=.8,
             **graphics_options
             ):
        from McUtils.Plots import Line, Disk
        if line_type is None:
            line_type = Line
        if origin is None:
            origin = (0, 0, 0)
        origin = np.asarray(origin)
        points1 = origin - size * self.axis
        points2 = origin + size * self.axis
        objects = []
        axis_line = line_type(
            [points1, points2],
            color=color,
            line_style=line_style,
            **graphics_options
            )
        objects.append(axis_line)
        if spoke_radius is not None:
            if disk_color is not None:
                if disk_type is None:
                    disk_type = Disk
                disk = disk_type(
                    origin,
                    radius=spoke_radius,
                    normal=self.axis,
                    color=disk_color,
                    transparency=disk_transparency
                )
                objects.append(disk)

            tf = self.get_transformation()
            baseline = nput.view_matrix(self.axis)[:, 0] * spoke_radius
            for r in range(self.order):
                new_line = line_type(
                    [origin, origin+baseline],
                    color=color,
                    line_style=line_style,
                    **graphics_options
                )
                objects.append(new_line)
                baseline = tf @ baseline

        return [
            o.plot(figure)
            for o in objects
        ]

# def enumerate_symmetry_operations(generators:'list[SymmetryElement]', max_order=None):
#     if max_order is None:
#         max_order = 60 # support up to C60 rotations...
#
#     real_ops = [g for g in generators if not isinstance(g, IdentityElement)]
#     tf_cache = {}
#     def get_tf(g):
#         tf = tf_cache.get(g)
#         if tf is None:
#             tf = g.get_transformation()
#             tf_cache[g] = tf
#         return tf
#
#     generator_orbits = []
#     for g in real_ops:
#         x = get_tf(g)
#         g2 = g
#         orbit = [g]
#         for i in range(max_order):
#             g3 = g @ g2
#             x2 = g3.get_transformation()
#             if np.allclose(x2, x):
#                 break
#             else:
#                 orbit.append(g3)
#
#     return generator_orbits

