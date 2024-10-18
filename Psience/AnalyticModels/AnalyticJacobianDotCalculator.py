"""
Temporary storage for an attempt to get G matrix elements purely symbolically
but which lacked sufficient information about the coordinates involved
"""

import numpy as np
from .Helpers import sym, AnalyticModelBase

class OrientationVector:
    def __init__(self, x, y, z, basis):
        self.vec = [x, y, z]
        self.basis = basis
    def __neg__(self):
        return type(self)(*(-a for a in self.vec), self.basis)
    def dot(self, other:'OrientationVector'):
        if other.basis != self.basis:
            raise ValueError("bases not the same: {} vs. {}".format(self.basis, other.basis))
        else:
            return sum(x*y for x,y in zip(self.vec, other.vec))
    @classmethod
    def rotation_matrix(cls, angle, axis:'OrientationVector'):
        # a = list(axis.vec)
        # if a == [1, 0, 0] or a == [-1, 0, 0]:
        #     mat = [
        #         OrientationVector(sym.cos(angle), -sym.sin(angle), 0, axis.basis),
        #         OrientationVector(sym.sin(angle),  sym.cos(angle), 0, axis.basis),
        #         OrientationVector(             0,               0, 1, axis.basis)
        #     ]
        #         elif i == 2:
        #             mat = [
        #                 OrientationVector(1,              0,               0, axis.basis),
        #                 OrientationVector(0, sym.cos(angle), -sym.sin(angle), axis.basis),
        #                 OrientationVector(0, sym.sin(angle),  sym.cos(angle), axis.basis)
        #             ]
        #         if v == -axis:
        #             mat = [-a for a in mat]
        # elif a == [1, 0, 0] or a == [-1, 0, 0]:
        #     mat = [
        #             OrientationVector(sym.cos(angle), 0, -sym.sin(angle), axis.basis),
        #             OrientationVector(             0, 1,               0, axis.basis),
        #             OrientationVector(sym.sin(angle), 0,  sym.cos(angle), axis.basis)
        #         ]
        # else:
        #     raise Exception("...")
        a, b, c = axis.vec
        cos = sym.cos(angle)
        ncos = 1-cos
        sin = sym.sin(angle)
        # Base formula: vXv*(1-cos(a) + e3.a * sin(a)
        mat = RotationMatrix(
            OrientationVector(cos + (a**2)*ncos,  a*b*ncos - c*sin,  a*c*ncos + b*sin, axis.basis),
            OrientationVector( a*b*ncos + c*sin, cos + (b**2)*ncos,  b*c*ncos - a*sin, axis.basis),
            OrientationVector( a*c*ncos - b*sin,  b*c*ncos + a*sin, cos + (c**2)*ncos, axis.basis)
        )
        return mat
    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, self.vec, self.basis)
class RotationMatrix:
    def __init__(self, v1:OrientationVector, v2:OrientationVector, v3:OrientationVector):
        self.vectors = [v1, v2, v3]
    def dot(self, v):
        return OrientationVector(*(r.dot(v) for r in self.vectors), self.vectors[0].basis)
    def __neg__(self):
        return type(self)(*(-x for x in self.vectors))
class BondVector:
    from McUtils.Scaffolding import Logger, NullLogger
    # logger = Logger()
    logger = NullLogger()
    def __init__(self, i, j, norm=1, embedding=None):
        if i == j:
            raise ValueError("bond from {} to {} invalid".format(i, j))
        self.i = i
        self.j = j
        self.norm = norm
        self.embedding = embedding

    def one_common_atom(self, bond:'BondVector', mode='shared'):
        if self.i == bond.i or self.i == bond.j:
            if self.j == bond.i or self.j == bond.j:
                raise ValueError("all shared atoms")
            if mode == 'shared':
                return self.i
            else:
                return self.j
        elif self.j == bond.i or self.j == bond.j:
            if mode == 'shared':
                return self.j
            else:
                return self.i
        else:
            raise ValueError("no shared atoms")
    def shared_atom(self, bond:'BondVector'):
        return self.one_common_atom(bond, mode='shared')
    def unshared_atom(self, bond:'BondVector'):
        return self.one_common_atom(bond, mode='unshared')
    @property
    def direction(self):
        return type(self)(self.i, self.j, norm=1, embedding=self.embedding)
    def __repr__(self):
        if self.norm == 1:
            return '{}({}->{})'.format(type(self).__name__, self.i, self.j)
        else:
            return '{} {}({}->{})'.format(self.norm, type(self).__name__, self.i, self.j)
    def __neg__(self):
        return type(self)(self.j, self.i, norm=self.norm, embedding=self.embedding)
    def __mul__(self, other):
        return type(self)(self.i, self.j, norm=self.norm*other, embedding=self.embedding)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __add__(self, other):
        if isinstance(other, BondVectorSum):
            return other + self
        else:
            return BondVectorSum(self, other)
    def __sub__(self, other):
        return self + (-other)
    def __eq__(self, other):
        if not isinstance(other, BondVector):
            return self.dot(other) == 1
        else:
            return other.i == self.i and other.j == self.j
    def as_expr(self):
        return self.norm * (
                AnalyticModelBase.symbolic_x(self.j) - AnalyticModelBase.symbolic_x(self.i)
        ) / AnalyticModelBase.symbolic_r(self.j, self.i)
    @classmethod
    def get_embedding_angle_cos(self, i, j, k, l):
        # l is attached to k by definition
        return -sym.sin(AnalyticModelBase.symbolic_a(j, k, l)) * sym.sin(AnalyticModelBase.symbolic_t(i, j, k, l))
    @classmethod
    def get_positioning_cos(cls, i, j, k, l):
        # l is attached to k by definition
        a_ijk = AnalyticModelBase.symbolic_a(i,j,k)
        a_jkl = AnalyticModelBase.symbolic_a(j,k,l)
        t_ijkl = AnalyticModelBase.symbolic_t(i, j, k, l)
        cos = sym.cos; sin = sym.sin
        return cos(a_ijk)*cos(a_jkl) - sin(a_ijk)*sin(a_jkl)*cos(t_ijkl)
    def angle_cos(self, other, require_int=False):
        if isinstance(other, BondNormal):
            return other.angle_cos(self, require_int=require_int)
        elif isinstance(other, BondVectorSum):
            return other.distribute(lambda x:self.angle_cos(x, require_int=require_int))
        else:
            match_table = (
                self.i == other.i,
                self.i == other.j,
                self.j == other.i,
                self.j == other.j
            )
            match_sum = np.sum(match_table)
            if match_sum == 2: # two shared points
                if match_table[0]:
                    cos = 1
                else:
                    cos = -1
            elif match_sum == 0:
                return self.get_positioning_cos(self.i, self.j, other.i, other.j)
                # i = self.i
                # j = self.j
                # k = other.i
                # l = other.j
                #
                #
                # embedding = self.embedding
                # if embedding is None:
                #     embedding = other.embedding
                # if embedding is None:
                embedding = (self.i, self.j, other.i)
                self.logger.log_print(f"using embedding {embedding}")
                my_embedding = self.polar_components(embedding, self)
                other_embedding = other.polar_components(embedding, other)
                self.logger.log_print(f"emb1: {my_embedding}")
                self.logger.log_print(f"emb2: {other_embedding}")
                cos = my_embedding.dot(other_embedding).expand().simplify()
            else:
                if require_int:
                    return 2 + np.ravel_multi_index(match_table, (2, 2, 2, 2)) # binary encoding
                    # we'll build a composite binary integer so we can keep a single return type
                else:
                    if match_table[0]: # shared first atom
                        i = self.j
                        j = self.i
                        k = other.j
                        sign = 1
                    elif match_table[1]: # shared first and second
                        i = self.j
                        j = self.i
                        k = other.i
                        sign = -1
                    elif match_table[2]: # shared second and first
                        i = self.i
                        j = self.j
                        k = other.j
                        sign = -1
                    else: # shared second atom
                        i = self.i
                        j = self.j
                        k = other.i
                        sign = 1
                    cos = sign * sym.cos(AnalyticModelBase.symbolic_a(i, j, k))
            return cos
    def dot(self, other:'BondVector'):
        if isinstance(other, BondVectorSum):
            return other.distribute(
                lambda x:self.dot(x)
            ).simplify()
        else:
            with self.logger.block():
                self.logger.log_print("{s}", s=self)
                self.logger.log_print("{s}", s=other)
                acos = self.angle_cos(other)
                self.logger.log_print("acos: {acos}", acos=acos)
            return acos * self.norm * other.norm
    def cross(self, other):
        if isinstance(other, BondVectorSum):
            return other.distribute(
                lambda x:self.cross(x)
            )
        else:
            sangle = sym.sin(sym.acos(self.angle_cos(other))).trigsimp().subs(sym.Abs, sym.Id)
            return BondNormal(self, other,
                              norm=self.norm * other.norm * sangle,
                              embedding=self.embedding if self.embedding is not None else other.embedding
                              )

    @classmethod
    def polar_components(self, embedding_atoms, vector):
        raise NotImplementedError("dead code path")
        # returns the components and the basis for the embdedding
        # of atom l with respect to reference atoms i, j, and k
        i, j, k = embedding_atoms#sorted(embedding_atoms)
        basis = BondNormal.embedding_basis((i, j, k))
        # x-axis: i->j
        # z-axis: i->j X j->k
        n = vector.i
        m = vector.j
        sign = 1
        # We have two cases:
        #  1. all stuff inside i,j,k (simple rotations)
        #  2. i->l, j->l, or k->l (rotation + out of plane transform)
        matches = [n in (i, j, k), m in (i, j, k)]
        match_sum = np.sum(matches)

        if match_sum == 2: # simple rotation of main
            if m == j: # pointing towards center, should be flipped
                n, m = m, n
                sign *= -1
            if n == j:
                if m == i: # j->i is a flip of the basis vector
                    sign *= -1
                    angle = None
                else: # j->k, need i->j into j->k
                    angle = AnalyticModelBase.symbolic_a(i, j, k)
            else:
                if n == k: sign *= - 1
                # vector is i->k, need to rotation i->j into i->k
                angle = AnalyticModelBase.symbolic_a(j, i, k)

            if angle is None:
                components = OrientationVector(1, 0, 0, basis)
            else:
                self.logger.log_print("{n}->{m} in ({i},{j},{k}) = {a}", n=n,m=m,i=i,j=j,k=k,a=angle)
                rmat = OrientationVector.rotation_matrix(angle, OrientationVector(0, 0, 1, basis))
                v = OrientationVector(1, 0, 0, basis)
                components = rmat.dot(v)
        elif match_sum == 0:
            self.logger.log_print(f"{n}->{m} in ({i},{j},{k})")
            raise NotImplementedError("this shouldn't happen...")
        else:
            # rotation in plane and rotation out of plane
            if matches[1]: # want first point shared with embedding
                n, m = m, n
                sign *= -1

            if n == i: # simple dihedral
                # m<-i->j->k
                chob = None
                angle = AnalyticModelBase.symbolic_a(m, i, j)
                polar = AnalyticModelBase.symbolic_t(k, j, i, m)
            elif n == j:
                # i->j(->m)->k
                chob = None
                angle = AnalyticModelBase.symbolic_a(i, j, m)
                polar = AnalyticModelBase.symbolic_t(k, i, j, m)
            else:
                sign *= -1
                # need to rotate i->j->k->m so that ij -> kj
                chob = AnalyticModelBase.symbolic_a(i, j, k)
                angle = AnalyticModelBase.symbolic_a(j, k, m)
                polar = AnalyticModelBase.symbolic_t(i, j, k, m)

            self.logger.log_print("{n}->{m} in ({i},{j},{k}) = {a} | {p}", n=n,m=m,i=i,j=j,k=k,a=angle, p=polar)

            ax = OrientationVector(1, 0, 0, basis)
            v = ax
            if chob is not None:
                cr = OrientationVector.rotation_matrix(chob, OrientationVector(0, 0, 1, basis))
                v = cr.dot(v)
            amat = OrientationVector.rotation_matrix(angle, OrientationVector(0, 0, 1, basis))
            v = amat.dot(v)
            tmat = OrientationVector.rotation_matrix(polar, ax)
            v = tmat.dot(v)
            components = v
        self.logger.log_print(f"{n}->{m} in ({i},{j},{k}) -> {sign}")
        if sign < 0:
            components = -components
        return components
class BondNormal:
    def __init__(self, a, b, norm=1, embedding=None):
        self.a = a
        self.b = b
        if abs(self.a.angle_cos(self.b, require_int=True)) == 1:
            raise NotImplementedError('cross product of vector with itself ill defined')
        self.norm = norm
        if embedding is None:
            embedding = a.embedding
        if embedding is None:
            embedding = b.embedding
        if embedding is None:
            if isinstance(a, BondVector):
                i = a.i
                j = a.j
                k = b.i if b.i not in (i, j) else b.j
                embedding = (i, j, k)
        self.embedding = embedding
    @property
    def direction(self):
        return type(self)(self.b, self.a, norm=1, embedding=self.embedding)
    def __repr__(self):
        if self.norm == 1:
            return '{}({}x{})'.format(type(self).__name__, self.a, self.b)
        else:
            return '{} {}({}x{})'.format(self.norm, type(self).__name__, self.a, self.b)
    def __neg__(self):
        return type(self)(self.b, self.a, norm=self.norm, embedding=self.embedding)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        return type(self)(self.a, self.b, norm=self.norm*other, embedding=self.embedding)
    def __add__(self, other):
        if isinstance(other, BondVectorSum):
            return other + self
        else:
            return BondVectorSum(self, other)
    def __sub__(self, other):
        if isinstance(other, BondVectorSum):
            return self + (-other)
        else:
            return BondVectorSum(self, -other)
    def __eq__(self, other):
        if not isinstance(other, BondNormal):
            return self.dot(other) == 1
        else:
            return self.a == other.a and self.b == other.b
    @classmethod
    def embedding_basis(cls, embedding):
        ei = BondVector(embedding[0], embedding[1])
        x = BondVector(embedding[1], embedding[2])
        ej = BondNormal(ei, x)
        ek = BondNormal(ei, ej)
        return [ei, ek, ej]
    @classmethod
    def polar_components(cls, embedding, self):
        raise NotImplementedError("dead code path")
        a = self.a.polar_components(embedding, self.a)
        a1, a2, a3 = a.vec
        b = self.b.polar_components(embedding, self.b)
        b1, b2, b3 = b.vec
        ang = sym.acos(self.a.angle_cos(self.b))
        sangle = sym.sin(ang).trigsimp().subs(sym.Abs, sym.Id)
        return OrientationVector( # gotta normalize...
            (a2*b3 - b2*a3)/sangle,
            (b1*a3 - a1*b3)/sangle,
            (a1*b2 - b1*a2)/sangle,
            a.basis
        )
    def angle_cos(self, other, require_int=True):
        if isinstance(other, BondVectorSum):
            return other.distribute(lambda v:self.angle_cos(v))
        if isinstance(other, BondVector) and isinstance(self.a, BondVector) and isinstance(self.b, BondVector):
            alignment_1 = self.a.angle_cos(other, require_int=True)
            alignment_2 = self.b.angle_cos(other, require_int=True)
            if abs(alignment_1) == 1 or abs(alignment_2) == 1: # cross product is perpendicular to contained vector
                cos = 0
            else:
                i = self.a.unshared_atom(self.b)
                j = self.a.shared_atom(self.b)
                k = self.b.unshared_atom(self.a)

                if other.i == i:
                    l = other.j
                    i,j,k = k,j,i
                    sign = -1
                elif other.j == i:
                    l = other.i
                    i,j,k = k,j,i
                    sign = 1
                elif other.i == j:
                    l = other.j
                    i,j,k = i,k,j
                    sign = -1
                elif other.j == j:
                    l = other.i
                    i,j,k = i,k,j
                    sign = 1
                elif other.i == k:
                    l = other.j
                    sign = 1
                elif other.j == k:
                    l = other.i
                    sign = -1
                else:
                    raise ValueError("5 atom angles unsupported")

                return sign * BondVector.get_embedding_angle_cos(i, j, k, l)

                # embedding = self.embedding
                # # if embedding is None:
                # #     embedding = other.embedding
                # if embedding is not None:
                #     my_embedding = self.polar_components(embedding, self)
                #     other_embedding = other.polar_components(embedding, other)
                #     BondVector.logger.log_print(f"emb1: {my_embedding}")
                #     BondVector.logger.log_print(f"emb2: {other_embedding}")
                #     cos = my_embedding.dot(other_embedding)
                # else:
                #     raise NotImplementedError("{}".format(self))
        else:
            # (axb).(cxd) = (a.c)(b.d) - (a.d)(b.c)
            #             = acos(a,c)acos(b,d) - acos(a,d)acos(b,c)
            aa = self.a.angle_cos(other.a)
            bb = self.b.angle_cos(other.b)
            ab = self.a.angle_cos(other.b)
            ba = self.b.angle_cos(other.a)
            # BondVector.logger.log_print(f"??? aa{aa} bb{bb} ab{ab} ba{ba}")

            ang_ab = sym.acos(self.a.angle_cos(self.b))
            sangle_ab = sym.sin(ang_ab).trigsimp().subs(sym.Abs, sym.Id)
            ang = sym.acos(other.a.angle_cos(other.b))
            sangle_cd = sym.sin(ang).trigsimp().subs(sym.Abs, sym.Id)
            return (aa*bb - ab*ba) / (sangle_ab*sangle_cd)

            alignments = [
                self.a.angle_cos(other.a, require_int=True),
                self.a.angle_cos(other.b, require_int=True),
                self.b.angle_cos(other.a, require_int=True),
                self.b.angle_cos(other.b, require_int=True)
            ]
            matches = [
                abs(x) == 1 for x in alignments
            ]
            match_sum = sum(matches)
            if match_sum == 0:
                raise NotImplementedError('angle between vector cross over different atoms not symbolically useful')
            elif match_sum == 2:
                # just need to determine the sign of alignment...
                if matches[0]: # same vectors basically, but we still need to account for signs
                    cos = alignments[0] * alignments[3]
                else: # flipped ordering so we pick up a negative sign
                    cos = -alignments[1] * alignments[2]
            else:
                # only one of the vectors matched, which means either we have some kind of dihedral angle
                # or we have zero, which we can check by checking to total cross product alignment
                if matches[0] or matches[2]:
                    subalignment = self.angle_cos(other.b, require_int=True)
                else:
                    subalignment = self.angle_cos(other.a, require_int=True)

                if abs(subalignment) == 1: # gotta be perpendicular...
                    cos = 0
                else: # need to pull out components to return the correct dihedral...
                    if matches[0]: # shared first bond, making these atoms j and k
                        i = self.b.unshared_atom(self.a)
                        j = self.b.shared_atom(self.a)
                        k = self.a.unshared_atom(self.b)
                        l = other.b.unshared_atom(other.a)
                        sign = alignments[0]
                    elif matches[1]:  # shared first and second
                        i = self.b.unshared_atom(self.a)
                        j = self.b.shared_atom(self.a)
                        k = self.a.unshared_atom(self.b)
                        l = other.a.unshared_atom(other.b)
                        sign = alignments[1]
                    elif matches[2]:  # shared second and first
                        i = self.a.unshared_atom(self.b)
                        j = self.a.shared_atom(self.b)
                        k = self.b.unshared_atom(self.a)
                        l = other.b.unshared_atom(other.a)
                        sign = alignments[2]
                    else: # shared second atom
                        i = self.a.unshared_atom(self.b)
                        j = self.a.shared_atom(self.b)
                        k = self.b.unshared_atom(self.a)
                        l = other.a.unshared_atom(other.b)
                        sign = alignments[3]
                    cos = sign * sym.cos(AnalyticModelBase.symbolic_t(i, j, k, l))

        return cos
    def as_expr(self):
        i = self.a.unshared_atom(self.b)
        j = self.a.shared_atom(self.b)
        k = self.b.unshared_atom(self.a)
        return self.norm * self.a.norm * self.b.norm * AnalyticModelBase.symbolic_n(k, j, i)
    def cross(self, other):
        if isinstance(other, BondVectorSum):
            return other.distribute(
                lambda x: self.cross(x)
            )
        else:
            sangle = sym.sin(sym.acos(self.angle_cos(other))).trigsimp().subs(sym.Abs, sym.Id)
            return BondNormal(self, other,
                              norm=self.norm * other.norm * sangle,
                              embedding=self.embedding if self.embedding is not None else other.embedding
            )
    def dot(self, other):
        if isinstance(other, BondVectorSum):
            return other.distribute(lambda v:self.dot(v))
        else:
            return self.angle_cos(other) * self.norm * other.norm
class BondVectorSum:
    def __init__(self, *terms):
        self.terms = terms
    def __add__(self, other):
        if isinstance(other, BondVectorSum):
            return type(self)(*self.terms, *other.terms)
        else:
            return type(self)(*self.terms, other)
    def __radd__(self, other):
        return self + other
    def __rmul__(self, other):
        return self.distribute(lambda a:other*a)
    def __mul__(self, other):
        return self.distribute(lambda a:other*a)
    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.terms)
    def __neg__(self):
        return type(self)(*(-t for t in self.terms))
    def simplify(self):
        return sum(self.terms)

    def as_expr(self):
        return sum(
            b.as_expr() if hasattr(b, 'as_expr') else b
            for b in self.terms
        )
    def distribute(self, f):
        all_terms = [f(a) for a in self.terms]
        new_terms = ()
        for t in all_terms:
            new_terms = new_terms + (t.terms if isinstance(t, BondVectorSum) else (t,))
        return type(self)(*new_terms)
    def dot(self, other):
        return self.distribute(lambda a:a.dot(other))

class InternalJacobianDisplacements:

    @staticmethod
    def symbolic_m(i):
        return sym.Symbol('m[{i}]'.format(i=i), positive=True)

    @staticmethod
    def symbolic_r(i, j):
        if j < i:
            i, j = j, i
        return sym.Symbol('r[{i},{j}]'.format(i=i, j=j), positive=True)

    @staticmethod
    def symbolic_a(i, j, k):
        if k < i:
            i, k = k, i
        return sym.Symbol('a[{i},{j},{k}]'.format(i=i, j=j, k=k), real=True)

    @staticmethod
    def symbolic_t(i, j, k, l):
        if l < i:
            i, j, k, l = l, k, j, i
        return sym.Symbol('t[{i},{j},{k},{l}]'.format(i=i, j=j, k=k, l=l), real=True)

    @staticmethod
    def symbolic_y(i, j, k, l):
        return sym.Symbol('y[{i},{j},{k},{l}]'.format(i=i, j=j, k=k, l=l), real=True)

    @classmethod
    def lam(cls, i, j, k):
        r12 = cls.symbolic_r(i, j)
        r23 = cls.symbolic_r(j, k)
        a = cls.symbolic_a(i, j, k)
        return 1 / sym.sin(a) * (1 / r12 - sym.cos(a) / r23)

    @classmethod
    def dr(cls, i, j):
        return [
            BondVector(j, i),
            BondVector(i, j)
        ]

    @classmethod
    def da(cls, i, j, k):
        e12 = BondVector(i, j, embedding=(i, j, k))
        r12 = cls.symbolic_r(i, j)
        e23 = BondVector(j, k, embedding=(i, j, k))
        r23 = cls.symbolic_r(j, k)
        a = cls.symbolic_a(i, j, k)
        l321 = cls.lam(k, j, i)
        l123 = cls.lam(i, j, k)
        return [
            -1 / (r12 * sym.sin(a)) * (sym.cos(a) * e12 + e23),
            -l321 * e12 + l123 * e23,
            1 / (r23 * sym.sin(a)) * (e12 + sym.cos(a) * e23)
        ]

    @classmethod
    def dt(cls, i, j, k, l):
        e12 = BondVector(i, j, embedding=(i, j, k))
        r12 = cls.symbolic_r(i, j)
        e23 = BondVector(j, k, embedding=(i, j, k))
        r23 = cls.symbolic_r(j, k)
        e34 = BondVector(k, l, embedding=(i, j, k))
        r34 = cls.symbolic_r(k, l)
        a1 = cls.symbolic_a(i, j, k)
        a2 = cls.symbolic_a(j, k, l)
        l123 = cls.lam(i, j, k)
        l432 = cls.lam(l, k, j)
        x123 = BondNormal(e12, e23, embedding=(i, j, k))
        x234 = BondNormal(e23, e34, embedding=(j, k, l))
        return [
           -1 / (r12 * sym.sin(a1)) * x123,
            l123 * x123 - sym.cot(a2) / r23 * x234,
           -l432 * x234 + sym.cot(a1) / r23 * x123,
            1 / (r34 * sym.sin(a2)) * x234
        ]

    @classmethod
    def dy(cls, i, j, k, l):
        raise NotImplementedError("oops")

    @classmethod
    def displacement_vectors(cls, inds, coord_type=None):
        if coord_type is None:
            if len(inds) == 2:
                return cls.dr(*inds)
            elif len(inds) == 3:
                return cls.da(*inds)
            elif len(inds) == 4:
                return cls.dt(*inds)
            else:
                return cls.dy(*inds)
        else:
            return {
                "r":cls.dr,
                "a":cls.da,
                "t":cls.dt,
                "y":cls.dy
            }[coord_type](*inds)