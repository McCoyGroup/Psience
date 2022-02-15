
import numpy as np

__all__ = [
    'AnalyticPotentialConstructor',
    'AnalyticGMatrixConstructor',
    "AnalyticPseudopotentialConstructor"
]

try:
    import sympy as sym
except ImportError:
    class SympyNotLoadedError(ImportError):
        ...
    class SympyNotLoadedStub:
        def __getattr__(self, item):
            raise SympyNotLoadedError("Analytic models depend on sympy, which couldn't be imported")
    sym = SympyNotLoadedStub()

class AnaylticModelBase:
    @classmethod
    def take_derivs(cls, expr, vars):
        if isinstance(expr, list):
            return [cls.take_derivs(e, vars) for e in expr]
        else:
            return [sym.diff(expr, v) for v in vars]

    @classmethod
    def eval_exprs(cls, expr, subs):
        if isinstance(expr, list):
            return [cls.eval_exprs(e, subs) for e in expr]
        elif isinstance(expr, (int, float, np.integer, np.floating)):
            return expr
        else:
            return expr.subs(subs)


    @classmethod
    def symbol_list(cls, names, instance=None):
        instance_tag = "" if instance is None else str(instance)
        return sym.symbols(" ".join(n+instance_tag for n in names))


class AnalyticPotentialConstructor(AnaylticModelBase):
    @classmethod
    def symbolic_morse(cls, instance=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.morse(*cls.symbol_list(["De", "a", "r", "re"], instance=instance))
    @staticmethod
    def morse(De, a, r, re):
        """
        :param De:
        :type De:
        :param a:
        :type a:
        :param r:
        :type r:
        :return:
        :rtype: sym.expr
        """
        return De * (1 - sym.exp(-a * (r-re)))**2
    @staticmethod
    def harm(k, x, x_e):
        """

        :param De:
        :type De:
        :param a:
        :type a:
        :param r:
        :type r:
        :return:
        :rtype: sym.expr
        """
        return k*(x-x_e)**2
    @classmethod
    def symbolic_harmonic(cls, instance=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.harm(*cls.symbol_list(["k", "q", "qe"], instance=instance))
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
        for i,v in enumerate(axis.vec):
            if v == axis or v == -axis:
                if i == 0:
                    mat = [
                        OrientationVector(sym.cos(angle), -sym.sin(angle), 0, axis.basis),
                        OrientationVector(sym.sin(angle),  sym.cos(angle), 0, axis.basis),
                        OrientationVector(             0,               0, 1, axis.basis)
                    ]
                elif i == 1:
                    mat = [
                        OrientationVector(sym.cos(angle), 0, -sym.sin(angle), axis.basis),
                        OrientationVector(             0, 1,               0, axis.basis),
                        OrientationVector(sym.sin(angle), 0,  sym.cos(angle), axis.basis)
                    ]
                elif i == 2:
                    mat = [
                        OrientationVector(1,              0,               0, axis.basis),
                        OrientationVector(0, sym.cos(angle), -sym.sin(angle), axis.basis),
                        OrientationVector(0, sym.sin(angle),  sym.cos(angle), axis.basis)
                    ]
                if v == -axis:
                    mat = [-a for a in mat]
        else:
            a, b, c = axis.vec
            cos = sym.cos(angle)
            ncos = 1-cos
            sin = sym.sin(angle)
            # Base formula: vXv*(1-cos(a) + e3.a * sin(a)
            mat = RotationMatrix(
                OrientationVector(cos + (a**2)*ncos,  a*b*ncos - c*sin,  a*c*ncos + b*sin, axis.basis),
                OrientationVector( a*b*ncos - c*sin, cos + (b**2)*ncos,  b*c*ncos - a*sin, axis.basis),
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
    def __init__(self, i, j, norm=1, embedding=None):
        if i == j:
            raise ValueError("bond from {} to {} invalid".format(i, j))
        self.i = i
        self.j = j
        self.norm = norm
        self.embedding = embedding
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
        return BondVectorSum(self, other)
    def __eq__(self, other):
        return self.dot(other) == 1
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
            else: # one shared point

                if require_int:
                    return 2 + np.ravel_multi_index(match_table, (2, 2, 2, 2)) # binary encoding
                else:
                    embedding = self.embedding
                    if embedding is None:
                        embedding = other.embedding
                    if embedding is not None:
                        my_embedding = self.polar_components(embedding, self)
                        other_embedding = other.polar_components(embedding, other)
                        return my_embedding.dot(other_embedding)
                    else:
                        if match_sum == 0:
                            raise NotImplementedError("...")
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
                            cos = sign * sym.cos(AnalyticKineticTerm.symbolic_a(i, j, k))
                # we'll build a composite binary integer so we can keep a single return type
            return cos
    def dot(self, other:'BondVector'):
        if isinstance(other, BondVectorSum):
            return other.distribute(
                lambda x:self.angle_cos(x) * self.norm * x.norm
            ).simplify()
        else:
            return self.angle_cos(other) * self.norm * other.norm
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
        # returns the components and the basis for the embdedding
        # of atom l with respect to reference atoms i, j, and k
        i, j, k = embedding_atoms
        # embedding = [
        #         BondVector(k, j),
        #         BondVector(j, k).cross(BondVector(i, j)),
        #         BondVector(j, k).cross(BondVector(i, j)).cross(BondVector(j, k))
        #     ]
        basis = [
            BondVector(k, j),
            BondNormal(BondVector(k, j), BondVector(i, j)),
            BondNormal(BondNormal(BondVector(k, j), BondVector(i, j)), BondVector(j, k))
            ]
        kk = vector.i
        l = vector.j
        sign = 1
        # given that our x is k->j we have two cases:
        #  1. all stuff inside i,j,k (simple rotations)
        #  2. i->l, j->l, or k->l (rotation + out of plane transform)
        if ( # simple canonicalizations to cut down on the number of cases
                l == k or
                kk == i or
                kk not in (i, j, k)
        ): # flip the vector and multiply all components by negative one
            kk, l = l, kk
            sign = -1
        matches = [l in (i, j, k), kk in (i, j, k)]
        match_sum = np.sum(matches)

        if match_sum == 2: # simple rotation of main
            if kk == k: # vector is k->i
                if l == i:
                    angle = AnalyticKineticTerm.symbolic_a(j, k, i)
                else:  # k->j
                    angle = None
            else:  # l == i and  k == j
                angle = sym.pi - AnalyticKineticTerm.symbolic_a(i, j, k)
            if angle is None:
                components = OrientationVector(1, 0, 0, basis)
            else:
                rmat = OrientationVector.rotation_matrix(angle, OrientationVector(0, 0, 1, basis))
                v = OrientationVector(1, 0, 0, basis)
                components = rmat.dot(v)
        elif match_sum == 0:
            raise NotImplementedError("this shouldn't happen...")
        else: # rotation in plane and rotation out of plane
            if kk == i:
                asign = -1
                chob = sym.pi - AnalyticKineticTerm.symbolic_a(i, j, k)
                angle = AnalyticKineticTerm.symbolic_a(j, i, l)
                polar = AnalyticKineticTerm.symbolic_t(k, j, i, l)
            elif kk == k:
                asign = 1
                chob = None
                angle = AnalyticKineticTerm.symbolic_a(j, k, l)
                polar = AnalyticKineticTerm.symbolic_t(i, j, k, l)
            elif kk == j: # this feels wrong...
                raise Exception(embedding_atoms, (i, j, k))
                asign = -1
                chob = None
                angle = AnalyticKineticTerm.symbolic_a(j, k, l)
                polar = AnalyticKineticTerm.symbolic_t(i, k, j, l)
            else:
                raise NotImplementedError("case: {}->{} in {}".format(kk, l, embedding_atoms))

            v = OrientationVector(1, 0, 0, basis)
            # if chob is not None:
            #     cr = OrientationVector.rotation_matrix(chob, OrientationVector(0, 0, 1, basis))
            #     v = cr.dot(v)
            amat = OrientationVector.rotation_matrix(angle, OrientationVector(0, 0, 1, basis))
            if asign == -1:
                amat = -amat
            v = amat.dot(v)
            tmat = OrientationVector.rotation_matrix(polar, OrientationVector(1, 0, 0, basis))
            v = tmat.dot(v)
            if chob is not None:
                cr = OrientationVector.rotation_matrix(-chob, OrientationVector(0, 0, 1, basis))
                v = cr.dot(v)
            components = v
            # raise Exception((kk, l), embedding_atoms, components)
        #
        # if match_sum == 2: # totally embedded, meaning simple rotation
        #     if kk == k:
        #         if l == i:
        #             components = [
        #                 sym.cos(AnalyticKineticTerm.symbolic_a(j, k, i)),
        #                 sym.sin(AnalyticKineticTerm.symbolic_a(j, k, i)),
        #                 0
        #             ]
        #         else: # l == j
        #             components = [
        #                 1,
        #                 0,
        #                 0
        #             ]
        #     else: # l == i and  k == j
        #         components = [
        #             -sym.cos(AnalyticKineticTerm.symbolic_a(i, j, k)),
        #             -sym.sin(AnalyticKineticTerm.symbolic_a(i, j, k)),
        #             0
        #         ]
        # elif match_sum == 0: # totally detached
        #     raise NotImplementedError("this shouldn't happen...")
        # else: # only one component is in the group
        #     if kk == i: # means l not in the group
        #         ...
        #     elif kk == j:
        #         ...
        #     elif kk == k:
        #         components = [
        #             sym.cos(AnalyticKineticTerm.symbolic_a(j, k, l)),
        #             sym.cos(AnalyticKineticTerm.symbolic_t(i, j, k, l))
        #               * sym.sin(AnalyticKineticTerm.symbolic_a(j, k, l)),
        #             sym.sin(AnalyticKineticTerm.symbolic_t(i, j, k, l))
        #         ]
        #     else:
        #         raise NotImplementedError("this shouldn't happen... ({} and {} in {})".format(kk, l, (i, j, k)))
        #
        #
        # #
        # #     if matches[0] and matches[1]:
        # #         ...
        # # if kk == k:
        # #     if l not in (i, j, k): #ordering i, j, k, l
        # #         components = [
        # #             sym.cos(AnalyticKineticTerm.symbolic_a(j, k, l)),
        # #             sym.cos(AnalyticKineticTerm.symbolic_t(i, j, k, l)) * sym.sin(AnalyticKineticTerm.symbolic_a(j, k, l)),
        # #             sym.sin(AnalyticKineticTerm.symbolic_t(i, j, k, l))
        # #         ]
        # #     elif l == i:
        # #         components = [
        # #             sym.cos(AnalyticKineticTerm.symbolic_a(j, k, i)),
        # #             sym.sin(AnalyticKineticTerm.symbolic_a(j, k, i)),
        # #             0
        # #         ]
        # #     elif l == j:
        # #         components = [
        # #             1,
        # #             0,
        # #             0
        # #         ]
        # #     else:
        # #         raise NotImplementedError("this shouldn't happen...")
        # # elif l == i:
        # #     if kk not in (i, j):  # ordering kk, i, j, k, l
        # #         components = [
        # #             sym.cos(AnalyticKineticTerm.symbolic_a(j, k, l)),
        # #             sym.cos(AnalyticKineticTerm.symbolic_t(i, j, k, l)) * sym.sin(
        # #                 AnalyticKineticTerm.symbolic_a(j, k, l)),
        # #             sym.sin(AnalyticKineticTerm.symbolic_t(i, j, k, l))
        # #         ]
        # #     elif l == i:
        # #         components = [
        # #             sym.cos(AnalyticKineticTerm.symbolic_a(j, k, i)),
        # #             sym.sin(AnalyticKineticTerm.symbolic_a(j, k, i)),
        # #             0
        # #         ]
        # #     elif l == j:
        # #         components = [
        # #             1,
        # #             0,
        # #             0
        # #         ]
        # #     else:
        # #         raise NotImplementedError("this shouldn't happen...")
        # # else:
        # #     # need to consider all cases where only one of the embedding atoms is anchored
        # #     matches = [l in (i, j), kk in (i, j)]
        # #     match_sum = np.sum(matches)
        # #     elif matches[0]:
        # #         raise Exception(l, kk, i, j, k)
        # #         if l == i:
        # #             ...
        # #         else:
        # #             ...
        # #     elif matches[1]:
        # #         raise NotImplementedError("we'll cover this when we get here...")
        # #     else:
        # #         raise NotImplementedError("this shouldn't happen...")
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
        return BondVectorSum(self, other)
    def __sub__(self, other):
        return BondVectorSum(self, -other)
    def __eq__(self, other):
        return self.dot(other) == 1
    @classmethod
    def polar_components(cls, embedding, self):
        a = self.a.polar_components(embedding, self.a)
        a1, a2, a3 = a.vec
        b = self.b.polar_components(embedding, self.b)
        b1, b2, b3 = b.vec
        sangle = sym.sin(sym.acos(self.a.angle_cos(self.b))).trigsimp().subs(sym.Abs, sym.Id)
        return OrientationVector( # gotta normalize...
            (a2*b3 - b2*a3)/sangle,
            (b1*a3 - a1*b3)/sangle,
            (a1*b2 - b1*a2)/sangle,
            a.basis
        )
    def angle_cos(self, other, require_int=True):
        if isinstance(other, BondVector):
            alignment_1 = self.a.angle_cos(other, require_int=True)
            alignment_2 = self.b.angle_cos(other, require_int=True)
            if abs(alignment_1) == 1 or abs(alignment_2) == 1: # cross product is perpendicular to contained vector
                return 0
            else:
                embedding = self.embedding
                if embedding is None:
                    embedding = other.embedding
                if embedding is not None:
                    my_embedding = self.polar_components(embedding, self)
                    other_embedding = other.polar_components(embedding, other)
                    return my_embedding.dot(other_embedding)
                else:
                    raise NotImplementedError("{}".format(self))
                # raise NotImplementedError('angle between vector and vector cross not symbolically useful ({} and {})'.format(
                #     self,
                #     other
                # ))
        else:
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

                if abs(subalignment) == 1:
                    cos = 0
                else: # need to pull out components to return the correct dihedral...
                    if matches[0]: # shared first bond, making these atoms j and k
                        i = self.b.unshared_atom(self.a)
                        j = self.b.shared_atom(self.a)
                        k = self.a.unshared_atom(self.b)
                        l = other.b.unshared_atom(other.a)
                    elif matches[1]:  # shared first and second
                        i = self.b.unshared_atom(self.a)
                        j = self.b.shared_atom(self.a)
                        k = self.a.unshared_atom(self.b)
                        l = other.a.unshared_atom(other.b)
                    elif matches[2]:  # shared second and first
                        i = self.a.unshared_atom(self.b)
                        j = self.a.shared_atom(self.b)
                        k = self.b.unshared_atom(self.a)
                        l = other.b.unshared_atom(other.a)
                    else: # shared second atom
                        i = self.a.unshared_atom(self.b)
                        j = self.a.shared_atom(self.b)
                        k = self.b.unshared_atom(self.a)
                        l = other.a.unshared_atom(other.b)
                    cos = sym.cos(AnalyticKineticTerm.symbolic_t(i, j, k, l))

            return cos
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
        return self.angle_cos(other) * self.norm * other.norm
class BondVectorSum:
    def __init__(self, *terms):
        self.terms = terms
    def __add__(self, other):
        if isinstance(other, BondVectorSum):
            return type(self)(*self.terms, *other.terms)
        else:
            return type(self)(*self.terms, other)
    def __rmul__(self, other):
        return self.distribute(lambda a:other*a)
    def __mul__(self, other):
        return self.distribute(lambda a:other*a)
    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.terms)
    def simplify(self):
        return sum(self.terms)
    def distribute(self, f):
        return type(self)(*(f(a) for a in self.terms))
    def dot(self, other):
        return self.distribute(lambda a:a.dot(other))

class AnalyticKineticTerm(AnaylticModelBase):
    """
    Provides the core components needed to evaluate Wilson G-matrix element terms
    using the approach of Frederick and Woywood
    """

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
        return 1/sym.sin(a)*(1/r12 - sym.cos(a)/r23)
    @classmethod
    def dr(cls, i, j):
        return [
            -BondVector(i, j),
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
        # raise Exception((-l1* + e12.dot(l2*e23)).simplify())
        return [
            1/(r12*sym.sin(a)) * (sym.cos(a)*e12 + e23),
            -l321*e12 + l123*e23,
            1/(r23*sym.sin(a)) * (e12 + sym.cos(a)*e23)
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
        l234 = cls.lam(l, k, j)
        x123 = BondNormal(e12, e23)
        x234 = BondNormal(e23, e34)
        # raise Exception(l123)
        return [-1*x for x in [
            -1/(r12*sym.sin(a1)) * x123,
             l123*x123 - sym.cot(a2)/r23 * x234,
            -l234*x234 + sym.cot(a1)/r23 * x123,
             1/(r34*sym.sin(a2)) * x234
        ]]
    @classmethod
    def dy(cls, i, j, k, l):
        raise NotImplementedError("oops")

    @classmethod
    def displacement_vectors(cls, inds):
        if len(inds) == 2:
            return cls.dr(*inds)
        elif len(inds) == 3:
            return cls.da(*inds)
        elif len(inds) == 4:
            return cls.dt(*inds)
        else:
            return cls.dy(*inds)

class AnalyticGMatrixConstructor(AnalyticKineticTerm):
    @classmethod
    def g(cls, inds1:'Iterable[int]', inds2:'Iterable[int]'):
        inds1 = list(inds1)
        inds2 = list(inds2)
        if len(inds1) < len(inds2):
            iter_inds = inds1
            rem_inds = inds2
        else:
            iter_inds = inds2
            rem_inds = inds1
        d1 = cls.displacement_vectors(iter_inds)
        d2 = cls.displacement_vectors(rem_inds)
        g_contribs = []
        for n,i in enumerate(iter_inds):
            try:
                m = rem_inds.index(i)
            except ValueError:
                pass
            else:
                # raise Exception(d1[n], d2[m])
                g_contribs.append(1/cls.symbolic_m(i)*d1[n].dot(d2[m]))
        g = sum(g_contribs)
        if isinstance(g, sym.Expr):
            g = g.expand().simplify().expand()#.subs(sym.Abs, sym.Id).simplify()
        return g

    @staticmethod
    def grr(m1, m2):
        return (1 / m1 + 1 / m2)
    @classmethod
    def symbolic_grr(cls, instance=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.grr(*cls.symbol_list(["m1", "m2"], instance=instance))
    @staticmethod
    def gr1r2(m1, q):
        return 1/m1 * sym.cos(q)
    @staticmethod
    def grq(m1, r, q):
        return -1/(m1*r)*sym.sin(q)
    @staticmethod
    def gqq(m1, m2, m3, r1, r2, q):
        return 1/(m2*r1**2) + 1/(m3*r2**2) + 1/(m1*r1**2) + 1/(m1*r2**2) - 2*sym.cos(q)/(m1 * r1 * r2)

class AnalyticPseudopotentialConstructor:
    @classmethod
    def u_rr(cls, m1, m2):
        return 0
    @classmethod
    def u_r1r2(cls, m, q, r1, r2):
        return AnalyticGMatrixConstructor.gr1r2(m, q) / (r1 * r2)
    @classmethod
    def u_rq(cls, m, r1, r2, q):
        return -sym.cos(q) / (m * r1 * r2)
    @classmethod
    def u_qq(cls, m1, m2, m3, r1, r2, q):
        return sym.cos(q) / (2 * m1 * r1 * r2) - 1 / 8 * AnalyticGMatrixConstructor.gqq(m1, m2, m3, r1, r2, q) * (2 + sym.cot(q) ** 2)