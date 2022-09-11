
import numpy as np, scipy
import sympy as sym

from McUtils.Data import AtomData, UnitsData
from .Helpers import AnalyticModelBase
from ..Data import KEData

__all__ = [
    'AnalyticPotentialConstructor',
    'AnalyticKineticEnergyConstructor',
    'AnalyticModel'
]
__reload_hook__ = ["..Data"]

class AnalyticPotentialConstructor(AnalyticModelBase):
    """
    Provides a set of symbolic potentials for use in models

    :related:AnalyticModel, AnalyticKineticEnergyConstructor
    """
    @classmethod
    def morse(cls, *args, De=None, a=None, re=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.calc_morse(
            AnalyticModelBase.symbol("De", *args) if De is None else De,
            AnalyticModelBase.symbol("ap", *args) if a is None else a,
            AnalyticModelBase.symbolic_r(*args),
            AnalyticModelBase.symbol("re", *args) if re is None else re
        )
    @staticmethod
    def calc_morse(De, a, r, re):
        return De * (1 - AnalyticModelBase.sym.exp(-a * (r-re)))**2

    @staticmethod
    def harm(k, x, x_e):
        return k*(x-x_e)**2
    @classmethod
    def harmonic(cls, *args, k=None, qe=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.harm(
            AnalyticModelBase.symbol("k", *args) if k is None else k,
            AnalyticModelBase.var(*args),
            AnalyticModelBase.symbol("qe", *args) if qe is None else qe,
        )

    @staticmethod
    def lin(k, x, x_e):
        return k * (x - x_e)
    @classmethod
    def linear(cls, *args, k=1, xe=None):
        """
        Returns a fully symbolic form of a linear function
        :return:
        :rtype:
        """
        return cls.lin(
            AnalyticModelBase.symbol("k", *args) if k is None else k,
            AnalyticModelBase.var(*args),
            AnalyticModelBase.symbol("xe", *args) if xe is None else xe
        )

    @classmethod
    def multiwell(cls, *args, turning_points=None, origin=0, minimum=0, depth=None):
        """

        :param args:
        :type args:
        :param turning_points:
        :type turning_points:
        :param depth:
        :type depth:
        :return:
        :rtype:
        """
        if turning_points is None:
            raise ValueError("need turning points for multiwell")
        turning_points = [origin+t for t in turning_points]
        x = AnalyticModelBase.var(*args)
        D = AnalyticModelBase.symbol("D", *args) if depth is None else depth
        terms = [x-t for t in turning_points]
        integrand = terms[0]
        for t in terms[1:]:
            integrand *= t
        poly = sym.integrate(integrand, x)
        min_val = np.inf
        min_pos = None
        turning_vals = [poly.subs([[x, t]]) for t in turning_points]
        for n, v in enumerate(turning_vals):
            if v < min_val:
                min_val = v
                min_pos = n
        if min_pos == 0:
            barrier_height = turning_vals[min_pos+1] - turning_vals[min_pos]
        elif min_pos == len(turning_vals) - 1:
            barrier_height = turning_vals[min_pos-1] - turning_vals[min_pos]
        else:
            barrier_height = min(turning_vals[min_pos+1], turning_vals[min_pos-1]) - turning_vals[min_pos]
        return D/barrier_height * (poly - min_val) + minimum

class AnalyticKineticEnergyConstructor(AnalyticModelBase):
    """
    Provides G and V' elements from Frederick and Woywood

    :related:AnalyticModel, AnalyticPotentialConstructor
    """

    @classmethod
    def _get_coord_key(cls, inds1:'Iterable[int]', inds2:'Iterable[int]', coord_types=None):
        # we infer which type of element to use and return that
        # with appropriate substitutions
        if coord_types is None:
            coord_types = [None, None]
        coord_types = list(coord_types)
        if coord_types[0] is None:
            coord_types[0] = cls.infer_coord_type(inds1)
        if coord_types[1] is None:
            coord_types[1] = cls.infer_coord_type(inds2)
        type1, type2 = coord_types
        sorting = ['r', 'a', 't', 'y']
        if sorting.index(coord_types[0]) > sorting.index(coord_types[1]):
            inds1, inds2 = inds2, inds1
            type1, type2 = type2, type1
        shared_indices = np.intersect1d(inds1, inds2)
        mapping = {i:k for i,k in zip(shared_indices, tuple(range(1, len(shared_indices) + 1)))}
        n = len(shared_indices) + 1
        diff_inds = np.concatenate([np.setdiff1d(inds1, shared_indices),np.setdiff1d(inds2, shared_indices)])
        for i in diff_inds:
            mapping[i] = n
            n += 1

        return ((type1, type2), tuple(mapping[i] for i in inds1), tuple(mapping[i] for i in inds2)), mapping

    @classmethod
    def g(cls, inds1:'Iterable[int]', inds2:'Iterable[int]', coord_types=None, target_symbols=None):
        key, mapping = cls._get_coord_key(inds1, inds2, coord_types=coord_types)
        (expr, _), perm = KEData.find_expressions(key, return_permutation=True)
        # print(expr, perm, mapping)
        if perm is not None:
            shared_indices = np.intersect1d(inds1, inds2)
            p1, p2 = perm
            # print(perm, mapping, inds1)
            if p1 is not None:
                updates = {}
                for i,p in enumerate(p1):
                    old = inds1[i]
                    new = inds1[p]
                    if (
                            old not in shared_indices and new not in shared_indices
                            or old in shared_indices and new in shared_indices
                    ):
                        updates[new] = mapping[old]
                mapping.update(updates)
            if p2 is not None:
                updates = {}
                for i,p in enumerate(p2):
                    old = inds2[i]
                    new = inds2[p]
                    if (
                            old not in shared_indices and new not in shared_indices
                            or old in shared_indices and new in shared_indices
                    ):
                        updates[new] = mapping[old]
                mapping.update(updates)
            # print(perm)
            # print(mapping)
            # print(expr)

        if not isinstance(expr, (int, np.integer)):
            rev_mapping = {v:k for k,v in mapping.items()}
            subs = tuple((s, cls.reindex_symbol(s, rev_mapping, target_symbols=target_symbols)) for s in expr.free_symbols)
            expr = expr.subs([(s, q) for s,(q, _) in subs])
            expr = expr.subs([(q, f) for s,(q, f) in subs])
        return expr

    @classmethod
    def vp(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None):
        key, mapping = cls._get_coord_key(inds1, inds2, coord_types=coord_types)
        (_, expr), perm = KEData.find_expressions(key, return_permutation=True)
        # print(expr, perm, mapping)
        if perm is not None:
            shared_indices = np.intersect1d(inds1, inds2)
            p1, p2 = perm
            # print(perm, mapping, inds1)
            if p1 is not None:
                updates = {}
                for i, p in enumerate(p1):
                    old = inds1[i]
                    new = inds1[p]
                    if (
                            old not in shared_indices and new not in shared_indices
                            or old in shared_indices and new in shared_indices
                    ):
                        updates[new] = mapping[old]
                mapping.update(updates)
            if p2 is not None:
                updates = {}
                for i, p in enumerate(p2):
                    old = inds2[i]
                    new = inds2[p]
                    if (
                            old not in shared_indices and new not in shared_indices
                            or old in shared_indices and new in shared_indices
                    ):
                        updates[new] = mapping[old]
                mapping.update(updates)
            # print(perm)
            # print(mapping)
            # print(expr)

        if not isinstance(expr, (int, np.integer)):
            rev_mapping = {v: k for k, v in mapping.items()}
            subs = tuple(
                (s, cls.reindex_symbol(s, rev_mapping, target_symbols=target_symbols)) for s in expr.free_symbols)
            expr = expr.subs([(s, q) for s, (q, _) in subs])
            expr = expr.subs([(q, f) for s, (q, f) in subs])
        return expr

    @classmethod
    def infer_coord_type(cls, inds):
        if len(inds) == 2:
            return "r"
        elif len(inds) == 3:
            return "a"
        elif len(inds) == 4:
            return "t"
        else:
            raise ValueError("too many indices")

    @classmethod
    def _g(cls, inds1:'Iterable[int]', inds2:'Iterable[int]'):
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
                a = d1[n]
                b = d2[m]
                # print(a, b)
                g_contribs.append(1/cls.symbolic_m(i)*a.dot(b))

        g = sum(g_contribs)
        if isinstance(g, AnalyticModelBase.sym.Expr):
            g = g.expand().simplify().expand()#.subs(sym.Abs, sym.Id).simplify()
        return g

class AnalyticModel:
    """
    Provides a symbolic representation of an analytically evaluatable Hamiltonian
    which can be used to get derived expressions to evaluate.
    Supplies methods to automatically run DVR and VPT calculations from the model
    specifications as well.

    :related:AnalyticPotentialConstructor, AnalyticPotentialConstructor
    """

    def __init__(self, coordinates, potential, dipole=None, values=None, rotation=None):
        self.coords = coordinates
        self.vals = values
        self._syms = None
        self._base_g = None
        self._g = None
        self._u = None
        self.pot = potential
        self.dip = dipole
        if rotation is not None:
            if len(rotation) == 2:
                r, i = rotation
                r = np.asanyarray(r)
                i = np.asanyarray(i)
                if r.ndim == 2:
                    rotation = r
                    inverse = i
                else:
                    inverse = np.linalg.inv(rotation)
            else:
                inverse = np.linalg.inv(rotation)
        else:
            inverse = None
        self.rotation = rotation
        self.inverse = inverse
        self._rot_coords = None

    @property
    def internal_coordinates(self):
        if self._syms is None:
            self._load_symbols()
        return self._syms[0]
    @property
    def constants(self):
        return {k:self.vals[k] for k in self.vals.keys() - set(self.internal_coordinates)}

    def normal_modes(self, dimensionless=True):
        v = self.v(evaluate=True)[-1]
        g = self.g(evaluate=True)[-1]
        freqs2, mode_inv = scipy.linalg.eigh(v, g, type=2) # modes come out as a (internals, normal_mode) array
        freqs = np.sqrt(freqs2)

        # normalization = np.broadcast_to(1 / np.linalg.norm(modes, axis=0), modes.shape)
        # mode_inv = mode_inv  # now as (normal_mode, internals) with frequency dimension removed
        modes = np.linalg.inv(mode_inv)
        if dimensionless:
            mode_inv = mode_inv * np.sqrt(freqs[np.newaxis, :])
            modes = modes / np.sqrt(freqs[:, np.newaxis])
        else:
            modes = modes / np.linalg.norm(modes, axis=1)[:, np.newaxis]
            mode_inv = mode_inv / np.linalg.norm(mode_inv, axis=0)[np.newaxis, :]


        # coords = [AnalyticModelBase.dot(v, self.coords) for v in modes]

        return freqs, modes, mode_inv

    def to_normal_modes(self, dimensionless=True):
        import copy

        if self.rotation is not None:
            raise ValueError("already rotated")
        else:
            freqs, modes, mode_inv = self.normal_modes(dimensionless=dimensionless)
            new = copy.copy(self)
            new._base_g = new._g = new._u = None
            new.rotation = modes
            new.inverse = mode_inv
            return new, freqs

    def get_VPT_expansions(self,
                           order=2,
                           expansion_order=None,
                           include_potential=None,
                           include_gmatrix=None,
                           include_pseudopotential=None,
                           evaluate=True
                           ):
        if expansion_order is None:
            expansion_order = {}
        if 'potential' in expansion_order:
            pot_order = expansion_order['potential']
        elif 'default' in expansion_order:
            pot_order = expansion_order['default']
        else:
            pot_order = order
        if 'kinetic' in expansion_order:
            ke_order = expansion_order['kinetic']
        elif 'default' in expansion_order:
            ke_order = expansion_order['default']
        else:
            ke_order = order
        if 'pseudopotential' in expansion_order:
            u_order = expansion_order['pseudopotential']
        elif 'default' in expansion_order:
            u_order = expansion_order['default']
        else:
            u_order = order

        terms = {}
        if include_gmatrix is not False:
            terms['kinetic_terms'] = self.g(ke_order, evaluate=evaluate)
        if include_potential is not False:
            terms['potential_terms'] = self.v(pot_order+2, evaluate=evaluate)[2:]
        if include_pseudopotential is not False:
            terms['pseudopotential_terms'] = self.vp(u_order, evaluate=evaluate)


        if 'dipole' in expansion_order:
            d_order = expansion_order['dipole']
        elif 'default' in expansion_order:
            d_order = expansion_order['default']
        else:
            d_order = order
        if self.dip is not None:
            terms['dipole_terms'] = self.mu(d_order+1, evaluate=evaluate)
        return terms
    def run_VPT(self,
                order=2,
                states=2,
                return_analyzer=True,
                expansion_order=None,
                include_potential=None,
                include_gmatrix=None,
                include_pseudopotential=None,
                **kwargs
    ):
        from ..VPT2 import VPTAnalyzer, VPTRunner

        if self.rotation is None:
            self, _ = self.to_normal_modes()
        expansions = self.get_VPT_expansions(order=order,
                                             expansion_order=expansion_order,
                                             include_potential=include_potential,
                                             include_gmatrix=include_gmatrix,
                                             include_pseudopotential=include_pseudopotential,
                                             evaluate=True
                                             )
        freqs = np.diag(expansions['potential_terms'][0])
        ncoords = len(freqs)
        natoms = (ncoords + 6) // 3
        if return_analyzer:
            runner = VPTAnalyzer.run_VPT
        else:
            runner = VPTRunner.run_simple
        analyzer = runner(
            [["H"] * natoms, np.zeros((natoms, 3))],  # dummy data since this won't be used at all
            states,
            modes={
                "freqs": freqs,
                "matrix": np.zeros((3 * natoms, len(freqs)))  # dummy data since this won't be used at all
            },
            order=order,
            **dict(expansions, **kwargs),
            expansion_order=expansion_order,
            include_potential=include_potential,
            include_gmatrix=include_gmatrix,
            include_pseudopotential=include_pseudopotential,
            include_coriolis_coupling=False,  # tell the job not to try to calculate the Coriolis coupling tensor
            calculate_intensities=self.dip is not None # no dipole so we can't get intensities
        )
        if return_analyzer:
            analyzer.print_output_tables(print_intensities=self.dip is not None, print_energies=self.dip is None)
        return analyzer

    def _get_rotated_coordinates(self):
        if self._rot_coords is None:
            self._rot_coords = AnalyticModelBase.transform_coordinates(AnalyticModelBase.transpose(self.rotation))
        return self._rot_coords
    def _get_inverse_coordinates(self):
        if self._inv_coords is None:
            self._inv_coords = AnalyticModelBase.transform_coordinates(self.inverse, coord_vec=self.coords)
        return self._inv_coords

    class SympyExpr:
        def __init__(self, expr, core, ndim):
            self.expr = expr
            self.lam = core
            self.ndim = ndim
        def __call__(self, grid, **kwargs):
            core = self.lam
            ndim = self.ndim
            if grid.ndim == 1:
                return core(grid)
            if grid.shape[-1] == ndim and grid.shape[0] != ndim:
                grid = grid.transpose(np.roll(np.arange(grid.ndim), 1))
            return core(*grid)
        def __repr__(self):
            return "{}({})".format(
                type(self).__name__,
                self.expr
            )
    def wrap_function(self, expr, transform_coordinates=True):
        if isinstance(expr, AnalyticModelBase.numeric_types):
            return float(expr)

        coord_vec = self.coords
        if transform_coordinates and self.rotation is not None:
            coord_vec, new_coords = self._get_rotated_coordinates()
            # subtract off equilibrium values since we have a rotation...
            shift_subs = [(s, s+self.vals[s]) for s in self.coords if s in self.vals]
            expr = expr.subs(shift_subs)
            expr = expr.subs([(old, new) for old,new in zip(self.coords, new_coords)])
        core = sym.lambdify(coord_vec, expr)
        ndim = len(coord_vec)

        return self.SympyExpr(expr, core, ndim)

    def expand_potential(self, order, lambdify=True, evaluate=True):
        potential_expansions = self.v(order, evaluate=evaluate)
        coord_vec = self.coords
        pots = []
        for i, d in enumerate(potential_expansions):
            for j in range(i):
                d = AnalyticModelBase.dot(coord_vec, d, axes=[0, -1])
            pots.append(1/np.math.factorial(i)*d)
        if lambdify:
            pots = self.wrap_function(sum(pots))

        return pots
    def get_DVR_parameters(self,
                           expansion_order=None,
                           lambdify=True,
                           evaluate='constants'
                           ):
        if expansion_order is not None:
            potential_function = self.expand_potential(expansion_order, lambdify=lambdify)
        else:
            potential_function = self.pot
            if lambdify:
                potential_function = self.wrap_function(potential_function)
        gds = self.g(order=2, evaluate=evaluate)
        g = gds[0]
        g_deriv = [gds[2][i][i][i][i] for i in range(len(self.coords))]

        # print(g)
        # print(g_deriv)
        if lambdify:
            g = [
                [
                    self.wrap_function(g[i][j])
                    for j in range(len(self.coords))
                ]
                for i in range(len(self.coords))
            ]
            g_deriv = [
                self.wrap_function(g_deriv[i])
                for i in range(len(self.coords))
            ]
        if len(g_deriv) == 1:
            g = g[0][0]
            g_deriv = g_deriv[0]
        return {
            'potential_function':potential_function,
            'g':g,
            'g_deriv':g_deriv
        }
        # raise Exception(potential_function, g)
    def setup_DVR(self,
                  domain=None, divs=None,
                  use_normal_modes=False,
                  expansion_order=None,
                  **params
                  ):
        if use_normal_modes:
            if self.rotation is None:
                self, _ = self.to_normal_modes(dimensionless=True)
            return self.setup_DVR(
                domain=domain, divs=divs,
                use_normal_modes=False,
                expansion_order=expansion_order,
                **params
            )

        from ..DVR import DVR
        params = dict(self.get_DVR_parameters(expansion_order=expansion_order), **params)
        return DVR(
            domain=domain,
            divs=divs,
            **params
        )

    def evaluate(self, expr, mode='all', numericize=False):
        if self.vals is None:
            raise ValueError("can't evaluate without values for parameters")
        if mode == True:
            mode = 'all'
        a = AnalyticModelBase.eval_exprs(expr, self.constants if mode == 'constants' else self.vals)
        if mode != 'constants':
            try:
                a = np.array(a, dtype=float)
            except ValueError:
                pass
            else:
                if a.ndim == 0:
                    a = a.tolist()
        return a

    def _load_symbols(self):
        sym_list = []
        sym_set = set()
        for x in self.coords:
            for s in x.free_symbols:
                if s not in sym_set:
                    sym_set.add(s)
                    sym_list.append(s)
        self._syms = (tuple(sym_list), tuple(self._parse_symbol(s) for s in sym_list))

    def _parse_symbol(self, sym):
        name = sym.name
        t, a = name.split("[")
        a = a.strip("]").split(",")
        return (t, tuple(int(x) for x in a))

    def jacobian(self, order=0, evaluate=False):
        ics = self.internal_coordinates
        crd = self.coords
        jac = AnalyticModelBase.sym.Matrix([AnalyticModelBase.take_derivs(c, ics) for c in crd]).transpose()
        if order > 0:
            if self.rotation is not None:
                jac = AnalyticModelBase.dot(jac, self.rotation)
            for i in range(order):
                jac = AnalyticModelBase.take_derivs(jac, ics)
        elif evaluate and self.rotation is not None:
                jac = AnalyticModelBase.dot(self.evaluate(jac, mode=evaluate), self.rotation)
        if evaluate:
            jac = self.evaluate(jac, mode=evaluate)
        else:
            jac = AnalyticModelBase.sym.Matrix(jac)
        return jac
    def jacobian_inverse(self, order=0, evaluate=False):
        ics = self.internal_coordinates
        crd = self.coords
        # try:
        #     jac = sym.Array([AnalyticModelBase.take_derivs(c, crd) for c in ics])
        #     for i in range(order):
        #         jac = AnalyticModelBase.take_derivs(jac, crd)
        # except ValueError:
        jac = self.jacobian().inv()
        if order > 0:
            if self.inverse is not None:
                jac = AnalyticModelBase.dot(self.inverse, jac)
            for i in range(order):
                jac = AnalyticModelBase.take_derivs(jac, crd)
        elif evaluate and self.inverse is not None:
                jac = AnalyticModelBase.dot(self.inverse, self.evaluate(jac, mode=evaluate))
        if evaluate:
            jac = self.evaluate(jac, mode=evaluate)
        else:
            jac = AnalyticModelBase.sym.Matrix(jac)
        # if self.inverse is not None:
        #     jac = AnalyticModelBase.dot(jac, self.inverse)
        # for i in range(order):
        #     jac = AnalyticModelBase.take_derivs(jac, crd)
        # if evaluate:
        #     jac = self.evaluate(jac)
        # else:
        #     jac = sym.Matrix(jac)
        return jac
    def _base_gmat(self):
        if self._base_g is None:
            if self._syms is None:
                self._load_symbols()
            symlist = self._syms[1]
            targets = {s.name for s in self._syms[0]}
            if self.vals is not None:
                targets.update(s.name for s in self.vals.keys())
            return [[AnalyticKineticEnergyConstructor.g(a[1], b[1], coord_types=[a[0], b[0]], target_symbols=targets) for b in symlist] for a in symlist]
        return self._base_g
    def g(self, order=0, evaluate=False):
        # Gmatrix elements will basically involve taking some
        # kind of direct product of coordinate reps
        J = self.jacobian_inverse(evaluate=evaluate)
        if self._g is None:
            G = self._base_gmat()
            J_t = J.transpose()
            self._g = AnalyticModelBase.dot(AnalyticModelBase.dot(J_t, G), J)
        all_gs = [self._g]
        if order > 0:
            J_inv = self.jacobian(evaluate=evaluate)
            for i in range(order):
                all_gs.append(AnalyticModelBase.take_derivs(all_gs[-1], self.internal_coordinates))
        if evaluate:
            for i,g in enumerate(all_gs):
                all_gs[i] = self.evaluate(g, mode=evaluate)
        for i in range(1, len(all_gs)):
            g = all_gs[i]
            for j in range(i):
                g = AnalyticModelBase.dot(J_inv, g, axes=[1, -3])
            all_gs[i] = g
        return all_gs
    def v(self, order=2, evaluate=False):
        # we provide a Taylor series expansion of the potential
        v = self.pot
        all_vs = [v]
        if order > 0:
            J_inv = self.jacobian(evaluate=evaluate)
        for i in range(order):
            all_vs.append(AnalyticModelBase.take_derivs(all_vs[-1], self.internal_coordinates))
        if evaluate:
            for i,v in enumerate(all_vs):
                all_vs[i] = self.evaluate(v, mode=evaluate)
        for i in range(1, len(all_vs)):
            v = all_vs[i]
            for j in range(i):
                v = AnalyticModelBase.dot(J_inv, v, axes=[1, -1])
            all_vs[i] = v
        return all_vs

    def _base_u(self):
        if self._syms is None:
            self._load_symbols()
        symlist = self._syms[1]
        targets = {s.name for s in self._syms[0]}
        if self.vals is not None:
            targets.update(s.name for s in self.vals.keys())
        return [[AnalyticKineticEnergyConstructor.vp(a[1], b[1], coord_types=[a[0], b[0]], target_symbols=targets) for b in symlist] for a in symlist]
        # return [[AnalyticKineticEnergyConstructor.vp(a[1], b[1], coord_types=[a[0], b[0]]) for b in symlist] for a in symlist]
    def vp(self, order=0, evaluate=False):
        J_t = self.jacobian()
        if self._u is None:
            U = self._base_u()
            J = AnalyticModelBase.transpose(J_t)
            self._u = 8*AnalyticModelBase.contract(AnalyticModelBase.dot(AnalyticModelBase.dot(J_t, U), J), (0, 1))
        u = self._u
        all_vs = [u]
        if order > 0:
            J_inv = self.jacobian(evaluate=evaluate).transpose()
        for i in range(order):
            all_vs.append(AnalyticModelBase.take_derivs(all_vs[-1], self.internal_coordinates))
        if evaluate:
            for i, v in enumerate(all_vs):
                all_vs[i] = self.evaluate(v, mode=evaluate)
        for i in range(1, len(all_vs)):
            v = all_vs[i]
            for j in range(i):
                v = AnalyticModelBase.dot(J_inv, v, axes=[1, -1])
            all_vs[i] = v
        return all_vs

    def mu(self, order=1, evaluate=False):
        # we provide a Taylor series expansion of the potential
        mu = []
        for v in self.dip:
            all_vs = [v]
            if order > 0:
                J_inv = self.jacobian(evaluate=evaluate)
            for i in range(order):
                all_vs.append(AnalyticModelBase.take_derivs(all_vs[-1], self.internal_coordinates))
            if evaluate:
                for i, v in enumerate(all_vs):
                    all_vs[i] = self.evaluate(v, mode=evaluate)
            for i in range(1, len(all_vs)):
                v = all_vs[i]
                for j in range(i):
                    v = AnalyticModelBase.dot(J_inv, v, axes=[1, -1])
                all_vs[i] = v
            mu.append(all_vs)
        return mu

    Potential = AnalyticPotentialConstructor
    morse = AnalyticPotentialConstructor.morse
    harmonic = AnalyticPotentialConstructor.harmonic
    linear = AnalyticPotentialConstructor.linear
    KE = AnalyticKineticEnergyConstructor

    @classmethod
    def sym(self, base, *args):
        return AnalyticModelBase.symbol(base, *args)
    @classmethod
    def m(self, i):
        return AnalyticModelBase.symbolic_m(i)
    @classmethod
    def r(self, i, j):
        return AnalyticModelBase.symbolic_r(i, j)
    @classmethod
    def a(self, i, j, k):
        return AnalyticModelBase.symbolic_a(i, j, k)
    @classmethod
    def t(self, i, j, k, l):
        return AnalyticModelBase.symbolic_t(i, j, k, l)
    @classmethod
    def y(self, i, j, k, l):
        return AnalyticModelBase.symbolic_t(i, j, k, l)

    convert = UnitsData.convert
    @staticmethod
    def mass(atom):
        return AtomData[atom]["Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")