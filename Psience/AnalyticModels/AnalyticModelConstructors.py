import itertools

import numpy as np, scipy
import sympy as sym

from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
from .Helpers import AnalyticModelBase
from ..Data import KEData

__all__ = [
    'AnalyticPotentialConstructor',
    'AnalyticKineticEnergyConstructor',
    'AnalyticModel',
    "MolecularModel"
]
__reload_hook__ = ["..Data"]

class AnalyticPotentialConstructor(AnalyticModelBase):
    """
    Provides a set of symbolic potentials for use in models

    :related:AnalyticModel, AnalyticKineticEnergyConstructor
    """
    @classmethod
    def morse(cls, i, j, De=None, a=None, re=None, eq=None, w=None, wx=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """

        if De is None and w is not None:
            if isinstance(w, str):
                w = AnalyticModelBase.symbol(w, i, j)
            if isinstance(wx, str):
                wx = AnalyticModelBase.symbol(wx, i, j)
            muv = (
                    1/AnalyticModelBase.symbolic_m(i)
                    + 1/AnalyticModelBase.symbolic_m(j)
            )
            De = (w ** 2) / (4 * wx)
            a = sym.sqrt(2 * wx / muv)

        if re is None:
            re = eq

        return cls.calc_morse(
            AnalyticModelBase.symbol("De", i, j) if De is None else De,
            AnalyticModelBase.symbol("ap", i, j) if a is None else a,
            AnalyticModelBase.symbolic_r(i, j),
            AnalyticModelBase.symbol("re", i, j) if re is None else re
        )
    @staticmethod
    def calc_morse(De, a, r, re):
        return De * (1 - AnalyticModelBase.sym.exp(-a * (r-re)))**2

    @staticmethod
    def harm(k, x, x_e):
        return k*(x-x_e)**2
    @classmethod
    def harmonic(cls, *args, k=None, eq=None, qe=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """

        if qe is None:
            qe = eq

        return cls.harm(
            AnalyticModelBase.symbol("k", *args) if k is None else k,
            AnalyticModelBase.var(*args),
            AnalyticModelBase.symbol("qe", *args) if qe is None else qe,
        )


    @staticmethod
    def lin(k, x, x_e):
        return k * (x - x_e)
    @classmethod
    def linear(cls, *args, k=1, eq=None, xe=None):
        """
        Returns a fully symbolic form of a linear function
        :return:
        :rtype:
        """
        if xe is None:
            xe = eq
        return cls.lin(
            AnalyticModelBase.symbol("k", *args) if k is None else k,
            AnalyticModelBase.var(*args),
            AnalyticModelBase.symbol("xe", *args) if xe is None else xe
        )

    @staticmethod
    def pow(k, x, x_e, n):
        return k * (x - x_e) ** n
    @classmethod
    def power(cls, *args, k=1, eq=None, n=None, xe=None):
        """
        Returns a fully symbolic form of a linear function
        :return:
        :rtype:
        """
        if xe is None:
            xe = eq
        return cls.pow(
            AnalyticModelBase.symbol("k", *args) if k is None else k,
            AnalyticModelBase.var(*args),
            AnalyticModelBase.symbol("xe", *args) if xe is None else xe,
            AnalyticModelBase.symbol("n", *args) if n is None else n
        )

    @classmethod
    def cos(cls, *args, eq=None, qe=None):
        """
        Returns a fully symbolic form of cos
        :return:
        :rtype:
        """
        if qe is None:
            qe = eq
        return AnalyticModelBase.sym.cos(
            AnalyticModelBase.var(*args) -
                (AnalyticModelBase.symbol("qe", *args) if qe is None else qe)
        )

    @classmethod
    def sin(cls, *args, eq=None, qe=None):
        """
        Returns a fully symbolic form of sin
        :return:
        :rtype:
        """
        if qe is None:
            qe = eq
        return AnalyticModelBase.sym.sin(
            AnalyticModelBase.var(*args) -
            (AnalyticModelBase.symbol("qe", *args) if qe is None else qe)
        )

    @classmethod
    def multiwell(cls, *args, turning_points=None, origin=None, eq=None, minimum=0, depth=None):
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
        if origin is None:
            origin = eq
        if origin is None:
            origin = 0
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
        swapped = sorting.index(coord_types[0]) > sorting.index(coord_types[1])
        if swapped:
            inds1, inds2 = inds2, inds1
            type1, type2 = type2, type1

        remapping = {}
        n = 1
        _ = []
        for i in inds1:
            if i not in remapping:
                remapping[i] = n
                _.append(n)
                n += 1
            else:
                _.append(remapping[i])
        i1 = tuple(_)
        _ = []
        for i in inds2:
            if i not in remapping:
                remapping[i] = n
                _.append(n)
                n += 1
            else:
                _.append(remapping[i])
        i2 = tuple(_)

        return ((type1, type2), i1, i2), remapping, swapped

    @classmethod
    def kinetic_exprs(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None):
        key, mapping, swapped = cls._get_coord_key(inds1, inds2, coord_types=coord_types)
        (expr_g, expr_vp), perm = KEData.find_expressions(key, return_permutation=True)
        if perm is not None:
            shared_indices = np.intersect1d(inds1, inds2)
            p1, p2 = perm
            if swapped:
                p1, p2 = p2, p1
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

        if not (isinstance(expr_g, (int, np.integer)) and isinstance(expr_vp, (int, np.integer))):
            rev_mapping = {v: k for k, v in mapping.items()}
            if not isinstance(expr_g, (int, np.integer)):
                subs = tuple(
                    (s, cls.reindex_symbol(s, rev_mapping, target_symbols=target_symbols)) for s in expr_g.free_symbols)
                expr_g = expr_g.subs([(s, q) for s, (q, _) in subs])
                expr_g = expr_g.subs([(q, f) for s, (q, f) in subs])
            if not isinstance(expr_vp, (int, np.integer)):
                subs = tuple(
                    (s, cls.reindex_symbol(s, rev_mapping, target_symbols=target_symbols)) for s in expr_vp.free_symbols)
                expr_vp = expr_vp.subs([(s, q) for s, (q, _) in subs])
                expr_vp = expr_vp.subs([(q, f) for s, (q, f) in subs])

        return expr_g, expr_vp

    @classmethod
    def g(cls, inds1:'Iterable[int]', inds2:'Iterable[int]', coord_types=None, target_symbols=None):
        return cls.kinetic_exprs(inds1, inds2, coord_types=coord_types, target_symbols=target_symbols)[0]
    @classmethod
    def vp(cls, inds1:'Iterable[int]', inds2:'Iterable[int]', coord_types=None, target_symbols=None):
        return cls.kinetic_exprs(inds1, inds2, coord_types=coord_types, target_symbols=target_symbols)[1]

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

    @classmethod
    def from_potential(cls, potential, dipole=None, values=None, rotation=None):
        # we infer the coordinates from the terms that appear in the potential

        coords = [
            x for x in
            potential.free_symbols
            if x.name.split("[", 1)[0] in {"r", "a", "t", "y"}
        ]
        return cls(
            coords,
            potential,
            dipole=dipole,
            values=values,
            rotation=rotation
        )


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
                atoms=None, coords=None,
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
        if atoms is None:
            atoms = ["H"] * natoms
        if coords is None:
            coords = np.zeros((natoms, 3))
        analyzer = runner(
            [atoms, coords],  # dummy data since this won't be used at all
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
        def _broadcast_tree(self, shape, expr_list):
            if not isinstance(expr_list, (np.ndarray,) + AnalyticModelBase.numeric_types):
                return [self._broadcast_tree(shape, l) for l in expr_list]
            elif isinstance(expr_list, np.ndarray):
                return expr_list
            else:
                return np.full(shape, expr_list)
        def __call__(self, grid, vector=None, **kwargs):
            core = self.lam
            ndim = self.ndim
            grid = np.asanyarray(grid)
            if grid.ndim == 1:
                vals = core(*grid)
                if not isinstance(vals, AnalyticModelBase.numeric_types):
                    vals = np.array(vals)
                return vals

            if vector is None:
                vector = grid.shape[-1] == ndim and grid.shape[0] != ndim
            if vector:
                grid = np.moveaxis(grid, -1, 0) #grid.transpose(np.roll(np.arange(grid.ndim), 1))

            vals = core(*grid)
            # broadcast appropriately
            if not isinstance(vals, np.ndarray):
                if isinstance(vals, AnalyticModelBase.numeric_types):
                    vals = np.full(grid.shape[1:], vals)
                else:
                    vals = np.array(self._broadcast_tree(grid.shape[1:], vals))
                    for _ in range(grid.ndim - 1):
                        vals = np.moveaxis(vals, -1, 0)
            return vals
        def __repr__(self):
            return "{}({})".format(
                type(self).__name__,
                self.expr
            )
    def wrap_function(self, expr, transform_coordinates=True, mode=None):
        coord_vec = self.coords
        ndim = len(coord_vec)
        if isinstance(expr, AnalyticModelBase.numeric_types):
            if isinstance(expr, np.ndarray):
                val = expr
            else:
                val = float(expr)
            return self.SympyExpr(val, lambda c,*r:np.full(c.shape, val), ndim)

        if transform_coordinates and self.rotation is not None:
            coord_vec, new_coords = self._get_rotated_coordinates()
            # subtract off equilibrium values since we have a rotation...
            shift_subs = [(s, s+self.vals[s]) for s in self.coords if s in self.vals]
            expr = expr.subs(shift_subs)
            expr = expr.subs([(old, new) for old,new in zip(self.coords, new_coords)])
        core = sym.lambdify(coord_vec, expr)


        return self.SympyExpr(expr, core, ndim)

    def expand_potential(self, order, lambdify=True, evaluate=True, contract=True):
        potential_expansions = self.v(order, evaluate=evaluate)
        coord_vec = self.coords
        if contract:
            pots = []
            for i, d in enumerate(potential_expansions):
                for j in range(i):
                    d = AnalyticModelBase.dot(coord_vec, d, axes=[0, -1])
                pots.append(1/np.math.factorial(i)*d)
        else:
            pots = potential_expansions
        if lambdify:
            pots = self.wrap_function(sum(pots)) if contract else [self.wrap_function(p) for p in pots]

        return pots
    def get_DVR_parameters(self,
                           expansion_order=None,
                           lambdify=True,
                           evaluate='constants'
                           ):
        if expansion_order is not None:
            potential_function = self.expand_potential(expansion_order, lambdify=lambdify)
        else:
            if evaluate:
                potential_function = self.evaluate(self.pot, mode=evaluate)
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
        # if len(g_deriv) == 1:
        #     g = g[0][0]
        #     g_deriv = g_deriv[0]
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

    def jacobian(self, order=0, evaluate=False, lambdify=False):
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
        if lambdify:
            jac = self.wrap_function(jac)
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
    def g(self, order=0, evaluate=False, lambdify=False):
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
        if lambdify:
            all_gs = [self.wrap_function(g) for g in all_gs]
        return all_gs
    def v(self, order=2, evaluate=False, lambdify=False):
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
        if lambdify:
            all_vs = [self.wrap_function(v) for v in all_vs]
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
    def vp(self, order=0, evaluate=False, lambdify=False):
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
        if lambdify:
            all_vs = [self.wrap_function(u) for u in all_vs]
        return all_vs

    def mu(self, order=1, evaluate=False, lambdify=False):
        # we provide a Taylor series expansion of the dipole
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
            if lambdify:
                all_vs = [self.wrap_function(u) for u in all_vs]
            mu.append(all_vs)
        return mu

    Potential = AnalyticPotentialConstructor
    morse = AnalyticPotentialConstructor.morse
    harmonic = AnalyticPotentialConstructor.harmonic
    linear = AnalyticPotentialConstructor.linear
    power = AnalyticPotentialConstructor.power
    sin = AnalyticPotentialConstructor.sin
    cos = AnalyticPotentialConstructor.cos
    KE = AnalyticKineticEnergyConstructor

    class NamespaceContext:
        def __init__(self, context=None):
            self._context = context
            self._additions = set()
        def _get_frame_vars(self):
            import inspect
            frame = inspect.currentframe()
            parent = frame.f_back.f_back.f_back
            # print(parent.f_locals)
            return parent.f_locals
        def insert_vars(self):
            globs = self._context
            if globs is None:
                globs = self._get_frame_vars()
            for k, v in dict(
                    r=AnalyticModel.r,
                    a=AnalyticModel.a,
                    cos=AnalyticModel.cos,
                    sin=AnalyticModel.sin,
                    morse=AnalyticModel.morse,
                    harmonic=AnalyticModel.harmonic,
                    sym=AnalyticModel.sym,
                    m=AnalyticModel.m,
            ).items():
                globs[k] = v
                self._additions.add(k)
        def prune_vars(self):
            globs = self._context
            if globs is None:
                globs = self._get_frame_vars()
            for x in self._additions:
                try:
                    del globs[x]
                except KeyError:
                    pass
        def __enter__(self):
            self.insert_vars()
        def __exit__(self, exc_type, exc_val, exc_tb):
            self.prune_vars()

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

    def molecular_potential(self, mol):
        if mol.internals is None:
            raise ValueError("Molecule {} needs internal coordinates to be used with a model potential")
        return MolecularModelPotentialFunction(self, mol)
    def molecular_dipole(self, mol):
        if mol.internals is None:
            raise ValueError("Molecule {} needs internal coordinates to be used with a model dipole")
        return MolecularModelDipoleFunction(self, mol)
    def molecular_gmatrix(self, mol):
        if mol.internals is None:
            raise ValueError("Molecule {} needs internal coordinates to be used with a model gmatrix")
        return MolecularModelGMatrixFunction(self, mol)

class MolecularModel(AnalyticModel):
    def __init__(self, mol, coords, potential, dipole=None, values=None, rotation=None):
        super().__init__(coords, potential, dipole=dipole, values=values, rotation=rotation)
        self.mol = mol
    @property
    def potential(self):
        return MolecularModelPotentialFunction(self, self.mol)
    @property
    def gmatrix(self):
        return MolecularModelGMatrixFunction(self, self.mol)
    @property
    def vprime(self):
        return MolecularModelVPrimeFunction(self, self.mol)
    @property
    def dipole(self):
        return MolecularModelDipoleFunction(self, self.mol)

class MolecularModelFunction:

    def __init__(self, deriv_evaluator, mol):
        self.mol = mol
        self.deriv_evaluator = deriv_evaluator
        self.derivs = []

    def evaluate(self, carts, deriv_order=None, internals=False, which=None, sel=None, axes=None, derivs=None):
        carts = np.asanyarray(carts)

        ads = 1 if deriv_order is None else deriv_order + 1
        if derivs is None:
            if len(self.derivs) < ads + 1:
                self.derivs = self.deriv_evaluator(order=ads-1, evaluate='constants', lambdify=True)
            derivs = self.derivs
        elif len(derivs) < ads:
            raise ValueError("need at least {} derivs".format(ads - 1))


        base_shape = carts.shape[:-2]
        carts = carts.reshape((-1,) + carts.shape[-2:])
        if (
                internals
                or which is not None
                or sel is not None
                or axes is not None
        ):
            carts = self.mol.get_displaced_coordinates(
                carts,
                which=which, sel=sel, axes=axes, internals=internals,
                shift=False
            )

        vals = self.mol.evaluate_at(
            lambda c, deriv_order=None: (
                derivs[0](c, vector=True)
                    if deriv_order is None else
                [
                    df(c, vector=True) for df in derivs[:deriv_order + 1]
                ]
            ),
            carts,
            deriv_order=None if deriv_order is not None and deriv_order == 0 else deriv_order,
            internals=True,
            strip_embedding=True
        )
        if deriv_order is not None:
            if deriv_order == 0:
                vals = [vals]

            take = None
            if which is not None:
                take = tuple(
                    np.ravel_multi_index(idx, (3, len(self.mol.masses)))
                    if not isinstance(idx, (int, np.integer)) else
                        idx
                    for idx in which
                )
            elif sel is not None or axes is not None:
                if sel is None:
                    sel = np.arange(len(self.mol.masses))
                if axes is None:
                    axes = np.arange(3)
                take = np.ravel_multi_index(np.array(list(itertools.product(sel, axes))).T, (3, len(self.mol.masses)))

            if take is not None:
                new = []
                for n, d in enumerate(vals):
                    for j in range(n):
                        d = np.take(d, take, axis=j + 1)
                    d = d.reshape(base_shape + d.shape[1:])
                    new.append(d)
                vals = new

            vals = [
                v.reshape(base_shape + v.shape[1:])
                for v in vals
            ]

        else:
            vals = vals.reshape(base_shape)

        return vals

    def __call__(self, carts, deriv_order=None, internals=False, which=None, sel=None, axes=None):
        return self.evaluate(carts, deriv_order=deriv_order, internals=internals, which=which, sel=sel, axes=axes)

class MolecularModelPotentialFunction(MolecularModelFunction):
    def __init__(self, model, mol):
        super().__init__(model.v, mol)

class MolecularModelDipoleFunction(MolecularModelFunction):
    def __init__(self, model, mol):
        super().__init__(model.mu, mol)
        self._subfuncs = []

    def evaluate(self, carts, deriv_order=None, internals=False, which=None, sel=None, axes=None):
        """
        This has the added complication of needing to dispatch over the axes...

        :param carts:
        :type carts:
        :param deriv_order:
        :type deriv_order:
        :param internals:
        :type internals:
        :param which:
        :type which:
        :param sel:
        :type sel:
        :param axes:
        :type axes:
        :return:
        :rtype:
        """
        ads = 1 if deriv_order is None else deriv_order + 1
        if len(self.derivs) < ads + 1:
            self.derivs = self.deriv_evaluator(order=ads - 1, evaluate='constants', lambdify=True)

        vals = []
        for d in self.derivs:
            vals.append(
                super().evaluate(carts, deriv_order=deriv_order, internals=internals, which=which, sel=sel, axes=axes, derivs=d)
            )

        if deriv_order is None:
            vals = np.moveaxis(np.array(vals), 0, -1)
        else:
            transpose_terms = [[] for _ in range(deriv_order + 1)]
            for term_list in vals:
                for i, t in enumerate(term_list):
                    transpose_terms[i].append(t)
            vals = [
                np.moveaxis(np.array(tl), 0, -1)
                for tl in transpose_terms
            ]

        return vals

class MolecularModelGMatrixFunction(MolecularModelFunction):
    def __init__(self, model, mol):
        super().__init__(model.g, mol)

class MolecularModelVPrimeFunction(MolecularModelFunction):
    def __init__(self, model, mol):
        super().__init__(model.vp, mol)