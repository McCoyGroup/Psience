import abc

import numpy as np

from McUtils.Coordinerds import CoordinateSet
from McUtils.Zachary import CoordinateInterpolator
import McUtils.Numputils as nput
import McUtils.Devutils as dev

from ..Molecools import Molecule
from ..Molecools.Evaluator import EnergyEvaluator

class ProfileGenerator:
    def __init__(self,
                 reactant_complex:Molecule
                 ):
        self.reactants = reactant_complex

    @abc.abstractmethod
    def generate(self, **opts):
        ...

    @classmethod
    def get_profile_generators(cls):
        return {
            'interpolate': InterpolatingProfileGenerator,
            'neb': NudgedElasticBand
        }
    _profile_dispatch = dev.uninitialized
    @classmethod
    def profile_generator_dispatch(cls):
        cls._profile_dispatch = dev.handle_uninitialized(
            cls._profile_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_profile_generators,),
            kwargs=dict(default_method='interpolate')
        )
        return cls._profile_dispatch

    @classmethod
    def resolve(cls, reactants, *args, profile_generator=dev.default, **kwargs) -> 'ProfileGenerator':
        new_cls, opts = cls.profile_generator_dispatch().resolve(profile_generator)
        return new_cls(
            reactants,
            *args,
            **dict(opts, **kwargs)
        )

class OptimizingProfileGenerator(ProfileGenerator):

    def generate(self, initial_search_dir, **opts):
        ...

class InterpolatingProfileGenerator(ProfileGenerator):
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 coordinate_interpolator=None,
                 ts_guesses:list[Molecule]=None,
                 internals=None,
                 num_images=10
                 ):
        super().__init__(reactant_complex)
        self.products = product_complex
        self.intermediates = ts_guesses
        if (
                internals is not None
                and not hasattr(internals, 'converter_options')
                and not callable(internals)
        ):
            spec = self.products.prep_internal_spec(internals)
            base_internals = self.products.modify(internals=spec).internal_coordinates
            internals = self.wrap_conversion(spec, base_internals)
        self.interpolator = self.get_interpolator(coordinate_interpolator, coordinate_system=internals)
        self.num_images = num_images

    def wrap_conversion(self, spec, base_internals):
        def convert(coords):
            new = self.products.modify(coords=coords, internals=spec)
            for k in [
                'origins',
                'axes'
            ]:
                if k in base_internals.converter_options:
                    new.coords.system.converter_options[k] = base_internals.converter_options[k]
            ints = new.internal_coordinates
            for k in [
                'reference_internals',
                'redundant_transformation',
                'origins',
                'axes'
            ]:
                if k in ints.converter_options:
                    ints.system.converter_options[k] = ints.converter_options[k]
            return ints
        return convert

    def get_interpolator(self, coordinate_interpolator, coordinate_system=None):
        if coordinate_interpolator is None:
            return CoordinateInterpolator(
                CoordinateSet(
                    [self.reactants.coords]
                    + [t.coords for t in (self.intermediates if self.intermediates is not None else [])]
                    + [self.products.coords],
                    self.products.coords.system
                ),
                coordinate_system=coordinate_system
            )
        else:
            return coordinate_interpolator

    def generate(self, num_images=None):
        if num_images is None:
            num_images = self.num_images
        return self.interpolator(np.linspace(0, 1, num_images))


class NEBEnergyEvaluator(EnergyEvaluator):
    def __init__(self, base_evaluator:EnergyEvaluator, spring_constant):
        self.base = base_evaluator
        self.energy_units = base_evaluator.energy_units
        self.distance_units = base_evaluator.distance_units
        self.batched_orders = base_evaluator.batched_orders
        self.analytic_derivative_order = base_evaluator.analytic_derivative_order
        self.spring_constant = spring_constant

    def evaluate_term(self, coords, order, **opts):
        return self.base.evaluate_term(coords, order, **opts)

    def image_pairwise_contribution(self, guess, cur, prev, next, order=0):
        if order > 2: return 0

        (prev_dist, prev_const), (next_dist, next_const) = self.get_spring_params(guess, cur, prev, next)
        contribution = 0
        if order == 0:
            if prev_dist is not None:
                contribution += (prev_dist ** 2) * prev_const
            if next_dist is not None:
                contribution += (next_dist ** 2) * next_const
        elif order == 1:
            if prev_dist is not None:
                contribution += 2 * prev_const * (guess[:, cur] - guess[:, prev])
            if next_dist is not None:
                contribution += 2 * next_const * (guess[:, cur] - guess[:, next])
        elif order == 2:
            if prev_dist is not None:
                contribution += 2 * prev_const * nput.identity_tensors(guess.shape[0], guess.shape[1])
            if next_dist is not None:
                contribution += 2 * next_const * nput.identity_tensors(guess.shape[0], guess.shape[1])
        return contribution

class NudgedElasticBand(InterpolatingProfileGenerator):
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 energy_evaluator: EnergyEvaluator,
                 coordinate_interpolator=None,
                 num_images=10,
                 spring_constant=.1
                 ):
        super().__init__(
            reactant_complex,
            product_complex,
            coordinate_interpolator=coordinate_interpolator,
            num_images=num_images
        )
        self.energy_evaluator = energy_evaluator
        self.spring_constant = spring_constant

    def _get_potential_expansion(self):
        ...

    def generate(self, num_images=None, spring_constant=None, energy_evaluator=None):
        base_images = super().generate(num_images=num_images)
        if spring_constant is None:
            spring_constant = self.spring_constant
        if energy_evaluator is None:
            energy_evaluator = self.energy_evaluator

        image_coords = [b.coords for b in base_images]
        return nput.iterative_chain_minimize(
            image_coords,


        )

