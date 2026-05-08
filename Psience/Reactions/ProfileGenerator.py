import abc

import numpy as np

from McUtils.Coordinerds import CoordinateSet
from McUtils.Zachary import CoordinateInterpolator
import McUtils.Numputils as nput
import McUtils.Devutils as dev
from McUtils.Data import UnitsData

from ..Molecools import Molecule
from ..Molecools.Evaluator import EnergyEvaluator

__all__ = [
    "ProfileGenerator"
]

class ProfileGenerator:
    def __init__(self,
                 reactant_complex:Molecule
                 ):
        self.reactants = reactant_complex

    @abc.abstractmethod
    def generate(self, **opts):
        ...

    def evaluate_profile_energies(self, profile:'list[Molecule]', **opts):
        raise NotImplementedError(f"{type(self).__name__} doesn't support profile energy evaluation")
    def evaluate_profile_distances(self, profile:'list[Molecule]', **opts):
        raise NotImplementedError(f"{type(self).__name__} doesn't support profile energy evaluation")

    profile_registry = {}
    @classmethod
    def register(cls, name, profile=None):
        if profile is not None:
            cls.profile_registry[name] = profile
            return profile
        else:
            def register(profile):
                return cls.register(name, profile)
            return register
    @classmethod
    def get_profile_generators(cls):
        return cls.profile_registry | {
            'interpolate': InterpolatingProfileGenerator,
            'neb': NudgedElasticBand,
            'ase-neb': ASENEBGenerator,
            'ase-dimer': ASEDimerGenerator,
            'pys-neb': PysisNEBGenerator,
            'pys-gsm': PysisGSMGenerator,
            'pys-fsm': PysisFSMGenerator,
            'pys-cos': PysisCOSGenerator,
        }
    _profile_dispatch = dev.uninitialized
    @classmethod
    def profile_generator_dispatch(cls):
        cls._profile_dispatch = dev.handle_uninitialized(
            cls._profile_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_profile_generators,),
            # kwargs=dict(default_method='interpolate')
        )
        return cls._profile_dispatch

    @classmethod
    def resolve(cls, reactants, *args, profile_generator=dev.default, **kwargs) -> 'ProfileGenerator':
        new_cls, opts = cls.profile_generator_dispatch().resolve(profile_generator)
        if new_cls is None:
            raise ValueError(f"couldn't resolve profile generator {profile_generator}")
        return new_cls(
            reactants,
            *args,
            **dict(opts, **kwargs)
        )



# class PyGSMProfileGenerator(ProfileGenerator):
#
#     def __init__(self,
#                  reactant_complex:Molecule,
#                  constraints=None
#                  ):
#         super().__init__(reactant_complex)

# class OptimizingProfileGenerator(ProfileGenerator):
#
#     def generate(self, initial_search_dir, **opts):
#         ...

class InterpolatingProfileGenerator(ProfileGenerator):
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 coordinate_interpolator=None,
                 ts_guesses:list[Molecule]=None,
                 internals=None,
                 num_images=10,
                 initial_image_positions=None,
                 max_displacement_step=None
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
        self.interpolator = self.get_interpolator(coordinate_interpolator,
                                                  coordinate_system=internals,
                                                  max_displacement_step=max_displacement_step
                                                  )
        self.num_images = num_images
        self.initial_image_positions = initial_image_positions

    def prep_cart_coordinate_system(self, cart_sys, base_internals):
        for k in [
            'origins',
            'axes',
            'masses'
        ]:
            if k in base_internals.converter_options:
                cart_sys.converter_options[k] = base_internals.converter_options[k]
        return cart_sys
    def prep_int_coordinate_system(self, ints, new_coords, base_internals):
        for k in [
            'reference_internals',
            'redundant_transformation',
            'reference_coordinates',
            'origins',
            'axes',
            'masses'
        ]:
            if k in ints.converter_options and k not in ints.system.converter_options:
                ints.system.converter_options[k] = ints.converter_options[k]
        return ints
    def wrap_conversion(self, spec, base_internals):
        def convert(coords):
            new = self.products.modify(coords=coords, internals=spec)
            self.prep_cart_coordinate_system(new.coords.system, base_internals)
            ints = new.internal_coordinates
            self.prep_int_coordinate_system(ints, new.coords, base_internals)
            return ints
        convert.system = base_internals.system
        return convert

    def get_interpolator(self,
                         coordinate_interpolator,
                         coordinate_system=None,
                         max_displacement_step=None
                         ):
        if coordinate_interpolator is None:
            if max_displacement_step is None:
                if hasattr(coordinate_system, 'system'):
                    sys = coordinate_system.system
                else:
                    sys = coordinate_system
                if sys is None:
                    max_displacement_step = 1.0
                elif hasattr(sys, 'name'):
                    if 'ZMatrix' in sys.name:
                        max_displacement_step = .5
                    elif 'Cartesian' in sys.name:
                        max_displacement_step = 1.0
                    else:
                        max_displacement_step = .5
                else:
                    max_displacement_step = 1.0
            return CoordinateInterpolator(
                CoordinateSet(
                    [self.reactants.coords]
                    + [t.coords for t in (self.intermediates if self.intermediates is not None else [])]
                    + [self.products.coords],
                    self.products.coords.system
                ),
                coordinate_system=coordinate_system,
                reembed=True,
                embedding_options={'masses':self.reactants.atomic_masses},
                max_displacement_step=max_displacement_step
            )
        else:
            return coordinate_interpolator

    def generate(self, num_images=None, initial_image_positions=None):
        if initial_image_positions is None:
            initial_image_positions = self.initial_image_positions
        if initial_image_positions is None:
            if num_images is None:
                num_images = self.num_images
            initial_image_positions = np.linspace(0, 1, num_images)
        crds = self.interpolator(initial_image_positions)
        images = [
            (self.reactants if p < .8 else self.products).modify(coords=c)
            for p,c in zip(initial_image_positions, crds)
        ]
        return images

class NudgedElasticBand(InterpolatingProfileGenerator):
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 *,
                 energy_evaluator: EnergyEvaluator,
                 coordinate_interpolator=None,
                 num_images=10,
                 initial_image_positions=None,
                 spring_constant=.01,
                 internals=None,
                 max_displacement_step=None,
                 interpolation_gradient_scaling=None
                 ):
        self.interpolation_gradient_scaling = interpolation_gradient_scaling
        self._energy_evaluator = energy_evaluator
        super().__init__(
            reactant_complex,
            product_complex,
            coordinate_interpolator=coordinate_interpolator,
            num_images=num_images,
            initial_image_positions=initial_image_positions,
            internals=internals,
            max_displacement_step=max_displacement_step
        )
        self.spring_constant = spring_constant
        self.internals = internals
    @property
    def energy_evaluator(self):
        return self.reactants.get_energy_evaluator(self._energy_evaluator)
    def _potential(self, energy_evaluator):
        def _func(coords, mask):
            coords = coords.reshape(coords.shape[:-1] + (-1, 3))
            coords = coords * UnitsData.convert("BohrRadius", energy_evaluator.distance_units)
            return energy_evaluator.evaluate(coords, order=0)[0]
        return _func
    def _jacobian(self, energy_evaluator, remove_translation_rotations=True):
        def _jac(coords, mask):
            coords = coords.reshape(coords.shape[:-1] + (-1, 3))
            coords = coords * UnitsData.convert("BohrRadius", energy_evaluator.distance_units)
            grad = energy_evaluator.evaluate(coords, order=[1])
            if remove_translation_rotations:
                grad = nput.remove_translation_rotations(grad,
                                                         coords,
                                                         masses=self.reactants.masses,
                                                         mass_weighted=False)
            return grad[0]
        return _jac

    # def prep_cart_coordinate_system(self, cart_sys, base_internals):
    #     super().prep_cart_coordinate_system(cart_sys, base_internals)
    #     if self.interpolation_gradient_scaling is not None:
    #         cart_sys.converter_options['gradient_function'] = self._jacobian(self.energy_evaluator)
    #         cart_sys.converter_options['gradient_scaling'] = self.interpolation_gradient_scaling
    #     return cart_sys
    def prep_int_coordinate_system(self, ints, new_coords, base_internals):
        super().prep_int_coordinate_system(ints, new_coords, base_internals)
        if self.interpolation_gradient_scaling is not None:
            ints.system.converter_options['gradient_function'] = self._jacobian(self.energy_evaluator)
            ints.system.converter_options['gradient_scaling'] = self.interpolation_gradient_scaling
        return ints

    def evaluate_profile_distances(self, profile:'list[Molecule]', normalize=True):
        step_finder = self.get_step_finder()
        dists = [
            step_finder.get_dist(
                p1.coords.flatten()[np.newaxis],
                p2.coords.flatten()[np.newaxis]
            )[0]
            for p1,p2 in zip(profile, profile[1:])
        ]
        d = np.cumsum([0] + dists)
        if normalize:
            d = d / d[-1]
        return d

    def evaluate_profile_energies(self, profile:'list[Molecule]', energy_evaluator=None):
        if energy_evaluator is None:
            energy_evaluator = self.energy_evaluator

        base_pot = self._potential(energy_evaluator)
        return [
            base_pot(i.coords.reshape(1, -1), None)[0]
            for i in profile
        ]

    def get_step_finder(self, spring_constant=None, energy_evaluator=None):
        if spring_constant is None:
            spring_constant = self.spring_constant
        if energy_evaluator is None:
            energy_evaluator = self.energy_evaluator
        return nput.NudgedElasticBandStepFinder(
                self._potential(energy_evaluator),
                self._jacobian(energy_evaluator),
                spring_constants=spring_constant
            )

    def generate(self, num_images=None, spring_constant=None, energy_evaluator=None, return_preopt=False,
                 embedding_options=None, base_images=None,
                 **opt_opts):
        if base_images is None:
            base_images = super().generate(num_images=num_images)
        step_finder = self.get_step_finder(spring_constant=spring_constant, energy_evaluator=energy_evaluator)

        if embedding_options is None:
            embedding_options = {'masses':self.reactants.masses}
        image_coords = [b.coords.flatten() for b in base_images]
        (res, nimg), _, _ = nput.iterative_chain_minimize(
            image_coords,
            step_finder,
            embedding_options=embedding_options,
            **opt_opts
        )

        new_img = [i.modify(coords=c.reshape(-1, 3)) for i,c in zip(base_images, res)]
        if return_preopt:
            return base_images, new_img
        else:
            return new_img

class GrowingString(ProfileGenerator):
    ...

class ASEProfileGenerator(InterpolatingProfileGenerator):
    default_method: str
    default_optimizer_method: str
    default_optimizer: str
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 *,
                 energy_evaluator: EnergyEvaluator,
                 coordinate_interpolator='ase',
                 num_images=10,
                 initial_image_positions=None,
                 internals=None,
                 max_displacement_step=None,
                 interpolation_gradient_scaling=None,
                 intermediates=None,
                 method=None,
                 optimizer_method=None,
                 optimizer=None,
                 **opts
                 ):
        self.interpolation_gradient_scaling = interpolation_gradient_scaling
        self._energy_evaluator = energy_evaluator
        if dev.str_is(coordinate_interpolator, 'ase'):
            if intermediates is not None:
                pre_traj = [reactant_complex] + list(intermediates) + [product_complex]
            else:
                pre_traj = [reactant_complex, product_complex]
            coordinate_interpolator = self.ASECoordinateInterpolator(pre_traj)
        super().__init__(
            reactant_complex,
            product_complex,
            coordinate_interpolator=coordinate_interpolator,
            num_images=num_images,
            initial_image_positions=initial_image_positions,
            internals=internals
        )
        self.opts = opts
        self.max_step = max_displacement_step
        self.num_images = num_images
        self.initial_image_positions = initial_image_positions
        if method is None:
            method = self.default_method
        self.method = method
        if optimizer_method is None:
            optimizer_method = self.default_optimizer_method
        self.optimizer_method = optimizer_method
        if optimizer is None:
            optimizer = self.default_optimizer
        self.optimizer = optimizer

    class ASECoordinateInterpolator:
        def __init__(self, initial_path):
            self.path = initial_path

    def prep_images(self,
                    num_images=None,
                    energy_evaluator=None,
                    base_images=None
                    ):

        if base_images is None:
            base_images = super().generate(num_images=num_images)

        if energy_evaluator is None:
            energy_evaluator = self._energy_evaluator
        if energy_evaluator is not None:
            base_images = [
                b.modify(energy_evaluator=energy_evaluator)
                for b in base_images
            ]

        base_images: list[Molecule]
        images = [
            img.to_ase()
            for img in base_images
        ]

        return images

    def evaluate_profile_distances(self, profile:'list[Molecule]', normalize=True):
        dists = [
            p1.get_rmsd(p2)
            for p1, p2 in zip(profile[:-1], profile[1:])
        ]
        d = np.cumsum([0] + dists)
        if normalize:
            d = d / d[-1]
        return d

    def evaluate_profile_energies(self, profile:'list[Molecule]', energy_evaluator=None):
        if energy_evaluator is None:
            energy_evaluator = self._energy_evaluator

        return [
            p.modify(energy_evaluator=energy_evaluator).calculate_energy()
            for p in profile
        ]

    def generate(self,
                 num_images=None,
                 energy_evaluator=None,
                 base_images=None,
                 *,
                 method=None,
                 method_options=None,
                 optimizer_method=None,
                 optimizer=None,
                 max_step=None,
                 **opt_opts):
        images = self.prep_images(
            num_images=num_images,
            energy_evaluator=energy_evaluator,
            base_images=base_images
        )


        if method is None:
            method = self.method
        if optimizer_method is None:
            optimizer_method = self.optimizer_method
        if optimizer is None:
            optimizer = self.optimizer

        if isinstance(method, str):
            method = {
                "method":method,
            }
        if method_options is not None:
            method = method | method_options
        method = self.opts | method

        if max_step is None:
            max_step = self.max_step

        if max_step is not None:
            opt_opts["max_step"] = max_step

        _, images, _ = images[0].optimize_trajectory(images,
                                                     method,
                                                     optimizer=optimizer,
                                                     optimizer_method=optimizer_method,
                                                     return_coords=False,
                                                     **opt_opts
                                                     )
        images = [Molecule.from_ase(i) for i in images]

        return images

class ASENEBGenerator(ASEProfileGenerator):
    default_method = 'neb'
    default_optimizer_method = 'improvedtangent'
    default_optimizer = 'FIRE'
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 *,
                 energy_evaluator: EnergyEvaluator,
                 spring_constant=.1,
                 **opts
    ):
        super().__init__(
            reactant_complex,
            product_complex,
            energy_evaluator=energy_evaluator,
            **opts
        )
        self.spring_constant = spring_constant


    def generate(self,
                 num_images=None,
                 k=None,
                 spring_constant=None,
                 energy_evaluator=None,
                 return_preopt=False,
                 base_images=None,
                 method=None,
                 optimizer_method=None,
                 optimizer=None,
                 **opt_opts):
        if spring_constant is None:
            spring_constant = self.spring_constant
        if k is None:
            k = spring_constant
        return super().generate(
            num_images=num_images,
            energy_evaluator=energy_evaluator,
            base_images=base_images,
            method=method,
            method_options={'k':k},
            optimizer_method=optimizer_method,
            optimizer=optimizer,
            **opt_opts
        )

class ASEDimerGenerator(ASEProfileGenerator):
    default_method = 'dimer'
    default_optimizer_method = None
    default_optimizer = 'minmode'
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 *,
                 climb=True, # ignored
                 **opts
                 ):
        if not climb:
            raise ValueError(f"{type(self).__name__} only supports `climb=True`")
        super().__init__(
            reactant_complex,
            product_complex,
            **opts
        )

class PysisyphusProfileGenerator(InterpolatingProfileGenerator):
    default_coord_type = "cartesian"
    default_method: str
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 *,
                 energy_evaluator: EnergyEvaluator,
                 coordinate_interpolator='pysis',
                 num_images=10,
                 initial_image_positions=None,
                 internals=None,
                 max_displacement_step=None,
                 interpolation_gradient_scaling=None,
                 intermediates=None,
                 coord_type=None,
                 **opts
                 ):
        self.interpolation_gradient_scaling = interpolation_gradient_scaling
        self._energy_evaluator = energy_evaluator
        if dev.str_is(coordinate_interpolator, 'pysis'):
            if intermediates is not None:
                pre_traj = [reactant_complex] + list(intermediates) + [product_complex]
            else:
                pre_traj = [reactant_complex, product_complex]
            coordinate_interpolator = self.PysisCoordinateInterpolator(pre_traj)
        super().__init__(
            reactant_complex,
            product_complex,
            coordinate_interpolator=coordinate_interpolator,
            num_images=num_images,
            initial_image_positions=initial_image_positions,
            internals=internals,
            max_displacement_step=max_displacement_step
        )
        self.num_images = num_images
        self.initial_image_positions = initial_image_positions
        if coord_type is None:
            coord_type = self.default_coord_type
        self.coord_type = coord_type
        self.opts = opts

    class PysisCoordinateInterpolator:
        def __init__(self, initial_path):
            self.path = initial_path

    def prep_images(self,
                    num_images=None,
                    energy_evaluator=None,
                    base_images=None,
                    coord_type=None
                    ):
        from McUtils.ExternalPrograms import patch_pysis_logging
        patch_pysis_logging()

        from pysisyphus.Geometry import Geometry

        if base_images is None:
            base_images = super().generate(num_images=num_images)

        if coord_type is None:
            coord_type = self.coord_type

        if energy_evaluator is None:
            energy_evaluator = self._energy_evaluator
        if energy_evaluator is not None:
            base_images = [
                b.modify(energy_evaluator=energy_evaluator)
                for b in base_images
            ]

        geoms = [
            Geometry(
                mol.atoms,
                mol.coords,
                coord_type=coord_type
            )
            for mol in base_images
        ]
        for g,m in zip(geoms, base_images):
            g.set_calculator(m.get_energy_evaluator().to_pysis())

        return base_images, geoms

    def evaluate_profile_distances(self, profile: 'list[Molecule]', normalize=True):
        dists = [
            p1.get_rmsd(p2)
            for p1, p2 in zip(profile[:-1], profile[1:])
        ]
        d = np.cumsum([0] + dists)
        if normalize:
            d = d / d[-1]
        return d

    def evaluate_profile_energies(self, profile: 'list[Molecule]', energy_evaluator=None):
        if energy_evaluator is None:
            energy_evaluator = self._energy_evaluator

        return [
            p.modify(energy_evaluator=energy_evaluator).calculate_energy()
            for p in profile
        ]

    def generate(self,
                 num_images=None,
                 energy_evaluator=None,
                 base_images=None,
                 *,
                 method=None,
                 optimizer=None,
                 **opt_opts):

        from McUtils.ExternalPrograms import patch_pysis_logging, run_pysisyphus
        patch_pysis_logging()

        if method is None:
            method = self.default_method

        base_images, images = self.prep_images(
            num_images=num_images,
            energy_evaluator=energy_evaluator,
            base_images=base_images
        )

        opt_data = run_pysisyphus(
            energy_evaluator,
            method,
            images=images,
            optimizer=optimizer,
            return_logs=False,
            **(self.opts | opt_opts)
        )

        return [
            b.modify(coords=i.coords.reshape(-1, 3))
            for b,i in zip(base_images, images)
        ]

class PysisNEBGenerator(PysisyphusProfileGenerator):
    def __init__(self,
                 reactant_complex: Molecule,
                 product_complex: Molecule,
                 *,
                 energy_evaluator: EnergyEvaluator,
                 coordinate_interpolator='pysis',
                 num_images=10,
                 initial_image_positions=None,
                 internals=None,
                 max_displacement_step=None,
                 intermediates=None,
                 coord_type=None,
                 spring_constant=.1
    ):
        super().__init__(
            reactant_complex,
            product_complex,
            energy_evaluator=energy_evaluator,
            coordinate_interpolator=coordinate_interpolator,
            num_images=num_images,
            initial_image_positions=initial_image_positions,
            internals=internals,
            max_displacement_step=max_displacement_step,
            intermediates=intermediates,
            coord_type=coord_type
        )
        self.spring_constant = spring_constant

    def generate(self,
                 num_images=None,
                 k_min=None,
                 spring_constant=None,
                 energy_evaluator=None,
                 return_preopt=False,
                 base_images=None,
                 method='neb',
                 optimizer=None,
                 **opt_opts):
        if spring_constant is None:
            spring_constant = self.spring_constant
        if k_min is None:
            k_min = spring_constant
        return super().generate(
            num_images=num_images,
            k_min=k_min,
            energy_evaluator=energy_evaluator,
            base_images=base_images,
            method=method,
            optimizer=optimizer,
            **opt_opts
        )

@ProfileGenerator.register('pys-gsm')
class PysisGSMGenerator(PysisyphusProfileGenerator):
    default_method = 'gsm'

@ProfileGenerator.register('pys-fsm')
class PysisFSMGenerator(PysisyphusProfileGenerator):
    default_method = 'fsm'

@ProfileGenerator.register('pys-cos')
class PysisCOSGenerator(PysisyphusProfileGenerator):
    default_method = 'cos'

@ProfileGenerator.register('pys-string')
class PysisZTSGenerator(PysisyphusProfileGenerator):
    default_method = 'zts'

@ProfileGenerator.register('pys-dimer')
class PysisDimerGenerator(PysisyphusProfileGenerator):
    default_method = 'dimer'