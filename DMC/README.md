# RynDMC

A package for doing DMC in a black-box, flexible way, implementing the tricks and utilities we use as a group

Definitely still in the very beta of beta stages, but will build out piece by piece

## Roadmap

### DMC

- [X] Implement an AbstractDMC class to be subclassed by smarter implementations
- [ ] Implement a DiscreteWeightingDMC subclass to handle discrete weighting DMC simulations
- [ ] Implement a ContinuousWeightingDMC subclass to handle continuous weighting DMC simulations
- [ ] Implement an ImportanceSampledDMC subclass of ContinuousWeightingDMC to add in importance sampling

### Potential Handlers

- [X] Create C-level potential handler library that can take just a pointer to a potential
- [ ] Write better communication paradigm for doing more work at the C level before passing back to python
- [ ] Create multiprocessing based potential handler for automatic parallelization of potential calls
- [ ] Write little wrapper to Fortan and C-types based hook-ins for simple potentials

### Customizations

- [ ] Add standard job config template for setting up a simulation
- [ ] Add 

## Dependencies

At the moment this depends on the [McUtils](https://github.com/McCoyGroup/McUtils) and [Coordinerds](https://github.com/McCoyGroup/Coordinerds)
These might find their way in a git module level dependencies, though.