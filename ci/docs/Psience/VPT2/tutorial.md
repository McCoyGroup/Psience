# Getting Started with VPT

I've set up this test environment so we can play with things more easily.
This kind of environment can easily be set up on Mox or in a container too
so that we can run stuffthere.

## Setting up the Environment

The most basic thing we need to do is make sure we have a python environment
we can work with.
To do that we'll set up what python calls a [virtual environment](https://docs.python.org/3.8/tutorial/venv.html)
which just means python makes a copy of its dependency directories so you can
install things without worrying about global conflicts.

First, we need to install Python 3.9+,
which just means going to the [python downloads page](https://www.python.org/downloads/)
and using the correct installer.

Once we have python installed, we can `cd` into this directory and create an
environment (which I'm calling `mcenv`) like

```lang-shell
python3.9 -m venv mcenv
```

Next, we'll activate this environment.

```lang-shell
. mcenv/bin/activate
```

Finally, we need to install the wrapper package for the development
version of my VPT stuff, which I chose to call `Psience`.

```lang-shell
pip install mccoygroup-psience
```

And that should be it!

When we need to update the package, we can do this with `pip`, too like

```lang-shell
pip install --upgrade mccoygroup-psience
```

## Using JupyterLab

Lately, I've been doing a lot of my work in JupyterLab, which is basically
python's version of a Mathematica-like interface.
It's a little bit clunky, but you get used to it.
I set up a demo notebook so we can see some examples of running jobs.
To use this, just run

```lang-shell
jupyter lab
```

and then use the browser window to open `VPT.ipynb`

It's worth checking out the [JupyterLab documentation](https://jupyterlab.readthedocs.io/en/stable/) as well.

## Running Scripts

Most of my calculations are run through python scripts.
Since the `VPT2` package is just a part of the larger `Psience` library, it is straightforward to create arbitrarily complicated scripts.
This is how I tend to run my results for papers, when the number of
systems and options for the runs changes frequently.

## Basic Examples

The simplest way to run jobs is through [`VPTRunner.run_simple`](https://mccoygroup.github.io/Psience/Psience/VPT2/Runner/VPTRunner.html).
The signature for that is

```python
VPTRunner.run_simple(
    system,
    states,
    target_property=None, 
    corrected_fundamental_frequencies=None, 
    calculate_intensities=True,
    **opts
    )
```

where the properties are as follows

- `system`: `list|str|Molecule`
    >the system spec, either as a `Molecule`, molecule spec (atoms, coords, opts) or a file to construct a `Molecule`
- `states`: `int|list`
    >the states to get corrections for either an `int` (up to that many quanta) or an explicit state list
- `target_property`: `str`
    >the target property to get corrections for (one of 'frequencies', 'intensities', 'wavefunctions')
- `corrected_fundamental_frequencies`: `Iterable[float]|None`
    >a set of fundamental frequencies to use to get new zero-order energies
- `calculate_intensities`: `bool default:True`
    >whether or not to calculate energies
- `opts`: `Any`
    >options that work for a `VPTSystem`, `VPTStateSpace`, `VPTRuntimeOptions`, `VPTSolverOptions`, or `VPTHamiltonianOptions` object which will be filtered automatically

This can be used quite simply if one has an `fchk` file containing a partial quartic expansion of a potential (e.g. from a `freq=Anh` job in Gaussian)

<div class="card in-out-block" markdown="1" id="Markdown_code">

```python
VPTRunner.run_simple(
        "HOH_freq.fchk",
        3 # up through three quanta of excitation
    )
```

<div class="card-body out-block" markdown="1">

```lang-none
...
                   Harmonic                  Anharmonic
State       Frequency    Intensity       Frequency    Intensity
  0 0 1    3937.52466     67.02051      3744.74223     64.17167
  0 1 0    3803.29960      4.14283      3621.97931      3.11401
  1 0 0    1622.30302     67.45626      1572.70734     68.32367
  0 0 2    7875.04932      0.00000      7391.41648      0.01483
  0 2 0    7606.59919      0.00000      7155.85397      0.31496
  2 0 0    3244.60604      0.00000      3117.39090      0.55473
  0 1 1    7740.82426      0.00000      7200.36337      2.20979
  1 0 1    5559.82768      0.00000      5294.37886      3.76254
  1 1 0    5425.60262      0.00000      5174.61359      0.06232
  0 0 3   11812.57398      0.00000     10940.02275      0.04985
  0 3 0   11409.89879      0.00000     10601.62396      0.00898
  3 0 0    4866.90906      0.00000      4634.05068      0.00350
  0 1 2   11678.34892      0.00000     10680.67944      0.00001
  1 0 2    9497.35234      0.00000      8917.98240      0.00333
  0 2 1   11544.12385      0.00000     10567.87984      0.08362
  2 0 1    7182.13070      0.00000      6815.99171      0.16303
  1 2 0    9228.90221      0.00000      8688.41518      0.00427
  2 1 0    7047.90564      0.00000      6699.22408      0.00661
  1 1 1    9363.12728      0.00000      8729.92693      0.09713
```

</div>
</div>

For more sophisticated calculations, many other flags are provided

### (from [`VPTStateSpace`](../VPT2/Runner/VPTStateSpace.md))

- `degeneracy_specs`: `'auto' | list | dict`
    >A specification of degeneracies, either as polyads, explicit groups of states, or parameters to a method. (see Details for more info)

### (from [`VPTHamiltonianOptions`](../VPT2/Runner/VPTHamiltonianOptions.md))

- `mode_selection`: `Iterable[int]|None`
    >the set of the supplied normal modes to do perturbation theory on (can also be used to rearrange modes to put them in ordering from Herzberg notation)
- `operator_coefficient_threshold`: `float|None`
    >the minimum size of a coefficient to keep when evaluating representation terms
- `pseudopotential_terms`: `Iterable[np.ndarray]`
    >explicit values for the psuedopotential terms
- `coriolis_terms`: `Iterable[np.ndarray]`
    >explicit values for the Coriolis terms
- `kinetic_terms`: `Iterable[np.ndarray]`
    >explicit values for the kinetic terms (e.g. from analytic models), same format as for the potential
- `potential_terms`: `Iterable[np.ndarray]`
    >explicit values for the potential terms (e.g. from analytic models), should be a list of tensors starting with the Hessian with each axis of length `nmodes`
- `include_pseudopotential`: `bool`
    >whether or not to include the pseudopotential/Watson term
- `include_coriolis_coupling`: `bool`
    >whether or not to include Coriolis coupling in Cartesian normal mode calculation
- `backpropagate_internals`: `bool`
    >whether or not to do Cartesian coordinate calculations with values backpropagated from internals


### (from [`VPTSolverOptions`](../VPT2/Runner/VPTSolverOptions.md))
- `check_overlap`: `bool default:True`
    >whether or not to ensure states are normalized in the VPT
- `zero_order_energy_corrections`: `dict`
    >energies to use for the zero-order states instead of the diagonal of `H(0)`
- `low_frequency_mode_cutoff`: `float (default:500 cm-1)`
    >the energy below which to consider a mode to be "low frequency"
- `state_space_filters`: `dict`
    >filters that can be used to cut down on the size of bases (see `VPTRunner.get_state_space_filter`)
- `order`: `int`
    >the order of perturbation theory to apply
- `expansion_order`: `int | dict`
    >the order to go to in the expansions of the perturbations, this can be supplied for different properties independently, like
    > ```python
    > expansion_order = {
    >  'potential':some_int,
    >  'kinetic':some_int,
    >  'dipole':some_int
    >  }
    >```
- `zero_element_warning`: `bool`
    >whether or not to warn if an element of the representations evaluated to zero (i.e. we wasted effort)

### (from [`VPTRuntimeOptions`](../VPT2/Runner/VPTRuntimeOptions.md))
- `results`: `str|Checkpointer|None default:None`
    >the `Checkpointer` to write corrections out to
- `logger`: `str|Logger|bool|None default:None`
    >the `Logger` object to use when logging the status of the calculation (`True` means log normally)
- `nondeg_hamiltonian_precision`: `int`
    >the precision with which to print out elements in the degenerate coupling Hamiltonians in the log file
- `matrix_element_threshold`: `float|None default:None`
    >the minimum size of matrix element to keep
- `operator_chunk_size`: `int|None default:None`
    >the number of representation matrix elements to calculate in at one time

## Anne Input Helpers

I've also attached an object called `AnneInputHelpers` to [`VPTRunner`](https://mccoygroup.github.io/Psience/Psience/VPT2/Runner/VPTRunner.html) object set up to make running VPT jobs a little bit simpler.
This object exists to enable people who aren't comfortable with python and standard python data loading to supply their custom job input data as text files.
It exposes one core function `run_anne_job`, that will take a folder with the appropriate files and run a VPT job.

The signature of the function (as of when this is being written) looks like

```python
VPTRunner.helpers.run_anne_job(
    base_dir,
    states=2,
    order=None,
    expansion_order=None,
    atoms_file='atom.dat',
    masses_file='mass.dat',
    coords_file='cart_ref.dat',
    modes_file='nm_int.dat',
    zmat_file='z_mat.dat',
    potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'),
    dipole_files=('lin_dip.dat', 'quad_dip.dat', "cub_dip.dat", "quart_dip.dat", 'quintic_dip.dat'),
    results_file=None,
    **opts
    )
```

The arguments are as follows
- `base_dir` (`str`): the directory where the files for running perturbation theory live
- `states` (`int` or `list`): the states to do PT for, an `int` means "do PT for all states with up through this many quanta
- `order` (`int`): the order to which to do PT, inferred automatically from the `expansion_order` if not provided
- `expansion_order` (`int` or `dict`): the order to which to take each part of the expansion. If not supplied, the order will be inferred from the potential and dipole files supplied. If an `int`, all terms are taken out to that order. Otherwise supply a `dict` like

```python
{
    'potential':int,
    'kinetic':int,
    'dipole':int,
}
```

- `atoms_file` (`str`): a file that contains the atomic numbers, like

```lang-none
Header Line
atom1 atom2 ... atomN
```

- `masses_file` (`str`): a file that contains the atomic masses, like

```lang-none
mass1 mass2 ... massN
```

Note that only one of `atoms` and `masses` needs to be supplied.
If both are provided, `masses` will take precedence.

- `modes_file` (`str`): a file that contains the frequencies of the modes as well as the transformation matrix from internals to normal modes (`L`) and the inverse transformation matrix (`L^-1`) like

```lang-none
freq_1 freq_2 freq_3 ... freq_m

L11 L12 L13 .. L21 L22 ... Lmn
K11 K12 K13 .. K21 K22 ... Knm
```

- `zmat_file` (`str`): a file containing the molecular Z-matrix for internal coordinate inputs
- `coords_file` (`str`): a file containing the molecular Cartesian coordinates, looking like

```lang-none
x1 y1 z1
x2 y2 z2
.  .  .
.  .  .
.  .  .
xN yN zN
```

- `potential_files` (`list[str]`): a list of files containing the potential tensors in dimensionless normal modes, starting with the Hessian and going up from there
- `dipole_files` (`list[str]`): a list of files containing the dipole tensors in dimensionless normal modes, starting with the linear dipole and going up from there

Tensor files are just lists of indices and corresponding values, like

```
i1_1 i1_2 i1_3  ... val1
i2_1 i2_2 i2_3  ... val2
.    .    .     .    .
.    .    .     .    .
.    .    .     .    .
iN_1 iN_2 iN_3  ... valN
```

For the `dipole_files`, the `x`/`y`/`z` axis index should be the left-most column

The job can of course be run _without_ these files, by simply using vanilla `numpy` and python to set up the job, but this has proven convenient.

### Extra Options

There are a very large number of knobs and levers one can turn when running a PT job.
Here are some of the standard ones

- `internals`: an internal coordinate spec, three forms of this are supported
  - `None`: use Cartesians
  - Z-matrix as a 4-column array of `int`s with each row being refs to `[atom, bond, angle, dihedral]` (where `atom` allows for rearrangement of the atoms)
  - `dict`: a spec with a Z-matrix and/or coordinate transform like

```python

{
    'zmatrix': [[atom1, bond1, angle1, dihed1], [atom2, bond2, angle2, dihed2], ...] or None,
    'conversion': 'a function to convert from Z-matrix coordinates to desired coordinates',
    'inverse': 'the inverse conversion'
}
```

- `coordinate_transformation` (`[conversion, inverse]`): a shorthand in `run_anne_job` that allows for just specifying a conversion function and inverse to apply to the Z-matrix coordinates
- `logger` (`str`): a file to print the job log info to
- `results` (`str`): a file to save results to (must have the extension `.hdf5` or `.json`)
- `degeneracy_specs`: The specification of degeneracies.

There are multiple possible values for this spec.
The simplest is to use the automatic approach, in which we supply a numeric type (`int`, `float`, etc.) to use as the `WFC` threshold.
The next simplest is to explicitly supply the groups we want, like

```python
[
    [ # the first resonant space
        state_1,
        state_2,
        state_3
    ],
    [ # the second
        state_5, state_11, ...
    ],
    ...
]
```

We can also supply pairs of relations for determining resonances, like

```python
[
    [state_1,  state_2], # A first relation
    [state_3,  state_4],  # another relation
    ...
]
```

To allow for extra options, you can also supply a `dict`. If you wanted to have a different `wfc_threshold` and you wanted to do the secondary resonant space splitting step with a very large threshold, you could do that by supplying

```python
{
    'wfc_threshold':.1,
    'energy_cutoff':1.0 # in Hartree
}
```

or you can explicitly add extra groups to the pairs of polyad rules by saying

```python
{
    'polyads':[
            [state_1,  state_2], # A first relation
            [state_3,  state_4],  # another relation
            ...
        ],
    'extra_groups': [
        [ # the first resonant space
            state_a,
            state_b,
            state_c
        ],
        [ # the second
            state_d, state_e, ...
        ],
        ...
    ]
}
```

This also allows us to define more resonance handling strategies.

The Martin Test is supported,
```python
{
    'martin_threshold':.1/219465, #in Hartree
}
```
As are total quanta vectors
```python
{
    'nT': [1, 1, 1, 0, 2, 2, 0] # e.g.
}
```
- `mode_selection` (`list[int]`): the subset of normal modes to use in the calculation as a `list` of `int`s corresponding to the desired modes (can also be used to rearrange from freq. ordering to Herzberg)
- `basis_postfilters` (`list[dict]`): a list of filters to apply sequentially to the basis of states used in the PT, each filter can look like one of the following
  - for excluding quanta

```python
{
    'max_quanta': [2, -1, 1, -1, ...] # the max number of quanta allowed in a given mode in the basis (-1 means infinity)
}
```

  - for excluding transitions

```python
{
    'excluded_transitions': [[0, 0, 1, 0, ...], [1, 0, 0, 0, ...], ...] # a set of transitions that are forbidden on the input states
}
```

  - for excluding transitions

```python
{
    'test': func # a function that takes the basis and tests if states should be allowed  
}
```
- `operator_coefficient_threshold` (`float` or `None`): the minimum coefficient size for a term in the expansion of the Hamiltonian (`None` means no threshold), jobs run faster with a larger threshold but the results are more approximate

See [the docs](https://mccoygroup.github.io/Psience/Psience/VPT2/) for the full list of options available to the [Hamiltonian](https://mccoygroup.github.io/Psience/Psience/VPT2/Runner/VPTHamiltonianOptions.html), [Solver](https://mccoygroup.github.io/Psience/Psience/VPT2/Runner/VPTSolverOptions.html), and [Runtime](https://mccoygroup.github.io/Psience/Psience/VPT2/Runner/VPTRuntimeOptions.html).

### Getting Expansions

We can also run a standard job with the `results` keyword to save the expansions (and wave functions and energies) to a file, which allows us to then extract the expansions.
Since it's easiest to run this all in a single directory, I added a little helper function `run_fchk_job`. This has the file signature


```python
VPTRunner.helpers.run_fchk_job(
    base_dir,
    states=2,
    fchk_file='fchk.fchk',
    zmat_file='z_mat.dat',
    results_file='output.hdf5',
    **opts
    )
```

which if given the file `fchk.fchk` in `base_dir` will run the job using the appropriate internals from `zmat_file` and save the outputs to `output.hdf5` and from there we can extract the potential and dipole expansions with the helpful functions `extract_potential` and `extract_dipole_expansion` like


```python
VPTRunner.helpers.extract_dipole_expansion(base_dir)
```

For example, we can run a job like this


```python
VPTRunner.helpers.run_fchk_job(
    "h2co/r",
    2,
    fchk_file='h2co.fchk'
)
```

and pull the expansions like this


```python
VPTRunner.helpers.extract_potential('h2co/r')
```

### Defining Custom Coordinate Systems

To use custom coordinate systems, we define a conversion function and an inverse


```python
def conv(r, t, f, **kwargs):
    """
    Takes the bond lengths (`r`), angles `(t)` and dihedrals `(f)`,
    and returns new coordinates that are functions of these coordinates
    """
    # ... means skip all indices until the specified ones (look up Numpy Fancy Indexing)
    cp1 = np.cos(f[..., 3])  # skip indices 0, 1, and 2 as these correspond to the embedding
    ct1 = np.cos(t[..., 2])  # skip 0 and 1 for embedding
    ct2 = np.cos(t[..., 3])
    st1 = np.sin(t[..., 2])
    st2 = np.sin(t[..., 3])
    f[..., 3] = np.arccos(st1 * st2 * cp1 + ct1 * ct2)
    return np.array([r, t, f])

def inv(r, t, f, **kwargs):
    cp1 = np.cos(f[..., 3])
    ct1 = np.cos(t[..., 2])
    ct2 = np.cos(t[..., 3])
    st1 = np.sin(t[..., 2])
    st2 = np.sin(t[..., 3])
    f[..., 3] = np.arccos((cp1 - ct1 * ct2) / (st1 * st2))
    return np.array([r, t, f])
```

then we add the flag `coordinate_transformation=[conv, inv]` to get the runner to use the transformations and inverse


```python
VPTRunner.helpers.run_anne_job(
    'nh3_2/r',
    2,
    order=2,
    degeneracy_specs='auto',
    coordinate_transformation=[conv, inv]
)
```