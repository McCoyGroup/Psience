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
The main thing we need to do for scripts (and Jupyter) is add the line

```python
from Psience.VPT2 import *
```

and all of the helpers should be loaded.

## Anne Input Helpers

I've attached an object called `AnneInputHelpers` to the [`VPTRunner`](https://mccoygroup.github.io/Psience/Psience/VPT2/Runner/VPTRunner.html) object set up to make running VPT jobs a little bit simpler.
This object exists to enable people who aren't comfortable with python and standard python data loading to supply their custom job input data as text files.
It exposes one core function `run_anne_job`, that will take a folder with the appropriate files and run a VPT job.

The signature of the function [can be found on GitHub](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#:~:text=run_anne_job) and looks like

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
- `mode_selection` (`list[int]`): the subset of normal modes to use in the calculation as a `list` of `int`s corresponding to the desired modes. This can also be used to rearrange from freq. ordering to Herzberg by supplying a list like `[1, 2, 3, 4, 0, 5]` to swap some modes around.
- `basis_filters` (`list[dict]`): a list of filters to apply sequentially to the basis of states used in the PT, each filter can look like one of the following
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
