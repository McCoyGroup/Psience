
In the following we provide some basic examples.
More complex cases can be composed from the many settings provided in the Hamiltonian, solver, and runtime objects.

### A note on input formats

A very common use case for this package is to extract more information from a VPT2 calculation
performed by some electronic structure package. To support this use case the [`Molecule`](../Molecools/Molecule.md)
object can read in formatted checkpoint (`.fchk`) files and most of these examples will run off of those.

### Running Jobs

The way to run jobs with the least boilerplate is to use the `run_simple` method of the `VPTRunner` class.

Run second-order vibrational perturbation theory for all states with up to 
three quanta of excitation and print their infrared spectroscopic data

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

It is can do the same by using the driver objects directly, but it is less convenient.
If one is interested in doing further analysis on the data, it is also possible to run jobs directly from a `VPTAnalyzer` to capture the output, like

```python
analyzer = VPTAnalyzer.run_VPT("HOH_freq.fchk", 3)
``` 

The call structure is identical. The output is just routed through the analyzer.
The output tables can then be printed out using the analyzer's `print_output_tables` method and various other properties can be requested.
For the complete set of properties, check the [`VPTAnalyzer`](VPTAnalyzer.md) documentation.

### Internal Coordinates

It is possble to run the VPT in curvilinear internal coordinates. 
The most straightforward way to do this is to supply a Z-matrix (i.e. using polyspherical coordinates).
For this, we say

```python
VPTRunner.run_simple(
    ...,
    internals=zmat
    )
```
where `zmat` is an N-by-4 array of `int`s like
```python
[
    [atom1, bond1, angle1, dihed1], 
    [atom2, bond2, angle2, dihed2], 
    ...
]
```
where the `atom` specifies which of the atoms from the original set of Cartesian coordinates we're describing (`0` for the first atom, `1` for the second, etc.), `bond` specifies which atom it is bonded to, `angle` provides the third atom for making an angle, and `dihed` provides the dihedral angle.
It is important to note that the first bond distance, the first two angles, and the first three dihedrals aren't well defined and are generally used to store information about how the system is embedded in Cartesian space.

It is also possible to use functions of the internal coordinates.
To do this we supply 
```python
internals={
    'zmatrix': zmat,
    'conversion': conv,
    'inverse': inv
}
```

where `conv` looks like
```python
def conv(r, t, f, **kwargs):
    """
    Takes the bond lengths (`r`), angles `(t)` and dihedrals `(f)`,
    and returns new coordinates that are functions of these coordinates
    """
    ... # convert the coordinates
    return np.array([r, t, f])
```

and `inv` will take the output of `conv` and return the original Z-matrix/polyspherical coordinates.