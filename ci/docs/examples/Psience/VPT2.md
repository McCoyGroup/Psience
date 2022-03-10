
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
        3,
        logger=True
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

You can do the same by using the provided objects directly like

<div class="card in-out-block" markdown="1" id="Markdown_code">

```python
system = VPTSystem("HOH_freq.fchk")
states = VPTStateSpace.from_system_and_quanta(system, 3)
pt_opts = VPTSolverOptions(state_space_filters=states.get_filter("intensities"))
run_opts = VPTRuntimeOptions(logger=True)
runner = VPTRunner(system, states, runtime_options=run_opts, solver_options=pt_opts)
runner.print_tables()
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

