# PyVPT

The Python Vibrational Perturbation Theory package (PyVPT) is intended to make it easy to do vibrational perturbation theory in internal coordinates in python.

## Dependencies

A bunch of packages have been built out in python to support this effort (and other research too)
The list of them will change over time so it doesn't make much sense to keep a running tally here, but they're all [here](https://github.com/McCoyGroup).

## Road Map

### Coordinate Basics
 
 - [x] Set up basic project structure
 - [x] Create coordinate system representation
 - [x] Create coordinate transformation system
 - [x] Create coordinate conversion system
 - [x] Implement Cartesian <-> ZMatrix conversions
 - [ ] Generalize conversion code?
 - [ ] Implement Cartesian <-> Spherical conversions
 - [ ] Implement ZMatrix <-> Spherical conversions?

### Molecules
 - [x] Write efficient Gaussian Log File importer
 - [x] Write efficient Gaussian fchk file importer
 - [x] Write convenience type representing molecular modes
 - [x] Create efficient importer of all relevant molecule data from Gaussian job data
 - [ ] Create efficient constructor of all molecule data from scan (for use with psi4)
 
### Expansions
 - [x] Implement arbitrary order finite differencing (higher-dimensional extension needs testing still)
 - [ ] Provide utility type representing function over a coordinate set (useful for potentials)
 - [ ] Implement actual expression for a fourth-order potential expansion (make C-compileable...?)

### Perturbation Theory
 - [ ] Get potential expansions with Gaussian potential / force consts
 - [ ] Write code to handle directly from potential scan
 - [ ] Implement proper exprs for perturbation expansions
 
### Wavefunctions / Spectra
 - [x] Write flexible Wavefunction metaclass that can be subclassed for DVR, DMC, and VPT2
 - [ ] Write Spectrum class that can plot and analyze itself in intuitive ways like I already have in Mathematica (got the plotting...)

## Suggestions / Issues / Warnings

This is really not intended for outside use, but if there are features that should be included they can be suggested on the repo [issues page](https://github.com/McCoyGroup/PyVPT/issues)

