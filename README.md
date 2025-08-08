# TheGeneralizedSBM.jl: A julia-based implementation of the Generalized Shifted Boundary Method using Gridap.jl

[![DOI](https://zenodo.org/badge/1029609698.svg)](https://doi.org/10.5281/zenodo.16764244)

## Abstract
The Generalized Shifted Boundary Method (GSBM) is a geometry-agnostic finite element framework for solving PDEs on complex, parametrized, and evolving domains.  This repository provides the reference implementation in [julia](https://julialang.org/) using the [Gridap.jl](https://github.com/gridap/Gridap.jl) ecosystem, with examples for the Poisson equation, linear elasticity, and transient Stokes flow.

## What is GSBM?

The **Generalized Shifted Boundary Method (GSBM)** is a geometry-agnostic finite element formulation that:
- Generalizes the Shifted Boundary Method (SBM) and the Weighted Shifted Boundary Method (WSBM)
- Avoids the need for redefinition of integration domains, FE spaces and geometry-specific data structures
- Is robust for problems with topological changes and moving domains

## Installation
The installation of this package is straight forward using the [Julia's package manager](https://julialang.github.io/Pkg.jl/v1/). Open the Julia REPL, type `]` to enter package mode, and install as follows
```julia
pkg> add https://github.com/CMOE-TUDelft/TheGeneralizedSBM.jl
```
or
```julia
using Pkg; Pkg.add(url="https://github.com/CMOE-TUDelft/TheGeneralizedSBM.jl")
```

## Usage
To run all the test cases in the GSBM journal paper do:
```julia
using TheGeneralizedSBM
run_tests("all")
```

To run only a specific test, for example the Poisson manufactured solution test, do:
```julia
using TheGeneralizedSBM
run_tests("ManufacturedSolutionsPoisson.jl")
```
All the tests in `scripts` folder can be executed.

After execution, the data will be stored in the respective folder `data/<problem-type>/<test-name>`. If the flag to generate VTK files is active, the VTK output will be stored in `data/<problem-type>/<test-name>/VTK_<geometry-name>`. The plots shown in the manuscript are stored in `plots/<problem-type>/<test-name>`.

This repository uses DrWatson package, the data will only be generated the first time the tests are executed. If the data is already stored, the scripts will only regenerate the figures.


## Authors
- [Oriol Colomés](www.oriolcolomes.com)$^*$
- [Jan Modderman]()$^*$

$^*$ Faculty of Civil Engineering and Geosciences, Delft University of Technology, Delft, 2628 CN, The Netherlands

