[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7382685.svg)](https://doi.org/10.5281/zenodo.7382685)

## What is optimwidths?

optimwidths is a tool to calculate the Gaussian width parameters used in *ab initio* multiple spawning (AIMS). It essentially implements a version of the width optimization algorithm for AIMS [^1] first described in the supplementary information of Ref. [^2] At the moment it only accepts Gaussian16 output files from a force calculation, but in the future more codes will be supported. 

## Installing and using otpimwidths

The only dependencies of optimwidths are the `BLAS` and `LAPACK` libraries, and a Fortran compiler. 
To compile the code go into the optimwidths folder, change the compiler and 
library flags in `make.vars` to the ones you want and type `make`. 

The binary `optimwidths.x` requires an outputfile generated from a Gaussian16 frequency calculation and an input file `input.in` that contains
the following parameters:

|  **Input variable name**  |       **Description**           |   **Possible options** | 
|:--------------------------|:--------------------------------|:----------------------:|
| `nrAtoms`                 |   Number of atoms               |      `integer`         |
| `xtol`                    |   Simplex convergence treshold  |       `real`           |
| `ftol`                    |   Function convergence treshold |       `real`           |
| `maxiter`                 |   Maximum number of iterations  |      `integer`         |
| `constrained`             |   Constrained optimazation?     |      `logical`         |
| `prog`                    |   Electronic structure program  |      `gaussian`        |
| `fName`                   |   Name of output file           |        `str`           |



[^1]: [ Thompson, A. L.; Punwong, C.; Mart&#237;nez, T. J. *Chem. Phys.* 2010, 370, 70-77 ]( https://doi.org/10.1016/j.chemphys.2010.03.020 )
[^2]: [ Esch, M. P.; Shu, Y.; Levine, B. G. *J. Phys. Chem. A* 2019, 123, 13, 2661–2673 ]( https://doi.org/10.1021/acs.jpca.9b00952 )


