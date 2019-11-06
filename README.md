# SG-FLEXI

[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000)]()

SG-FLEXI is a framework for solving the uncertain compressible 
Navier-Stokes and Euler equations using a Stochastic Galerkin 
expansion in the stochastic space and the Discontinuous Galerkin 
Spectral Element Method in the physical space. This allows for a
high order in all dimensions. 

SG-FLEXI can be used on fully unstructured hexahedral meshes.
The solver is parallelized very efficiently and scales up
to hundreds of thousand cores.

# FLEXI

SG-FLEXI is based on the high-order flow solver FLEXI.

FLEXI has been developed by the [Numerics Research Group (NRG)][nrg]
lead by Prof. Claus-Dieter Munz at the Institute of Aerodynamics
and Gasdynamics at the University of Stuttgart, Germany.

It is open source and can be found on the [FLEXI homepage][flexi].
There is extensive documentation for FLEXI, including a user guide 
and tutorials. For a deeper understanding of SG-FLEXI, we refer to 
the FLEXI tutorials. 

# License
FLEXI is Copyright (C) 2019, Prof. Claus-Dieter Munz and is 
released under the terms of the
GNU General Public License v3.0. For the full license terms see
the included license file [LICENSE.md](LICENSE).

## Citation
This is a scientific project. If you use SG-FLEXI for publications or
presentations in science, please support the project by citing
our publication given in [REFERENCE.md](REFERENCE.md).

## List of Contributors
Several people have worked on and with SG-FLEXI over the last years.
We would like to thank all these [CONTRIBUTORS.md](contributors)
for their efforts they spent on building SG-FLEXI.

## Funding
This code was developed in the course of the project “SEAL” funded by the Baden-Württemberg Stiftung. 

# Installation

For installation instruction see [INSTALL.md](INSTALL.md).

In case you have question regarding FLEXI, want to report bugs
or contribute to the project you can use the mailing list
<flexi@listserv.uni-stuttgart.de>.
You can also subscribe to the mailing list [here][list].

## Used libraries

FLEXI uses several external libraries as well as auxilliary functions from open source projects, including:
* [HDF5](https://www.hdfgroup.org/)
* [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](http://www.netlib.org/lapack/)
* [PAPI](http://icl.cs.utk.edu/papi/)
* [FOUL](http://foul.sourceforge.net/)
* [OpenMP](http://www.openmp.org/)
* [FFTW](http://www.fftw.org/)

[nrg]:  https://www.iag.uni-stuttgart.de/en/working-groups/numerical-methods/
[flexi]: https://www.flexi-project.org/
[list]: https://listserv.uni-stuttgart.de/mailman/listinfo/flexi

# Branches 

Please note that there are two branches of SG-FLEXI. Next to the master branch, there is a branch named viscosity, where the viscosity is uncertain and depends on the random input.
