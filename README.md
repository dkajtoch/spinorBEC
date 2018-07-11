spinorBEC
==
(C) 2018 Dariusz Kajtoch

This library was developed to study quantum metrological properties of lowest energy/thermal quantum states in spinor F=1 Bose-Einstein condensates. Results were included in my PhD thesis and a paper [PhysRevA.97.023616](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.023616).

Installation
--
The package is written entirely in Fortran and requires `Intel MKL` library and `ARPACK`. Intel routines handle sparse matrices whereas ARPACK is used for efficient computation of lowest energy states.

Currently, there is only `Makefile` and `Make.inc` which build the static library `libspinor.a` and generate module files `*.mod`. For proper build edit the `Make.inc` file.

Structure
--
* `data_struct.f90`: contains basic definitions of data structures.
* `sparse.f90`: wrapper for ARPACK eigensolver designed for real upper-symmetric sparse matrices.
* `hamiltonian.f90`: definition of the sparse matrix Hamiltonian.
* `density_matrix.f90`: definition of the density matrices, computed based on the Hamiltonian.
* `measures.f90`: functions calculating various properties of quantum states.
