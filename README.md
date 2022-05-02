# EnsembleMC: Ensemble Monte Carlo Framework


**EnsembleMC** is a C++ framework which implements an ensemble Monte Carlo code and uses 
the scalfmm library to calculate Coulomb interactions with the Fast Multipole Method.


## Requirements

  - CMake v3.10.0 or later
  - C++ compiler that supports
    - C++14 [compiler support list](http://en.cppreference.com/w/cpp/compiler_support)
    - [OpenMP](http://www.openmp.org/resources/openmp-compilers/)
  - ScalFMM C++ library: Installation instruction on https://gitlab.inria.fr/solverstack/ScalFMM

<!-- 
#The following are optional:
#
#  - [Doxygen](http://www.stack.nl/~dimitri/doxygen/) to build the documentation.
#  - An MPI implementation to build the distributed files.
#  - Custom BLAS, FFT implementations.
#  - [StarPU](http://starpu.gforge.inria.fr/) for the relevant FMM implementations.
-->

## Build EnsembleMC


mkdir build
cd build
# Use cmake
cmake .. 
make
```
