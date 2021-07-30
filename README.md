# GEM

## Compiling
Compiling on a local machine with the GNU compilers requires installing LAPACK and OpenMPI as well as setting the following environment variables by exporting them (in ~/.bashrc for instance):
* LAPACK_LIB_DIR
* OMP_LIB_DIR
* OMP_INC_DIR

Compiling on NERSC/Cori requires loading the correct modules for the following compilers:</br>
**Note:** Use `module swap PrgEnv-intel PrgEnvl-cray` if you need to swap rather than load. This also works with `module swap impi openmpi` for instance.
* Cray:  module load PrgEnv-cray
* Intel: module load PrgEnv-intel impi
* GNU:   module load PrgEnv-gnu openmpi

## Debugging
Debugging flags for the compilers can be enabled by setting the `DEBUG` flag to 1 in the makefile. The makefile is source controlled though so maybe there is a better way to handle this that won't get tracked.
