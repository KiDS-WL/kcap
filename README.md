# KiDS Cosmology Analysis Pipeline

Pipeline for the cosmology analysis of KiDS 1000.

The pipeline is built on CosmoSIS, albeit a modified version that doesn't rely on environmental variables.

The different modules (CosmoSIS standard library, HMx, etc) are included as git subtree. Users don't have to worry about this detail but if you make changes to any of the modules it helps to structure your commits such that they only touch on one module at a time, such that these changes can be easily backported to the individual repositories.

## Installation

It's strongly recommended to use some kind of encapsulated environment to install `kcap`, e.g., using `conda`.

To get the pipeline working 
```
python build.py
```
should be sufficient. This installs the standalone CosmoSIS version with `pip` and builds the CosmoSIS standard library and other modules.

If MPI support is required, set the environmental variable `MPIFC` before building CosmoSIS. This requires the python package `mpi4py`. The conda version seems to be broken at this point. Using pip (i.e., `pip install mpi4py`) is working however.

To test that everything is working, run the tests (todo...) examples in `examples/`:
```
cosmosis examples/example_b.ini
```
For MPI:
```
mpirun -n 4 cosmosis --mpi examples/example_b.ini
```