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
should be sufficient. On macOS the current standard compiler suite is too old for building CosmoSIS. To alleviate this, install an up-to-date version of gcc (e.g., with `homebrew`) and set the environmental variables `CC`, `CXX`, (and `FC`, in case of a Fortran compiler other than `gfortran` is needed). This installs the standalone CosmoSIS version with `pip` and builds the CosmoSIS standard library and other modules.

If MPI support is required, set the environmental variable `MPIFC` before building CosmoSIS. This requires the python package `mpi4py`. The conda version seems to be broken at this point. Using pip (i.e., `pip install mpi4py`) is working however.

To test that everything is working, run the tests (todo...) examples in `examples/`:
```
cosmosis examples/example_b.ini
```
For MPI:
```
mpirun -n 4 cosmosis --mpi examples/example_b.ini
```

## Development

The different modules are organised as git subtrees.

### Pull updates from a specific module repository

Using `git subtree`:
```
git subtree pull --prefix=HMx --squash HMx-remote kcap
```
where `HMx` is the module to be updated, `HMx-remote` is the remote for the module, and `kcap` the remote branch. The option `--squash` collapses the histories of the modules. Especially for the CosmoSIS standard library this is useful, since we don't want its whole history in kcap.

Using the `subtree` merge strategy:
```
git merge -s subtree --squash --allow-unrelated-histories HMx-remote/kcap
```
The `--allow-unrelated-histories` seems to be necessary (probably because of the `--squash` option used earlier).

### Push changes of modules in the kcap repository to the module repository
Using `git subtree`:
```
git subtree push --prefix=HMx HMx-remote remote_branch
```
where `HMx` is again the module to has been updated, `HMx-remote` is the remote for the module, and `remote_branch` the remote branch that the updates will get push to.

Using the `subtree` merge strategy:
```
git checkout -b backport HMx-remote/kcap
git cherry-pick -x --strategy=subtree commits_to_push
git push
```
This is a bit more involved but allows for more control. First create a branch (`backport`) tracking the remote branch that is being targeted (`HMx-remote/kcap`). Then cherry pick the commits (`commits_to_push`) that should be pushed to the module repository.