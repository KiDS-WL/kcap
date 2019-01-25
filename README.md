# KiDS Cosmology Analysis Pipeline

Pipeline for the cosmology analysis of KiDS 1000.

The pipeline is built on CosmoSIS, albeit a modified version that doesn't rely on environmental variables.

The different modules (CosmoSIS standard library, HMx, etc) are included as git subtree. Users don't have to worry about this detail but if you make changes to any of the modules it helps to structure your commits such that they only touch on one module at a time, such that these changes can be easily backported to the individual repositories.

## Installation

Clone the repository:
```
git clone git@github.com:KiDS-WL/kcap.git
cd kcap
```

It's strongly recommended to use some kind of encapsulated environment to install `kcap`, e.g., using `conda`. Here we assume that there is a anaconda installation available, that we need MPI support, and that we're on a machine with up-to-date GCC compilers. Notes on installations on macOS and details on how to set up things manually are [here](#installation-on-macos-and-other-details).

On machines with `module` support (e.g., cuillin), load the anaconda and openmpi modules first:
```
module load anaconda
module load openmpi
```
If there's no automated way to load these modules, make sure `conda` and MPI executables (`mpif90`) are on your `PATH`. For instructions on how to set up your own conda installation, see [Install conda](#install-conda).

Now set up the conda environment using the provided `conda_env.yaml` file:
```
conda env create -f conda_env.yaml
```
This creates a `kcap_env` environment that should have all the necessary dependencies. Activate the environment with `source activate kcap_test`.

We can now build kcap (which installs a standalone version of CosmoSIS):
```
python build.py
```

To test that everything is working, run the tests (todo...) examples in `examples/`:
```
mkdir examples/output
cosmosis examples/example_b.ini
```
For MPI:
```
mpirun -n 4 cosmosis --mpi examples/example_b.ini
```

To uninstall CosmoSIS (for example if you need to get the newest version), run `pip uninstall cosmosis_standalone`. To make a fresh installation of kcap, run `python build.py --clean`.

### Installation on macOS and other details

On macOS the current standard compiler suite is too old for building CosmoSIS. To alleviate this, install an up-to-date version of gcc (e.g., with `homebrew` by running `brew install gcc`) and set the environmental variables `CC`, `CXX`, (and `FC` and `MPIFC`, in case Fortran compilers other than `gfortran` and `mpif90` are needed). 

If no MPI support is required, run `python build.py --no-mpi`.

### Install conda

Get [miniconda](https://conda.io/en/master/miniconda.html). For example on a Linux machine: 
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

The installation will ask whether you want to add this conda installation to your `PATH` by adding some lines to your `.bashrc`. If you choose not to (for example because you don't want it to be the default installation), make sure the installation is accessible when building and running kcap.

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