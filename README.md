# KiDS Cosmology Analysis Pipeline

This repository contains the cosmology inference pipeline that was used in the KiDS-1000 analyses:
 - Methodology: [Joachimi, Lin, Asgari, Tröster, Heymans et al. 2021](https://arxiv.org/abs/2007.01844)
 - Cosmic shear: [Asgari, Lin, Joachimi et al. 2021](https://arxiv.org/abs/2007.15633)
 - 3x2pt: [Heymans, Tröster et al. 2021](https://arxiv.org/abs/2007.15632)
 - Beyond flat ΛCDM: [Tröster et al. 2021](https://arxiv.org/abs/2010.16416)

The pipeline is built on CosmoSIS, albeit a [modified version](https://bitbucket.org/tilmantroester/cosmosis/src/kcap/) that is `pip`-installable and doesn't rely on environmental variables.

A MontePython likelihood that wraps the kcap functionality can be found at [here](https://github.com/BStoelzner/KiDS-1000_MontePython_likelihood). 
Note that the standard version of MontePython does not support non-flat priors yet, which is a problem for samplers that distiguish between likelihood and prior (such as MultiNest and PolyChord). 
A version that supports Gaussian priors with MultiNest can be found [here](https://github.com/BStoelzner/montepython_public/tree/gaussian_prior).

The different modules (CosmoSIS standard library, etc) are included as git subtree. Users don't have to worry about this detail but if you make changes to any of the modules it helps to structure your commits such that they only touch on one module at a time, such that these changes can be easily backported to the individual repositories.

For a fiducial KV450 setup, have a look at `runs/config/KV450_fiducial.ini`.


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
If there's no automated way to load these modules, make sure `conda` and MPI executables (`mpif90`) are on your `PATH`. For instructions on how to set up your own conda installation, see [Install conda](#install-conda). If you're using your own anaconda installation, don't load the module as well, as this just causes conflicts.

Now set up the conda environment using the provided `conda_env.yaml` file:
```
conda env create -f conda_env.yaml
```
This creates a `kcap_env` environment that should have all the necessary dependencies. Activate the environment with `source activate kcap_env`.

We need to install CAMB because we use the new python interface for it. If `kcap` is to be used on a local machine, `pip install camb` is all there is to do. On a heterogenous cluster like `cuillin`, we need to build CAMB ourselves, however. To do so, run
```
git clone --recursive git@github.com:cmbant/CAMB.git
cd CAMB
python setup.py build_cluster
python setup.py install
```

We can now build kcap (which installs a standalone version of CosmoSIS):
```
python build.py
```

To uninstall CosmoSIS (for example if you need to get the newest version), run `pip uninstall cosmosis_standalone`. To make a fresh installation of kcap, run `python build.py --clean`.

### Installation on macOS and other details

The default macOS compilers are supported now but `gfortran` still needs to be installed. This can be done with `homebrew` by running `brew install gcc`. Note that `gcc 9.2` seems to be incompatible with the `PolyChord` samplers included in cosmosis, so use a different version (e.g., 9.1).

If no MPI support is required, run `python build.py --no-mpi`.

### Install conda

Get [miniconda](https://conda.io/en/master/miniconda.html). For example on a Linux machine: 
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

The installation will ask whether you want to add this conda installation to your `PATH` by adding some lines to your `.bashrc`. If you choose not to (for example because you don't want it to be the default installation), make sure the installation is accessible when building and running kcap.

## Usage

Make sure `conda` and MPI are accessible (e.g., by running `module load anaconda` and `module load openmpi`) and that `kcap_env` is activated (`source activate kcap_env`).
To test that everything is working, run the tests (todo...) and some of the configs in `runs/config`:
```
mkdir runs/output
cosmosis runs/config/KV450_no_sys.ini
```
For MPI:
```
mpirun -n 4 cosmosis --mpi runs/config/KV450_no_sys.ini
```

## Repository structure

- `cosebis`, `cosmosis-standard-library`: git subtrees that mirror the `kcap` branches of https://bitbucket.org/marika_a/cosebis_cosmosis and https://bitbucket.org/tilmantroester/cosmosis-standard-library, respectively.
- `data`: collection of data files. Does not include the KiDS-1000 data. KiDS-1000 data can be found [here](https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1) or a select subset in the `runs` directory.
- `kcap`: Python package for the pipeline. Deprecated.
- `modules`: non-trivial CosmoSIS modules used in kcap.
- `montepython`: old MontePython likelihood files. Deprecated in favour of https://github.com/BStoelzner/KiDS-1000_MontePython_likelihood.
- `runs`: collection of configurations and scripts to generate said configurations.
- `tests`: some tests to check consistency between kcap, CCL, and MontePython.
- `utils`: collection of tools and simple CosmoSIS modules used in kcap.


## Development

The different modules are organised as git subtrees.

### Pull updates from a specific module repository

Using `git subtree`:
```
git subtree pull --prefix=cosebis --squash cosebis-remote kcap
```
where `cosebis` is the module to be updated, `cosebis-remote` is the remote for the module, and `kcap` the remote branch. The option `--squash` collapses the histories of the modules. Especially for the CosmoSIS standard library this is useful, since we don't want its whole history in kcap.

Using the `subtree` merge strategy:
```
git merge -s subtree --squash --allow-unrelated-histories cosebis-remote/kcap
```
The `--allow-unrelated-histories` seems to be necessary (probably because of the `--squash` option used earlier).

### Push changes of modules in the kcap repository to the module repository
Using `git subtree`:
```
git subtree push --prefix=cosebis cosebis-remote remote_branch
```
where `cosebis` is again the module to has been updated, `cosebis-remote` is the remote for the module, and `remote_branch` the remote branch that the updates will get push to.

Using the `subtree` merge strategy:
```
git checkout -b backport cosebis-remote/kcap
git cherry-pick -x --strategy=subtree commits_to_push
git push
```
This is a bit more involved but allows for more control. First create a branch (`backport`) tracking the remote branch that is being targeted (`cosebis-remote/kcap`). Then cherry pick the commits (`commits_to_push`) that should be pushed to the module repository.
