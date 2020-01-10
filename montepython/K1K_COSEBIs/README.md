This folder contains the likelihood module for the KiDS-1000 (in short: K1K) COSEBIs measurements from [K1K et al. 2020 (arXiv:20XX.YYYY)](http://adsabs.harvard.edu/abs/2020arXiv20XXYYYY).
The module will be working 'out-of-the-box' within this KCAP setup and an additional install of [MontePython](https://github.com/brinckmann/montepython_public) and [CLASS](https://github.com/lesgourg/class_public) (version >= 2.8 including the HMcode module) setup. The required KiDS-1000 data files can be downloaded from the [KiDS science data webpage](http://kids.strw.leidenuniv.nl/sciencedata.php) and the parameter file for reproducing the fiducial run of [K1K et al. 2020 (arXiv:20XX.YYYY)](http://adsabs.harvard.edu/abs/2020arXiv20XXYYYY) is supplied in the subfolder `INPUT`.

Assuming that KCAP and MontePython (with CLASS version >= 2.8 including the HMcode module) are set up (we recommend to use nested sampling), please proceed as follows:

1) Copy `__init__.py` and `K1K_COSEBIs.data` from this folder into a folder named `K1K_COSEBIs` within `/your/path/to/montepython_public/montepython/likelihoods/`.

(you can rename the folder to whatever you like, but you must use this name then consistently for the whole likelihood which implies to rename the `*.data`-file, including the prefixes of the parameters defined in there, the name of the likelihood in the `__init__.py`-file and also in the `*.param`-file.)

2) Set the path to the data folder (i.e. `KiDS-1000_DATA_RELEASE` from the tarball available from the [KiDS science data webpage](http://kids.strw.leidenuniv.nl/sciencedata.php') in `K1K_COSEBIs.data` and modify parameters as you please (note that everything is set up to reproduce the fiducial run with `K1K_COSEBIs.param`).

3) Start your runs using e.g. the `K1K_COSEBIs.param` supplied in the subfolder `INPUT`.

4) If you publish your results based on using this likelihood, please cite [K1K et al. 2020 (arXiv:20XX.YYYY)](http://adsabs.harvard.edu/abs/2020arXiv20XXYYYY) and all further references for the KiDS-1000 data release (as listed on the [KiDS science data webpage](http://kids.strw.leidenuniv.nl/sciencedata.php)) and also all relevant references for KCAP, Monte Python and CLASS.

Refer to `mpirun_with_multinest.sh` or `mpirun_with_polychord.sh` within the subfolder `INPUT` for all MultiNest and PolyChord-related settings that were used for the fiducial runs.

WARNING: This likelihood only produces valid results for `\Omega_k = 0`, i.e. flat cosmologies!
