# Installation
## Create environment and install dependencies
```
mamba create -n kcapv2
mamba activate kcapv2
```
```
mamba install python cosmosis astropy fast-pt camb
```

## Download cosmosis standard library and data
```
git clone https://github.com/joezuntz/cosmosis-standard-library
git clone https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1
```
No need to compile anything!

## Get kcap v2 and some the python bandpower/cosebis code
```
git clone git@github.com:KiDS-WL/kcap kcapv2
cd kcapv2
git checkout v2
```

```
cd 2pt_transforms
pip install -e .
```

## Run some chains
```
cd tests/kcapv2_vs_kcapv1/bp/kcapv2
cosmosis pipeline.ini
```