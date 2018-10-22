#!/bin/bash

# Fitting mode
mode=11

# Number in chain
n=10

# Simulations
sims+=('AGN_TUNED_nu0')
sims+=('AGN_7p6_nu0')
sims+=('AGN_8p0_nu0')

# Redshifts
zs+=('0.0')
zs+=('0.5')
zs+=('1.0')
#zs+=('2.0')

# executable
bin=./bin/HMx_fitting

# Loop and run fitting code
for sim in "${sims[@]}"; do
    for z in "${zs[@]}"; do
	outbase=fitting/${sim}_z${z}_n${n}_m${mode}
	$bin $mode $n $outbase $sim $z
    done
done
