#!/bin/bash

tmin=0.50
tmax=300.00

python run_measure_cosebis_cats2stats.py -i /disk09/ma/Flinc/xi/xi_LOS1_nBins_5_Bin1_Bin1.ascii --cfoldername /disk09/ma/Flinc/cosebis/ -o nBins_1_Bin1_Bin1_lin_to_log --norm ./TLogsRootsAndNorms/Normalization_${tmin}-${tmax}.table -r ./TLogsRootsAndNorms/Root_${tmin}-${tmax}.table -b lin_to_log --thetamin ${tmin} --thetamax ${tmax} -n 5