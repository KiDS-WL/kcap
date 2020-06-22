import subprocess
import os
import argparse
import sys
sys.path.append("modules/scale_cuts/")
from wrapper_twopoint import TwoPointWrapper
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free mocks')
    #parser.add_argument('--noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy mocks indexed by i with BEGIN <= i < END')
    #parser.add_argument('--multi-start', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create random staring point from multinest prior indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # Always create 3x2pt data vector for mocks
    # Configs know how to use only part of them properly
    run_type = "EE_nE_w"
    data_name_root = 'VD'
    KiDS_twopoint_file = f"../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits" ## Provide cov
    KiDS_twopoint_file_base = f"runs/methodology/data/noisefree_fiducial/base_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits" ## Provide mean

    # Create noiseless mock data vector
    if args.noise_free:
        ## Scale by variable depth ratio
        VD_ratio_file = 'data2/mockFootprint/buceros/MFP_for_others/varDepth_ratio.npy'

        TP_cov     = TwoPointWrapper.from_fits(KiDS_twopoint_file, covmat_name='COVMAT')
        mean_base  = TwoPointWrapper.from_fits(KiDS_twopoint_file_base, covmat_name='COVMAT').makeMeanVector()
        ratio_mean = np.load(VD_ratio_file)
        TP_cov.replaceMeanVector(mean_base * ratio_mean)

        output_dir = "runs/methodology/data/noisefree_fiducial/"
        data_name = f"{data_name_root}_{run_type}"
        
        KiDS_twopoint_file_save  = f"{output_dir}{data_name}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits"
        os.makedirs(f'{output_dir}{data_name}/data/KiDS/', exist_ok=True)
        TP_cov.to_fits(KiDS_twopoint_file_save, overwrite=True)

