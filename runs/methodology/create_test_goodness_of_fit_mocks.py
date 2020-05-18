import subprocess
import os
import argparse
import sys
sys.path.append("modules/scale_cuts/")
from wrapper_twopoint import TwoPointWrapper
import numpy as np
import scipy.stats as stats


def sample_random_from_multinest(data, wgt):
    kernel = stats.gaussian_kde(data, 'silverman', weights=wgt)
    random_start = kernel.resample(1).flatten()
    
    priorDict = {
        'cosmological_parameters': {
            'omch2': [0.051, 0.2], 
            'ombh2': [0.019, 0.026], 
            'h0': [0.64, 0.82], 
            'n_s': [0.84, 1.1], 
            'S_8_input': [0.1, 1.3]
        }, 

        'halo_model_parameters': {
            'A': [2.0, 3.13]
        },

        'intrinsic_alignment_parameters': {
            'A': [-6.0, 6.0]
        },

        'nofz_shifts': {
            'p_1': [-5.0, 5.0], 
            'p_2': [-5.0, 5.0], 
            'p_3': [-5.0, 5.0], 
            'p_4': [-5.0, 5.0], 
            'p_5': [-5.0, 5.0]
        }, 

        'bias_parameters': {
            'b1_bin_1': [0.5, 9.0], 
            'b2_bin_1': [-4.0, 8.0], 
            'gamma3_bin_1': [-8.0, 8.0], 
            'a_vir_bin_1': [0.000001, 12.0], ## 1e-6 instead of 0 due to a bug
            'b1_bin_2': [0.5, 9.0], 
            'b2_bin_2': [-4.0, 8.0], 
            'gamma3_bin_2': [-8.0, 8.0], 
            'a_vir_bin_2': [0.000001, 12.0] ## 1e-6 instead of 0 due to a bug
        }
    }

    lower = []
    upper = []
    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            lower.append(value2[0])
            upper.append(value2[1])
    lower = np.array(lower)
    upper = np.array(upper)

    while np.any(random_start <= lower) or np.any(upper <= random_start):
        random_start = kernel.resample(1).flatten()

    starting_point_settings = []
    ind = 0

    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            p = random_start[ind]
            starting_point_settings.append("--set-parameters")
            starting_point_settings.append(key1)
            starting_point_settings.append(key2)
            starting_point_settings.append(f'{value2[0]} {p} {value2[1]}')
            ind += 1

    return starting_point_settings


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free mocks')
    parser.add_argument('--noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy mocks indexed by i with BEGIN <= i < END')
    parser.add_argument('--multi-start', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create random staring point from multinest prior indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    if args.noise_range[1] > args.noise_range[0]:
        print('You want to create noisy mocks.')
        print('Note that noisy mocks for test_goodness_of_fit require those for main_chains.')
        print('Make sure you have already created them.')

    script = "utils/run_kcap.py"

    # Always create 3x2pt data vector for mocks
    # Configs know how to use only part of them properly
    run_type = "EE_nE_w"
    data_name_root_2 = 'sabotaged2'
    data_name_root_5 = 'sabotaged5'
    data_name_root_S = 'sabotagedS'
    KiDS_twopoint_file = f"../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits"
    dz_cov_file = "data/KV450/nofz/SOM_cov_multiplied_diagonal.asc"

    # Create noiseless mock data vector
    if args.noise_free:
        output_dir = "runs/methodology/data/noisefree_fiducial/"
        data_name = f"{data_name_root_2}_{run_type}"
        dz_sabotage_settings_2 = [
            "--set-parameters", "nofz_shifts", "p_3", "2.0", ## This means 2-sigma.
            "--set-parameters", "nofz_shifts", "p_4", "2.0"
        ]

        cmd = ["--create-mocks", "--noiseless-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--dz-covariance-file", dz_cov_file,
                "--run-name", data_name,
                "--run-type", run_type,
                *dz_sabotage_settings_2]
        subprocess.run(["python", script] + cmd, check=True)

        data_name = f"{data_name_root_5}_{run_type}"
        dz_sabotage_settings_5 = [
            "--set-parameters", "nofz_shifts", "p_3", "5.0", ## This means 5-sigma.
            "--set-parameters", "nofz_shifts", "p_4", "5.0"
        ]

        cmd = ["--create-mocks", "--noiseless-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--dz-covariance-file", dz_cov_file,
                "--run-name", data_name,
                "--run-type", run_type,
                *dz_sabotage_settings_5]
        subprocess.run(["python", script] + cmd, check=True)

        data_name = f"{data_name_root_S}_{run_type}"
        dz_sabotage_settings_S = [
            "--set-parameters", "nofz_shifts", "p_1", "-2.0", ## This means 2-sigma.
            "--set-parameters", "nofz_shifts", "p_2", "2.0",
            "--set-parameters", "nofz_shifts", "p_3", "-2.0",
            "--set-parameters", "nofz_shifts", "p_4", "2.0",
            "--set-parameters", "nofz_shifts", "p_5", "-2.0"
        ]

        cmd = ["--create-mocks", "--noiseless-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--dz-covariance-file", dz_cov_file,
                "--run-name", data_name,
                "--run-type", run_type,
                *dz_sabotage_settings_S]
        subprocess.run(["python", script] + cmd, check=True)

    # Create noisy mocks
    # We want to use the same noisy mocks as `main_chains`
    # So we create new files by calculating the difference between noise-free data mean
    # and adding it to old mocks
    output_dir = "runs/methodology/data/noisy_fiducial/"
    KiDS_twopoint_file_base = f"runs/methodology/data/noisefree_fiducial/base_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits"
    KiDS_twopoint_file_sabo = f"runs/methodology/data/noisefree_fiducial/{data_name_root_2}_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits" 

    TP_cov     = TwoPointWrapper.from_fits(KiDS_twopoint_file, covmat_name='COVMAT')
    mean_base  = TwoPointWrapper.from_fits(KiDS_twopoint_file_base, covmat_name='COVMAT').makeMeanVector()
    mean_sabo  = TwoPointWrapper.from_fits(KiDS_twopoint_file_sabo, covmat_name='COVMAT').makeMeanVector()
    delta_mean = mean_sabo - mean_base

    for i in range(args.noise_range[0], args.noise_range[1]):
        data_name = f"{data_name_root_2}_{i}_{run_type}"

        KiDS_twopoint_file_noisy = f"{output_dir}base_{i}_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"
        KiDS_twopoint_file_save  = f"{output_dir}{data_name}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"

        os.makedirs(f'{output_dir}{data_name}/data/KiDS/', exist_ok=True)
        TP_mean = TwoPointWrapper.from_fits(KiDS_twopoint_file_noisy)
        mean    = TP_mean.makeMeanVector()
        TP_cov.replaceMeanVector(mean + delta_mean)
        TP_cov.to_fits(KiDS_twopoint_file_save, overwrite=True)

    output_dir = "runs/methodology/data/noisy_fiducial/"
    KiDS_twopoint_file_base = f"runs/methodology/data/noisefree_fiducial/base_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits"
    KiDS_twopoint_file_sabo = f"runs/methodology/data/noisefree_fiducial/{data_name_root_5}_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits" 

    TP_cov     = TwoPointWrapper.from_fits(KiDS_twopoint_file, covmat_name='COVMAT')
    mean_base  = TwoPointWrapper.from_fits(KiDS_twopoint_file_base, covmat_name='COVMAT').makeMeanVector()
    mean_sabo  = TwoPointWrapper.from_fits(KiDS_twopoint_file_sabo, covmat_name='COVMAT').makeMeanVector()
    delta_mean = mean_sabo - mean_base

    for i in range(args.noise_range[0], args.noise_range[1]):
        data_name = f"{data_name_root_5}_{i}_{run_type}"

        KiDS_twopoint_file_noisy = f"{output_dir}base_{i}_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"
        KiDS_twopoint_file_save  = f"{output_dir}{data_name}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"

        os.makedirs(f'{output_dir}{data_name}/data/KiDS/', exist_ok=True)
        TP_mean = TwoPointWrapper.from_fits(KiDS_twopoint_file_noisy)
        mean    = TP_mean.makeMeanVector()
        TP_cov.replaceMeanVector(mean + delta_mean)
        TP_cov.to_fits(KiDS_twopoint_file_save, overwrite=True)

    output_dir = "runs/methodology/data/noisy_fiducial/"
    KiDS_twopoint_file_base = f"runs/methodology/data/noisefree_fiducial/base_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits"
    KiDS_twopoint_file_sabo = f"runs/methodology/data/noisefree_fiducial/{data_name_root_S}_{run_type}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits" 

    TP_cov     = TwoPointWrapper.from_fits(KiDS_twopoint_file, covmat_name='COVMAT')
    mean_base  = TwoPointWrapper.from_fits(KiDS_twopoint_file_base, covmat_name='COVMAT').makeMeanVector()
    mean_sabo  = TwoPointWrapper.from_fits(KiDS_twopoint_file_sabo, covmat_name='COVMAT').makeMeanVector()
    delta_mean = mean_sabo - mean_base

    for i in range(args.noise_range[0], args.noise_range[1]):
        data_name = f"{data_name_root_S}_{i}_{run_type}"

        KiDS_twopoint_file_noisy = f"{output_dir}base_{i}_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"
        KiDS_twopoint_file_save  = f"{output_dir}{data_name}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"

        os.makedirs(f'{output_dir}{data_name}/data/KiDS/', exist_ok=True)
        TP_mean = TwoPointWrapper.from_fits(KiDS_twopoint_file_noisy)
        mean    = TP_mean.makeMeanVector()
        TP_cov.replaceMeanVector(mean + delta_mean)
        TP_cov.to_fits(KiDS_twopoint_file_save, overwrite=True)


