import subprocess
import os
import argparse
import sys
sys.path.append("modules/scale_cuts/")
from wrapper_twopoint import TwoPointWrapper


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free mocks')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy mocks indexed by i with BEGIN <= i < END')
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
            "--set-parameters", "nofz_shifts", "p_2", "2.0"
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

