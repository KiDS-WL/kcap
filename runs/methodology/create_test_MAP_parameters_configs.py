import subprocess
import os
import argparse
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    #parser.add_argument('--noise-free', action='store_true', help='create noise-free configs')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy configs indexed by i with BEGIN <= i < END')
    parser.add_argument('--do-fiducial', action='store_true', help='Change tolerance values without randomizing starting points')
    parser.add_argument('--random-start-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create configs with random starting points indexed by i with BEGIN <= i < END')
    parser.add_argument('--do-multinest', action='store_true', help='Create a specific multinest chain')
    
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # In `test_MAP_parameters`: 3x2pt, SOM dz covariance scaled by factor 4
    test_name   = "test_MAP_parameters"
    run_type    = "EE_nE_w"
    dz_cov_file = "data/KV450/nofz/SOM_cov_multiplied.asc"

    data_name_root = "base"

    # base is EE_nE_w in main chains; create symbolic links
    for i in range(args.noise_range[0], args.noise_range[1]):
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/MAP', exist_ok=True)
        
        destination = f"runs/methodology/{test_name}/{data_name_root}/MAP/MAP_{i}_{run_type}"
        if not os.path.islink(destination):
            os.symlink(f"../../../main_chains/MAP/MAP_{i}_{run_type}", destination)

    data_name_root_list = [
        "tolerance0.005",
        "tolerance0.003",
        "tolerance0.001"
    ]
    MAP_settings_list = [
        ["--sampler-config", "maxlike_tolerance", "0.005"],
        ["--sampler-config", "maxlike_tolerance", "0.003"],
        ["--sampler-config", "maxlike_tolerance", "0.001"]
    ]
    
    # Do nothing for noise-free cases

    # Noisy MAP runs; loop over noise realizations
    for i in range(args.noise_range[0], args.noise_range[1]):
        root_data_dir = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
        twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")
        boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_1.txt"),
                          os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_2.txt")]
        boss_cov_files  = [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                          os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

        # Original start, different MAP parameters
        if args.do_fiducial:
            for data_name_root, MAP_settings in zip(data_name_root_list, MAP_settings_list):
                output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP"
                run_name_root = "MAP"

                run_name = f"{run_name_root}_{i}_{run_type}"
                cmd = ["--root-dir", output_root_dir,
                        "--run-name", run_name,
                        "--run-type", run_type,
                        "--KiDS-data-file", twopoint_file,
                        "--dz-covariance-file", dz_cov_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--BOSS-covariance-files", *boss_cov_files,
                        "--sampler", "maxlike",
                        *MAP_settings]
                subprocess.run(["python", script] + cmd, check=True)

        for j in range(args.random_start_range[0], args.random_start_range[1]):
            output_dir = f"runs/methodology/data/noisy_fiducial/random_start{j}/"
            random_start_file = f'{output_dir}start{j}_noise{i}.npy'
            starting_point_settings = np.load(random_start_file)
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{i}_{run_type}"

            output_root_dir = f"runs/methodology/{test_name}/base_start{j}/MAP"
            MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"] ## Original settings from the main chains

            # Random start, original MAP parameters
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", "maxlike",
                    *MAP_settings,
                    *starting_point_settings]
            subprocess.run(["python", script] + cmd, check=True)

            # Random start, different MAP parameters
            for data_name_root, MAP_settings in zip(data_name_root_list, MAP_settings_list):
                output_root_dir = f"runs/methodology/{test_name}/{data_name_root}_start{j}/MAP"

                cmd = ["--root-dir", output_root_dir,
                        "--run-name", run_name,
                        "--run-type", run_type,
                        "--KiDS-data-file", twopoint_file,
                        "--dz-covariance-file", dz_cov_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--BOSS-covariance-files", *boss_cov_files,
                        "--sampler", "maxlike",
                        *MAP_settings,
                        *starting_point_settings]
                subprocess.run(["python", script] + cmd, check=True)

            # Multinest start
            output_dir = f"runs/methodology/data/noisy_fiducial/multinest_start{j}/"
            random_start_file = f'{output_dir}start{j}_noise{i}.npy'
            starting_point_settings = np.load(random_start_file)
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{i}_{run_type}"
            
            output_root_dir = f"runs/methodology/{test_name}/base_multiStart{j}/MAP"
            MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"] ## Original settings from the main chains

            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", "maxlike",
                    *MAP_settings,
                    *starting_point_settings]
            subprocess.run(["python", script] + cmd, check=True)


    if args.do_multinest:
        root_data_dir = f"runs/methodology/data/noisy_fiducial/base_38_EE_nE_w/data/"
        twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")
        boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_1.txt"),
                          os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_2.txt")]
        boss_cov_files  = [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                          os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

        output_dir = f"runs/methodology/data/noisy_fiducial/random_start2/"
        random_start_file = f'{output_dir}start2_noise38.npy'
        starting_point_settings = np.load(random_start_file)
        run_name_root = "multinest"
        run_name = f"{run_name_root}_38_{run_type}"

        output_root_dir = f"runs/methodology/{test_name}/base_start2/multinest"
        multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                              "--sampler-config", "nested_sampling_tolerance", "1.0e-2"]

        cmd = ["--root-dir", output_root_dir,
                "--run-name", run_name,
                "--run-type", run_type,
                "--KiDS-data-file", twopoint_file,
                "--dz-covariance-file", dz_cov_file,
                "--BOSS-data-files", *boss_data_files,
                "--BOSS-covariance-files", *boss_cov_files,
                "--sampler", "multinest",
                *multinest_settings,
                *starting_point_settings]
        subprocess.run(["python", script] + cmd, check=True)
