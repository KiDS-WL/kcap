import subprocess
import os
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free configs')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy configs indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # In `main_chains`: SOM dz covariance scaled by factor 4
    test_name   = "main_chains"
    dz_cov_file = "data/KV450/nofz/SOM_cov_multiplied.asc"

    # Run over all possible combinations
    run_type_list = ["EE_nE_w", "EE_nE", "EE_w", "nE_w", "EE", "nE", "w"]

    # Loop over different run_type
    for run_type in run_type_list:

        # Configs for noise-free cases
        if args.noise_free:
            
            root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits")
            boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_1.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_2.txt")]
            boss_cov_files =  [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

            # Multinest
            output_root_dir = "runs/methodology/{test_name}/multinest"
            multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                                  "--sampler-config", "nested_sampling_tolerance", "1.0e-2"]

            run_name_root = "multinest"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", "multinest",
                    *multinest_settings]
            subprocess.run(["python", script] + cmd, check=True)

            # Noiseless MAP sampler
            output_root_dir = "runs/methodology/{test_name}/MAP_noiseless"
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", "maxlike"]
            subprocess.run(["python", script] + cmd, check=True)

            # Test sampler
            output_root_dir = "runs/methodology/{test_name}/test"
            run_name_root = "test_sampler"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", "test"]
            subprocess.run(["python", script] + cmd, check=True)

        # Noisy MAP runs; loop over noise realizations
        output_root_dir = "runs/methodology/{test_name}/MAP"
        run_name_root = "MAP"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"]

        for i in range(args.noise_range[0], args.noise_range[1]):
            root_data_dir   = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
            twopoint_file   = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")
            boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_1.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_2.txt")]
            boss_cov_files  = [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

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

