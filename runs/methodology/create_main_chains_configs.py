import sys
import subprocess
import os

if __name__ == "__main__":

    if len(sys.argv) > 1:
        noise_begin = int(sys.argv[1])
        noise_end   = int(sys.argv[2])
    else:
        noise_begin = 0
        noise_end   = 1

    script = "utils/run_kcap.py"

    # Updated to the SOM covariance
    dz_cov_file = "data/KV450/nofz/SOM_cov_multiplied.asc"

    # To be expanded
    run_type_list = ["EE_nE_w"]

    # Loop over different run_type
    for run_type in run_type_list:

        # Configs for noiseless cases; only made if noise_begin = 0
        if noise_begin == 0:
            
            root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits")
            boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_1.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_2.txt")]
            boss_cov_files =  [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

            # Multinest
            output_root_dir = "runs/methodology/main_chains/multinest"
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
            output_root_dir = "runs/methodology/main_chains/MAP_noiseless"
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
            output_root_dir = "runs/methodology/main_chains/test"
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
        output_root_dir = "runs/methodology/main_chains/MAP"
        run_name_root = "MAP"

        for i in range(noise_begin, noise_end):
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
                  "--sampler-config", "maxlike_tolerance", "0.01"]
            subprocess.run(["python", script] + cmd, check=True)

