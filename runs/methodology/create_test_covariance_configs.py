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

    # No dz
    #dz_cov_file = "data/KV450/nofz/id_cov.asc"

    # Case for 2x2pt
    run_type = "EE_nE"

    KiDS_twopoint_tag_list = ['theoryBuceros', 'theoryEgretta', 'simBuceros', 'simEgretta']
    data_name_root_list    = ['theory_simple', 'theory_complex', 'mock_simple', 'mock_complex']
    cut_dz_modules         = ["--cut-modules", "correlated_dz_priors",
                              "--cut-modules", "source_photoz_bias"]

    # Loop over different covariances
    for KiDS_twopoint_tag, data_name_root in zip(KiDS_twopoint_tag_list, data_name_root_list):

        # Configs for noiseless cases; only made if noise_begin = 0
        if noise_begin == 0:

            root_data_dir = f"runs/methodology/data/noisefree_fiducial/{data_name_root}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, f"KiDS/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad_mock_noiseless.fits")

            # Multinest
            output_root_dir = f"runs/methodology/test_covariance/{data_name_root}/multinest"
            multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                                  "--sampler-config", "nested_sampling_tolerance", "1.0e-2"]
            run_name_root = "multinest"
            run_name = f"{run_name_root}_{run_type}"

            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--sampler", "multinest",
                    *multinest_settings,
                    *cut_dz_modules]
            subprocess.run(["python", script] + cmd, check=True)

            # Noiseless MAP sampler
            output_root_dir = f"runs/methodology/test_covariance/{data_name_root}/MAP_noiseless"
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--sampler", "maxlike",
                    *cut_dz_modules]
            subprocess.run(["python", script] + cmd, check=True)

            # Test sampler
            output_root_dir = f"runs/methodology/test_covariance/{data_name_root}/test"
            run_name_root = "test_sampler"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--sampler", "test",
                    *cut_dz_modules]
            subprocess.run(["python", script] + cmd, check=True)

        # Noisy MAP runs; loop over noise realizations
        output_root_dir = f"runs/methodology/test_covariance/{data_name_root}/MAP"
        run_name_root = "MAP"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"]

        for i in range(noise_begin, noise_end):
            root_data_dir = f"runs/methodology/data/noisy_fiducial/{data_name_root}_{i}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, f"KiDS/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad_mock_noisy.fits")
            run_name = f"{run_name_root}_{i}_{run_type}"

            cmd = ["--root-dir", output_root_dir,
                  "--run-name", run_name,
                  "--run-type", run_type,
                  "--KiDS-data-file", twopoint_file,
                  "--sampler", "maxlike",
                  *MAP_settings,
                  *cut_dz_modules]
            subprocess.run(["python", script] + cmd, check=True)

