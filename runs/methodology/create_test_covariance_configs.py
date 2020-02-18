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

    # Identity
    dz_cov_file = "data/KV450/nofz/id_cov.asc"

    # Case for 3x2pt. For the full setup we'd want to loop over the run_type
    run_type  = "EE_nE"

    cov_tag_list = ['theory_simple', 'theory_complex', 'mock_simple', 'mock_complex']
    
    for cov_tag in cov_tag_list:

        # MAP runs
        output_root_dir = f"runs/methodology/test_covariance/{cov_tag}/MAP"
        run_name_root = "MAP"

        for i in range(noise_begin, noise_end):
            root_data_dir = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")

            run_name = f"{run_name_root}_{i}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                  "--run-name", run_name,
                  "--run-type", run_type,
                  "--KiDS-data-file", twopoint_file,
                  "--dz-covariance-file", dz_cov_file,
                  "--sampler", "maxlike",
                  "--sampler-config", "maxlike_tolerance", "0.01"]
            subprocess.run(["python", script] + cmd, check=True)

        if noise_begin == 0:
            
            # Multinest and test sampler runs
            root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits")

            # Multinest
            output_root_dir = f"runs/methodology/test_covariance/{cov_tag}/multinest"
            multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                                  "--sampler-config", "nested_sampling_tolerance", "1.0e-2"]

            run_name_root = "multinest"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "multinest",
                    *multinest_settings]
            subprocess.run(["python", script] + cmd, check=True)

            # Noiseless MAP sampler
            output_root_dir = f"runs/methodology/test_covariance/{cov_tag}/MAP_noiseless"
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "maxlike"]
            subprocess.run(["python", script] + cmd, check=True)

            # Test sampler
            output_root_dir = f"runs/methodology/test_covariance/{cov_tag}/test"
            run_name_root = "test_sampler"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "test"]
            subprocess.run(["python", script] + cmd, check=True)
