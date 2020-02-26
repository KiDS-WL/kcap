import subprocess
import os
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free configs')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy configs indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # Let dz_cov fall back to the default value; we are fixing dz shift anyway
    #dz_cov_file = "data/KV450/nofz/id_cov.asc"

    # In`test_covariance`: 2x2pt, fixed dz to 0
    test_name = "test_covariance"
    run_type  = "EE_nE"
    fix_dz    = ["--fix-values", "nofz_shifts"]

    KiDS_twopoint_tag_list = ['theoryBuceros', 'theoryEgretta', 'simBuceros', 'simEgretta']
    data_name_root_list    = ['theory_simple', 'theory_complex', 'mock_simple', 'mock_complex']

    # Loop over different covariances
    for KiDS_twopoint_tag, data_name_root in zip(KiDS_twopoint_tag_list, data_name_root_list):

        # Configs for noise-free cases
        if args.noise_free:

            root_data_dir = f"runs/methodology/data/noisefree_fiducial/{data_name_root}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, f"KiDS/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad_mock_noiseless.fits")

            # Multinest
            output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/multinest"
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
                    *fix_dz]
            subprocess.run(["python", script] + cmd, check=True)

            # Noiseless MAP sampler
            output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP_noiseless"
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--sampler", "maxlike",
                    *fix_dz]
            subprocess.run(["python", script] + cmd, check=True)

            # Test sampler
            output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/test"
            run_name_root = "test_sampler"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--sampler", "test",
                    *fix_dz]
            subprocess.run(["python", script] + cmd, check=True)

        # Noisy MAP runs; loop over noise realizations
        output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP"
        run_name_root = "MAP"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"]

        for i in range(args.noise_range[0], args.noise_range[1]):
            root_data_dir = f"runs/methodology/data/noisy_fiducial/{data_name_root}_{i}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, f"KiDS/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad_mock_noisy.fits")
            run_name = f"{run_name_root}_{i}_{run_type}"

            cmd = ["--root-dir", output_root_dir,
                  "--run-name", run_name,
                  "--run-type", run_type,
                  "--KiDS-data-file", twopoint_file,
                  "--sampler", "maxlike",
                  *MAP_settings,
                  *fix_dz]
            subprocess.run(["python", script] + cmd, check=True)

