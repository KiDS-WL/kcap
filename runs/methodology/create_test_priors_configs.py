import subprocess
import os
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free configs')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy configs indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # In `test_priors`: 2x2pt
    test_name = "test_priors"
    run_type  = "EE_nE"
    uncorr_dz_cov_file = "data/KV450/nofz/id_cov.asc" # uncorrelated dz covariance
    corr_dz_cov_file   = "data/KV450/nofz/SOM_cov_multiplied.asc" # SOM dz covariance scaled by factor 4
    sampler_ln_As_settings = [
        "--enable-modules", "sample_ln_As",
        "--cut-modules", "sample_S8",
        "--cut-modules", "sigma8toAs"
    ]

    data_name_root = "S8_corr"

    # S8 corr is a repeated case; establish symbolic links
    if args.noise_free:
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}', exist_ok=True)
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/multinest', exist_ok=True)
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/MAP_noiseless', exist_ok=True)
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/test', exist_ok=True)
        os.symlink(f"../../../main_chains/multinest/multinest_{run_type}", f"runs/methodology/{test_name}/{data_name_root}/multinest/multinest_{run_type}")
        os.symlink(f"../../../main_chains/MAP_noiseless/MAP_{run_type}", f"runs/methodology/{test_name}/{data_name_root}/MAP_noiseless/MAP_{run_type}")
        os.symlink(f"../../../main_chains/test/test_sampler_{run_type}", f"runs/methodology/{test_name}/{data_name_root}/test/test_sampler_{run_type}")

    for i in range(args.noise_range[0], args.noise_range[1]):
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/MAP', exist_ok=True)
        os.symlink(f"../../../main_chains/MAP/MAP_{i}_{run_type}", f"runs/methodology/{test_name}/{data_name_root}/MAP/MAP_{i}_{run_type}")

    data_name_root_list = ["S8_uncorr", "lnAs_corr", "lnAs_uncorr"] # Skip "S8_corr" as it is the same as in a case from main chains

    # Loop over different prior
    for data_name_root in data_name_root_list:
        tag_list = data_name_root.split('_')

        if tag_list[0] == "lnAs":
            S8_lnAs_settings = sampler_ln_As_settings
        else:
            S8_lnAs_settings = []

        if tag_list[1] == "uncorr":
            dz_cov_file = uncorr_dz_cov_file
        else:
            dz_cov_file = corr_dz_cov_file

        # Configs for noise-free cases
        if args.noise_free:

            root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits")

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
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "multinest",
                    *multinest_settings,
                    *S8_lnAs_settings]
            subprocess.run(["python", script] + cmd, check=True)

            # Noiseless MAP sampler
            output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP_noiseless"
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "maxlike",
                    *S8_lnAs_settings]
            subprocess.run(["python", script] + cmd, check=True)

            # Test sampler
            output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/test"
            run_name_root = "test_sampler"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "test",
                    *S8_lnAs_settings]
            subprocess.run(["python", script] + cmd, check=True)

        # Noisy MAP runs; loop over noise realizations
        output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP"
        run_name_root = "MAP"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"]

        for i in range(args.noise_range[0], args.noise_range[1]):
            root_data_dir = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")
            run_name = f"{run_name_root}_{i}_{run_type}"

            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "maxlike",
                    *MAP_settings,
                    *S8_lnAs_settings]
            subprocess.run(["python", script] + cmd, check=True)

