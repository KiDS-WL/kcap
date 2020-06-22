import subprocess
import os
import argparse
import numpy as np


def getPriorDict(tag):
    if tag in ['S8U']:
        priorDict = {
            'cosmological_parameters': {
                'omch2': [0.051, 0.255], 
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
                #'a_vir_bin_1': [0.000001, 12.0], ## 1e-6 instead of 0 due to a bug
                'b1_bin_2': [0.5, 9.0], 
                'b2_bin_2': [-4.0, 8.0], 
                'gamma3_bin_2': [-8.0, 8.0], 
                #'a_vir_bin_2': [0.000001, 12.0] ## 1e-6 instead of 0 due to a bug
            }
        }
    elif tag in ['lnAsC', 'lnAsU']:
        priorDict = {
            'cosmological_parameters': {
                'omch2': [0.051, 0.255], 
                'ombh2': [0.019, 0.026], 
                'h0': [0.64, 0.82], 
                'n_s': [0.84, 1.1], 
                'ln_1e10_a_s': [1.5, 4.0]
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
                #'a_vir_bin_1': [0.000001, 12.0], ## 1e-6 instead of 0 due to a bug
                'b1_bin_2': [0.5, 9.0], 
                'b2_bin_2': [-4.0, 8.0], 
                'gamma3_bin_2': [-8.0, 8.0], 
                #'a_vir_bin_2': [0.000001, 12.0] ## 1e-6 instead of 0 due to a bug
            }
        }
    return priorDict

def get_multinest_max_as_start(data, post, tag):
    max_ind = np.argmax(post)
    best_fit = data[:, max_ind]
    priorDict = getPriorDict(tag)
    
    lower = []
    upper = []
    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            lower.append(value2[0])
            upper.append(value2[1])
    lower = np.array(lower)
    upper = np.array(upper)

    starting_point_settings = []
    ind = 0

    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            p = best_fit[ind]
            starting_point_settings.append("--set-parameters")
            starting_point_settings.append(key1)
            starting_point_settings.append(key2)
            starting_point_settings.append(f'{value2[0]} {p} {value2[1]}')
            ind += 1

    return starting_point_settings
  
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free configs')
    #parser.add_argument('--MAP-free', action='store_true', help='create MAP noise-free configs')
    #parser.add_argument('--noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy configs indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # In `test_priors`: 2x2pt
    test_name = "test_priors"
    run_type  = "EE_nE"
    uncorr_dz_cov_file = "data/KV450/nofz/SOM_cov_multiplied_diagonal.asc" # uncorrelated dz covariance
    corr_dz_cov_file   = "data/KV450/nofz/SOM_cov_multiplied.asc" # SOM dz covariance scaled by factor 4
    sampler_ln_As_settings = [
        "--enable-modules", "sample_ln_As",
        "--set-parameters", "cosmological_parameters", "S_8_input", "0.7458",
        "--cut-modules", "sample_S8",
        "--cut-modules", "sigma8toAs"
    ]
    sampler_ln_As_settings_2 = [
        "--set-parameters", "cosmological_parameters", "ln_1e10_a_s", "1.5 2.72 4.0"
    ]

    data_name_root = "S8_corr"

    # S8 corr is a repeated case; establish symbolic links
    if args.noise_free:
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}', exist_ok=True)
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/multinest', exist_ok=True)
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/test', exist_ok=True)

        destination = f"runs/methodology/{test_name}/{data_name_root}/multinest/multinest_{run_type}"
        if not os.path.islink(destination):
            os.symlink(f"../../../main_chains/multinest/multinest_{run_type}", destination)

        destination = f"runs/methodology/{test_name}/{data_name_root}/test/test_sampler_{run_type}"
        if not os.path.islink(destination):
            os.symlink(f"../../../main_chains/test/test_sampler_{run_type}", destination)

    if args.MAP_free:
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/MAP_noiseless', exist_ok=True)
        destination = f"runs/methodology/{test_name}/{data_name_root}/MAP_noiseless/MAP_{run_type}"
        if not os.path.islink(destination):
            os.symlink(f"../../../main_chains/MAP_noiseless/MAP_{run_type}", destination)

    for i in range(args.noise_range[0], args.noise_range[1]):
        os.makedirs(f'runs/methodology/{test_name}/{data_name_root}/MAP', exist_ok=True)

        destination = f"runs/methodology/{test_name}/{data_name_root}/MAP/MAP_{i}_{run_type}"
        if not os.path.islink(destination):
            os.symlink(f"../../../main_chains/MAP/MAP_{i}_{run_type}", destination)

    data_name_root_list = ["S8_uncorr", "lnAs_corr", "lnAs_uncorr"] # Skip "S8_corr" as it is the same as in a case from main chains
    data_name_tag_list = ['S8U', 'lnAsC', 'lnAsU']

    # Loop over different prior
    for data_name_root, data_name_tag in zip(data_name_root_list, data_name_tag_list):
        tag_list = data_name_root.split('_')

        if tag_list[0] == "lnAs":
            S8_lnAs_settings = sampler_ln_As_settings + sampler_ln_As_settings_2
        else:
            S8_lnAs_settings = []

        if tag_list[1] == "uncorr":
            dz_cov_file = uncorr_dz_cov_file
        else:
            dz_cov_file = corr_dz_cov_file


        root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
        twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits")

        # Configs for noise-free cases
        if args.noise_free:

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

        # Noiseless MAP sampler
        if args.MAP_free:
            output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP_noiseless"
            multi_max_dir = f"runs/methodology/data/multinest_prior/"
            data = np.load(f'{multi_max_dir}multinest_prior_{data_name_tag}_samples.npy')
            post = np.load(f'{multi_max_dir}multinest_prior_{data_name_tag}_post.npy')
            starting_point_settings = get_multinest_max_as_start(data, post, data_name_tag)
            
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "maxlike",
                    *sampler_ln_As_settings,
                    *starting_point_settings]
            subprocess.run(["python", script] + cmd, check=True)

        # Noisy MAP runs; loop over noise realizations
        output_root_dir = f"runs/methodology/{test_name}/{data_name_root}/MAP"
        random_start_dir = f"runs/methodology/data/multinest_start/{test_name}/{data_name_root}/"
        run_name_root = "MAP"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"]

        for i in range(args.noise_range[0], args.noise_range[1]):
            root_data_dir = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
            twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")

            # Multinest start
            random_start_file = f'{random_start_dir}start{i}.npy'
            starting_point_settings = np.load(random_start_file)

            run_name = f"{run_name_root}_{i}_{run_type}"
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--sampler", "maxlike",
                    *MAP_settings,
                    *sampler_ln_As_settings,
                    *starting_point_settings]
            subprocess.run(["python", script] + cmd, check=True)

