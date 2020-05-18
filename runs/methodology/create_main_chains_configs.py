import subprocess
import os
import argparse
import numpy as np


def getPriorDict(tag):
    if tag in ['EE_nE_w', 'EE_w']:
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
    elif tag in ['EE_nE']:
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
    elif tag in ['EE']:
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
            }
        }
    elif tag in ['nE_w']:
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

            #'intrinsic_alignment_parameters': {
                #'A': [-6.0, 6.0]
            #},

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
    elif tag in ['nE']:
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

            #'intrinsic_alignment_parameters': {
                #'A': [-6.0, 6.0]
            #},

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
    elif tag in ['w']:
        priorDict = {
            'cosmological_parameters': {
                'omch2': [0.051, 0.2], 
                'ombh2': [0.019, 0.026], 
                'h0': [0.64, 0.82], 
                'n_s': [0.84, 1.1], 
                'S_8_input': [0.1, 1.3]
            }, 

            #'halo_model_parameters': {
                #'A': [2.0, 3.13]
            #},

            #'intrinsic_alignment_parameters': {
                #'A': [-6.0, 6.0]
            #},

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
    parser.add_argument('--MAP-free', action='store_true', help='create MAP noise-free configs')
    parser.add_argument('--noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy configs indexed by i with BEGIN <= i < END')
    parser.add_argument('--m-noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy configs for multinest indexed by i with BEGIN <= i < END')
    parser.add_argument('--ma-noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy configs for Marika indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    # In `main_chains`: SOM dz covariance scaled by factor 4
    test_name   = "main_chains"
    dz_cov_file = "data/KV450/nofz/SOM_cov_multiplied.asc"

    # Run over all possible combinations
    run_type_list = ["EE_nE_w", "EE_nE", "EE_w", "nE_w", "EE", "nE", "w"]

    # Loop over different run_type
    for run_type in run_type_list:
        root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
        twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits")
        boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_1.txt"),
                          os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_2.txt")]
        boss_cov_files =  [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                          os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

        # Configs for noise-free cases
        if args.noise_free:

            # Multinest
            output_root_dir = f"runs/methodology/{test_name}/multinest"
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

            # Test sampler
            output_root_dir = f"runs/methodology/{test_name}/test"
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

        # Noiseless MAP sampler
        if args.MAP_free:
            output_root_dir = f"runs/methodology/{test_name}/MAP_noiseless"
            multi_max_dir = f"runs/methodology/data/multinest_prior/"
            data = np.load(f'{multi_max_dir}multinest_main_{run_type}_samples.npy')
            post = np.load(f'{multi_max_dir}multinest_main_{run_type}_post.npy')
            starting_point_settings = get_multinest_max_as_start(data, post, run_type)
            
            run_name_root = "MAP"
            run_name = f"{run_name_root}_{run_type}" 
            cmd = ["--root-dir", output_root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", "maxlike",
                    *starting_point_settings]
            subprocess.run(["python", script] + cmd, check=True)

    # Noisy MAP runs; loop over noise realizations; EE_nE_w only
    # Loop over different run_type
    for run_type in ['EE_nE_w', 'EE_nE']:
        output_root_dir  = f"runs/methodology/{test_name}/MAP"
        random_start_dir = f"runs/methodology/data/multinest_start/{test_name}/{run_type}/"
        run_name_root = "MAP"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01"]
        
        for i in range(args.noise_range[0], args.noise_range[1]):
            root_data_dir   = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
            twopoint_file   = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")
            boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_1.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_2.txt")]
            boss_cov_files  = [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

            # Multinest start
            random_start_file = f'{random_start_dir}start{i}.npy'
            starting_point_settings = np.load(random_start_file)

            run_name = f"{run_name_root}_{i}_{run_type}"
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

    # Noisy multinest runs; loop over noise realizations; EE_nE_w only
    # Loop over different run_type
    for run_type in ['EE_nE']:
        output_root_dir  = f"runs/methodology/{test_name}/multinest_noisy"
        run_name_root = "multinest"
        multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                              "--sampler-config", "nested_sampling_tolerance", "1.0e-2"]
        
        for i in range(args.m_noise_range[0], args.m_noise_range[1]):
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
                    "--sampler", "multinest",
                    *multinest_settings]
            subprocess.run(["python", script] + cmd, check=True)

    ## Marika's
    for run_type in ['EE_nE']:
        output_root_dir  = f"runs/methodology/{test_name}/MAP_Marika"
        multi_max_dir = f"runs/methodology/data/multinest_prior/"
        run_name_root = "MAP_Marika"
        MAP_settings = ["--sampler-config", "maxlike_tolerance", "0.01",
                        "--sampler-config", "max_posterior", '']
        
        for i in range(args.ma_noise_range[0], args.ma_noise_range[1]):
            root_data_dir   = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
            twopoint_file   = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits")
            boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_1.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_2.txt")]
            boss_cov_files  = [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                              os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

            data = np.load(f'{multi_max_dir}multinest_main_{run_type}_{i}_samples.npy')
            post = np.load(f'{multi_max_dir}multinest_main_{run_type}_{i}_post.npy')
            starting_point_settings = get_multinest_max_as_start(data, post, run_type)

            run_name = f"{run_name_root}_{i}_{run_type}"
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

