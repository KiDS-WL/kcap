import subprocess
import os
import argparse
import sys
sys.path.append("modules/scale_cuts/")
from wrapper_twopoint import TwoPointWrapper
import numpy as np
import scipy.stats as stats


def sample_random_from_multinest_EE_nE_S8(data, wgt):
    kernel = stats.gaussian_kde(data, 'silverman', weights=wgt)
    random_start = kernel.resample(1).flatten()
    
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
            #'a_vir_bin_1': [0.000001, 12.0], ## 1e-6 instead of 0 due to a bug
            'b1_bin_2': [0.5, 9.0], 
            'b2_bin_2': [-4.0, 8.0], 
            'gamma3_bin_2': [-8.0, 8.0], 
            #'a_vir_bin_2': [0.000001, 12.0] ## 1e-6 instead of 0 due to a bug
        }
    }

    lower = []
    upper = []
    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            lower.append(value2[0])
            upper.append(value2[1])
    lower = np.array(lower)
    upper = np.array(upper)

    while np.any(random_start <= lower) or np.any(upper <= random_start):
        random_start = kernel.resample(1).flatten()

    starting_point_settings = []
    ind = 0

    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            p = random_start[ind]
            starting_point_settings.append("--set-parameters")
            starting_point_settings.append(key1)
            starting_point_settings.append(key2)
            starting_point_settings.append(f'{value2[0]} {p} {value2[1]}')
            ind += 1

    return starting_point_settings
    
def sample_random_from_multinest_EE_nE_lnAs(data, wgt):
    kernel = stats.gaussian_kde(data, 'silverman', weights=wgt)
    random_start = kernel.resample(1).flatten()
    
    priorDict = {
        'cosmological_parameters': {
            'omch2': [0.051, 0.2], 
            'ombh2': [0.019, 0.026], 
            'h0': [0.64, 0.82], 
            'n_s': [0.84, 1.1], 
            #'S_8_input': [0.1, 1.3],
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

    lower = []
    upper = []
    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            lower.append(value2[0])
            upper.append(value2[1])
    lower = np.array(lower)
    upper = np.array(upper)

    while np.any(random_start <= lower) or np.any(upper <= random_start):
        random_start = kernel.resample(1).flatten()

    starting_point_settings = []
    ind = 0

    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            p = random_start[ind]
            starting_point_settings.append("--set-parameters")
            starting_point_settings.append(key1)
            starting_point_settings.append(key2)
            starting_point_settings.append(f'{value2[0]} {p} {value2[1]}')
            ind += 1

    return starting_point_settings
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--multi-start', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create random staring point from multinest prior indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    if args.noise_range[1] > args.noise_range[0]:
        print('You want to create noisy mocks.')
        print('Note that noisy mocks for test_covariance require those for main_chains.')
        print('Make sure you have already created them.')

    data_name_root_list = ["S8_uncorr", "lnAs_corr", "lnAs_uncorr"]
    key_list            = ['S8U', 'lnAsC', 'lnAsU']

    # Multinest random for test of goodness of fit
    for data_name_root, key in zip(data_name_root_list, key_list):
        if args.multi_start[0] < args.multi_start[1]:
            multinest_prior_name = f'runs/methodology/data/multinest_prior/multinest_prior_{key}_samples.npy'
            weight_name = f'runs/methodology/data/multinest_prior/multinest_prior_{key}_weight.npy'
            data = np.load(multinest_prior_name)
            wgt  = np.load(weight_name)
            output_dir = f"runs/methodology/data/multinest_start/test_priors/{data_name_root}/"
            os.makedirs(f'{output_dir}', exist_ok=True)

            for i in range(args.multi_start[0], args.multi_start[1]):
                starting_point_settings = sample_random_from_multinest_EE_nE(data, wgt)
                np.save(f'{output_dir}start{i}.npy', starting_point_settings)

