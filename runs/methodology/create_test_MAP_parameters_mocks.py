import os
import argparse
import numpy as np
import scipy.stats as stats


def sample_random_starting_point():
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

    starting_point_settings = []

    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            p = np.random.uniform(value2[0], value2[1], 1)[0]
            starting_point_settings.append("--set-parameters")
            starting_point_settings.append(key1)
            starting_point_settings.append(key2)
            starting_point_settings.append(f'{value2[0]} {p} {value2[1]}')

    return starting_point_settings

def sample_random_from_multinest(data, wgt):
    kernel = stats.gaussian_kde(data, 'silverman', weights=wgt)
    random_start = kernel.resample(1)
    
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

    lower = []
    upper = []
    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            lower.append(value2[0])
            upper.append(value2[1])
    lower = np.array(lower)
    upper = np.array(upper)

    while np.any(random_start <= lower) or np.any(upper <= random_start):
        random_start = kernel.resample(1)

    starting_point_settings = []
    ind = 0

    for key1, value1 in priorDict.items():
        for key2, value2 in value1.items():
            p = random_start[ind][0]
            starting_point_settings.append("--set-parameters")
            starting_point_settings.append(key1)
            starting_point_settings.append(key2)
            starting_point_settings.append(f'{value2[0]} {p} {value2[1]}')
            ind += 1

    return starting_point_settings


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    #parser.add_argument('--noise-free', action='store_true', help='create noise-free mocks')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy mocks indexed by i with BEGIN <= i < END')
    parser.add_argument('--random-start-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create configs with random starting points indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"

    for j in range(args.random_start_range[0], args.random_start_range[1]):

        # Do nothing for noise-free cases

        # Create noisy mocks
        output_dir = f"runs/methodology/data/noisy_fiducial/random_start{j}/"
        os.makedirs(f'{output_dir}', exist_ok=True)

        for i in range(args.noise_range[0], args.noise_range[1]):
            starting_point_settings = sample_random_starting_point()
            np.save(f'{output_dir}start{j}_noise{i}.npy', starting_point_settings)

        # Multinest
        multinest_prior_name = 'runs/methodology/data/multinest_prior/multinest_samples.npy'
        weight_name = 'runs/methodology/data/multinest_prior/multinest_random_weight.npy'
        data = np.load(multinest_prior_name)
        wgt  = np.load(weight_name)
        output_dir = f"runs/methodology/data/noisy_fiducial/multinest_start{j}/"
        os.makedirs(f'{output_dir}', exist_ok=True)

        for i in range(args.noise_range[0], args.noise_range[1]):
            starting_point_settings = sample_random_from_multinest(data, wgt)
            np.save(f'{output_dir}start{j}_noise{i}.npy', starting_point_settings)





