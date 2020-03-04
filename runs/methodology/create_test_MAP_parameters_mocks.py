import os
import argparse
import numpy as np


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
            'a_vir_bin_1': [0.0, 12.0], 
            'b1_bin_2': [0.5, 9.0], 
            'b2_bin_2': [-4.0, 8.0], 
            'gamma3_bin_2': [-8.0, 8.0], 
            'a_vir_bin_2': [0.0, 12.0]
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
        output_dir = "runs/methodology/data/noisy_fiducial/random_start{j}/"
        os.makedirs(f'{output_dir}', exist_ok=True)

        for i in range(args.noise_range[0], args.noise_range[1]):
            starting_point_settings = sample_random_starting_point()
            np.save('{output_dir}start{j}_noise{i}.npy', starting_point_settings)

