import os
import argparse
import numpy as np
import scipy.stats as stats
import pandas as pd


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

def loadAsciiChain(name, changeCol={}, lnAs=True, verbose=True):
    f = open(name, 'r')
    hdr = f.readline().split()
    hdr[0] = hdr[0][1:]
    f.close()
    
    data = pd.read_csv(name, dtype=float, sep='\t', comment='#', header=None)
    data.columns = hdr
    
    if verbose:
        print('Loaded \"%s\"' % name)
        print('shape = %s' % str(data.shape))
    
    data.rename(columns=changeCol, inplace=True)
    
    if lnAs:
        col_As   = 'COSMOLOGICAL_PARAMETERS--A_S'
        col_lnAs = 'COSMOLOGICAL_PARAMETERS--ln_1e10_A_S'
        data[col_As] = np.log(data[col_As]*1e+10)
        data.rename(columns={col_As: col_lnAs}, inplace=True)
    return data

def sample_random_from_multinest(data):
    colList = [
        'cosmological_parameters--omch2', 
        'cosmological_parameters--ombh2',
        'cosmological_parameters--h0', 
        'cosmological_parameters--n_s',
        'COSMOLOGICAL_PARAMETERS--S_8',
        
        'halo_model_parameters--a',
        
        'intrinsic_alignment_parameters--a', 
        
        'nofz_shifts--p_1',
        'nofz_shifts--p_2', 
        'nofz_shifts--p_3', 
        'nofz_shifts--p_4',
        'nofz_shifts--p_5', 
        
        'bias_parameters--b1_bin_1',
        'bias_parameters--b2_bin_1', 
        'bias_parameters--gamma3_bin_1',
        'bias_parameters--a_vir_bin_1', 
        'bias_parameters--b1_bin_2',
        'bias_parameters--b2_bin_2', 
        'bias_parameters--gamma3_bin_2',
        'bias_parameters--a_vir_bin_2', 
    ]
    wgt = data['weight'].values
    data2 = [data[col].values for col in colList]
    kernel = stats.gaussian_kde(data2, 'silverman', weights=wgt)
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
        name = 'runs/methodology/main_chains/multinest/multinest_EE_nE_w/chain/samples_EE_nE_w.txt'
        data = loadAsciiChain(name, changeCol={}, lnAs=False, verbose=True)
        output_dir = f"runs/methodology/data/noisy_fiducial/multinest_start{j}/"
        os.makedirs(f'{output_dir}', exist_ok=True)

        for i in range(args.noise_range[0], args.noise_range[1]):
            starting_point_settings = sample_random_from_multinest(data)
            np.save(f'{output_dir}start{j}_noise{i}.npy', starting_point_settings)





