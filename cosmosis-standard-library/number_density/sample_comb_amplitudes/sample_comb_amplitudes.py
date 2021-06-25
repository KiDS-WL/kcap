from builtins import range
from cosmosis.datablock import option_section, names
from scipy.interpolate import interp1d
import numpy as np
from scipy import special

def combonent(z, amp, z_mean, sigma):
    norm = np.sqrt(np.pi/2) * z_mean * sigma * special.erfc(-z_mean/(np.sqrt(2)*sigma)) + sigma**2 * np.exp(-z_mean**2/(2*sigma**2))
    return(np.exp(amp) * z * np.exp(-(z-z_mean)**2/(2*sigma**2)) / norm)
def comb_nz(z, amplitudes, means, sigma):
    return(np.sum([combonent(z, amplitudes[i], means[i], sigma) for i in range(len(amplitudes))], axis = 0))

def setup(options):
    idx_low = options[option_section, "idx_low"]
    idx_high = options[option_section, "idx_high"]
    amplitudes_section = options[option_section, "amplitudes_section"]
    return {"idx_low": idx_low, "idx_high": idx_high, "amplitudes_section": amplitudes_section}


def execute(block, config):
    idx_low = config['idx_low'] - 1
    idx_high = config['idx_high']
    nbin = block["NZ_SOURCE", "nbin"]
    z = block["NZ_SOURCE", "z"]
    amps = config['amplitudes_section']
    amplitudes = block["amplitudes", "amp"]
    means = block["comb", "means"]
    sigma = block["comb", "sigma"]
    n_tomo = amplitudes.shape[0]
    n_comp = amplitudes.shape[1]
    n = block["NZ_SOURCE", "nz"]

    for i in range(idx_high-idx_low):
        amp_name = "amp_%d" % (i+1)
        amp = block[amps, amp_name]
        z_bin = (idx_low+i) // n_comp
        idx_amp = (idx_low+i) % n_comp
        amplitudes[z_bin,idx_amp] = amp
        
    block["amplitudes", "amp"] = amplitudes
    b = block["NZ_SOURCE", "bin_2"]
    for i in range(n_tomo):
        nz = comb_nz(z, amplitudes[i], means, sigma)
        block["NZ_SOURCE", "bin_{0}".format(i + 1)] = nz
    a = block["NZ_SOURCE", "bin_2"]

    return 0


def cleanup(config):
    pass
