from __future__ import print_function
from cosmosis.datablock import names, option_section
import numpy as np
from scipy import special
import scipy

def combonent(z, amp, z_mean, sigma):
	norm = np.sqrt(np.pi/2) * z_mean * sigma * special.erfc(-z_mean/(np.sqrt(2)*sigma)) + sigma**2 * np.exp(-z_mean**2/(2*sigma**2))
	return(np.exp(amp) * z * np.exp(-(z-z_mean)**2/(2*sigma**2)) / norm)
def comb_nz(z, amplitudes, means, sigma):
	return(np.sum([combonent(z, amplitudes[i], means[i], sigma) for i in range(len(amplitudes))], axis = 0))

# alternate comb functions with integration over histogram bins
def combonent1(z, z_mean, sigma):	# this version takes no amplitude argument
	norm = np.sqrt(np.pi/2)*z_mean*sigma*special.erfc(-z_mean/(np.sqrt(2)*sigma))+sigma**2*np.exp(-z_mean**2/(2*sigma**2))
	return z * np.exp(-(z-z_mean)**2 / (2*sigma**2)) / norm
# z_bins_mid: central z values of z-bins in the data histgram
# for each combonent: construct a binned version by integrating between data histogram z-bin boundaries
def integrate_combonents(n_comp, z_bins_mid, means, sigma):
	width = z_bins_mid[1] - z_bins_mid[0]
	n_combonents = np.zeros((n_comp, len(z_bins_mid)))
	for i in range(n_combonents.shape[0]):
		for j in range(n_combonents.shape[1]):
			n_combonents[i,j] = scipy.integrate.quad(combonent1, z_bins_mid[j]-width/2, z_bins_mid[j]+width/2, args=(means[i], sigma))[0] / width
	return n_combonents
# This function sums up combonents for a given set of amplitudes and return the n(z) histogram
def combmodel(z, n_tomo, n_combonents, *ln_amps):
	amps = np.exp(np.split(np.array(ln_amps),n_tomo))
	n_z = np.concatenate([np.sum(n_combonents*amps[i][:,None], axis=0) for i in range(n_tomo)])
	return n_z

def setup(options):
	section = options.get_string(option_section, "section")
	z_min = options.get_double(option_section, "z_min")
	z_max = options.get_double(option_section, "z_max")
	sigma = options.get_double(option_section, "sigma")
	n_comp = options.get_int(option_section, "n_components")
	n_tomo = options.get_int(option_section, "n_tomo")
	nz_hist = options.get_int(option_section, "nz_hist")
	sum_amplitudes = options.get_bool(option_section, "sum_amplitudes")
	means = z_min + np.arange(n_comp)*sigma
	z = np.linspace(z_min, z_max, nz_hist)
	_n_combonents = integrate_combonents(n_comp, z, means, sigma)
	try:
		covariance_file = options.get_string(option_section, "covariance_file", None)
		covariance = np.loadtxt(covariance_file)
	except:
		covariance = 0.

	return {'section':section, 'z':z, 'means':means, 'sigma':sigma, 'covariance':covariance,
			'n_comp':n_comp, 'n_tomo':n_tomo, 'sum_amplitudes':sum_amplitudes, '_n_combonents':_n_combonents}


def execute(block, config):
	section, z, means, sigma, covariance, n_comp, n_tomo, sum_amplitudes, _n_combonents = \
		map(lambda x: config[x], ['section', 'z', 'means', 'sigma', 'covariance', 'n_comp', 'n_tomo', 'sum_amplitudes', '_n_combonents'])

	# read amplitudes from block
	amplitudes = np.zeros((n_tomo, n_comp))
	for nt in range(n_tomo):
		for nc in range(n_comp):
			amplitudes[nt, nc] = block["amplitudes", "a_{0}_{1}".format(str(nt+1).zfill(2), str(nc+1).zfill(2))]
		# per bin amplitudes must sum to unity
		amplitudes[nt] /= amplitudes[nt].sum()
	# take natural log for computations
	amplitudes = np.log(amplitudes)

	data = {}
	if sum_amplitudes:
		print("Sum amplitudes to generate redshift distributions of {0} tomographic bins.".format(n_tomo))
		#nz = np.zeros((n_tomo, len(z)))
		#for i in range(n_tomo):
		#	nz[i,:] = comb_nz(z, amplitudes[i], means, sigma)
		nz = combmodel(z, n_tomo, _n_combonents, *amplitudes)
		nz = nz.reshape(n_tomo, len(z))
	else:
		print("Generate redshift distributions of {0} comb components.".format(n_comp))
		nz = np.zeros((n_comp, len(z)))
		for i in range(n_comp):
			#nz[i,:] = combonent(z, 0, means[i], sigma)
			nz[i,:] = combonent1(z, means[i], sigma)
	data[section] = (z, nz)

	for name, data in list(data.items()):
		z, nz = data
		nbin = len(nz)
		ns = len(z)
		block[name, "nbin"] = nbin
		block[name, "nz"] = ns
		block[name, "z"] = z
		for i, n in enumerate(nz):
			block[name, "bin_{0}".format(i + 1)] = n
	block["amplitudes", "n_tomo"] = n_tomo
	block["amplitudes", "n_comp"] = n_comp
	block["amplitudes", "amp"] = amplitudes
	block["amplitudes", "cov"] = covariance
	block["comb", "means"] = means
	block["comb", "sigma"] = sigma
	return 0

