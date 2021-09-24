import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
	covariance_file = options[option_section, "covariance"]
	section = options[option_section, "section"]
	n_tomo = options[option_section, "n_tomo"]
	n_comp = options[option_section, "n_comp"]

	cov = np.loadtxt(covariance_file)
	L = np.linalg.cholesky(cov) 
	return section, n_tomo, n_comp, L

def execute(block, config):
	section, n_tomo, n_comp, L = config

	names = []
	for nt in range(n_tomo):
		for nc in range(n_comp):
			names.append('a_{0}_{1}'.format(str(nt+1).zfill(2), str(nc+1).zfill(2)))

	p = []
	for name in names:
		p.append(block[section, name])
	p = L @ np.array(p)
	p[p < 0.] = 0.
	p[p > 1.] = 1.
	for nt in range(n_tomo):
		p[nt*n_comp:(nt+1)*n_comp] /= p[nt*n_comp:(nt+1)*n_comp].sum()

	for i, name in enumerate(names):
		block[section, name] = p[i]

	return 0

def clean(config):
	pass
