import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors
from scipy import linalg

def one_dim_index(bin1, bin2, n_components):
    if bin1 <= bin2:
        return int(bin2 + n_components * bin1 - (bin1 * (bin1 + 1)) / 2)
    else:
        return int(bin1 + n_components * bin2 - (bin2 * (bin2 + 1)) / 2)
def one_dim_index_L_vector(z_bin, comb_bin, n_components):
    if comb_bin>=n_components:
        raise Exception("Bin 2 > n_components!")
    return int(comb_bin + n_components * z_bin)

def setup(options):
    like_name = options.get_string(option_section, "like_name")
    input_section_name = options.get_string(option_section, "input_section_name", default="likelihood")
    analytic_marginalisation = options.get_bool(option_section, "analytic_marginalisation")

    return like_name, input_section_name, analytic_marginalisation

def execute(block, config):
    like_name, input_section_name, analytic_marginalisation = config
    d = block[input_section_name, "data"]
    # Theory vector for all combinations of comb components
    mu_combonents = block[input_section_name, "theory"]
    block[input_section_name, "theory_combonents"] = mu_combonents
    # Set up empty theory vector for correlations between tomographic bins
    mu = np.zeros(d.shape)
    mu_ind = np.zeros((*d.shape, 3))
    # Read comb amplitudes and determine the number of comb components
    amplitudes = np.exp(block["amplitudes", "amp"])
    n_tomo = amplitudes.shape[0]
    n_comp = amplitudes.shape[1]
    for i in range(n_tomo):
        amplitudes[i] /= np.sum(amplitudes[i])
    # Split comb theory vector
    mu_combonents = np.split(np.array(mu_combonents), (n_comp**2+n_comp)/2)
    # Determine number of data points (should be equal to len(mu_combonents[0]) !)
    n_data = int(mu.shape[0]/((n_tomo**2+n_tomo)/2))
    # Fill the theory vector of tomographic bins. Nested loops here; there's probably a smarter method!
#   A_a, A_b, c_ij = [], [], []
    for bin1 in range(n_tomo):
        for bin2 in range(bin1, n_tomo):
            index = one_dim_index(bin1, bin2, n_tomo)
            mu_ind[index*n_data:(index+1)*n_data] = np.column_stack(([bin1]*n_data, [bin2]*n_data, range(n_data)))
            # Both sums sum over all gaussians!
            for comp1 in range(n_comp):
                for comp2 in range(n_comp):
                    index_comp = one_dim_index(comp1, comp2, n_comp)
                    mu[index*n_data:(index+1)*n_data] += amplitudes[bin1, comp1] * amplitudes[bin2, comp2] * mu_combonents[index_comp]
#                    if bin1==0 and bin2==1:
#                        A_a.append(amplitudes[bin1, comp1])
#                        A_b.append(amplitudes[bin2, comp2])
#                        c_ij.append(mu_combonents[index_comp])
#    import h5py
#    A_a, A_b, c_ij = (np.array(x).squeeze() for x in (A_a, A_b, c_ij))
#    with h5py.File('bs_ch_test.h5', 'w') as file:
#        file.create_dataset('A_a', data=A_a)
#        file.create_dataset('A_b', data=A_b)
#        file.create_dataset('c_ij', data=c_ij)
    # Compute fiducial chi-square
    inv_cov = block[input_section_name, "inv_covariance"]
    cov = block[input_section_name, "covariance"]
    bin2, bin1, ang = mu_ind.T + 1
    sorter = np.lexsort((ang, bin2, bin1))
    _mu = mu.copy()
    mu = mu[sorter]
    block[input_section_name, "theory"] = mu
    r = d - mu
    chi2_fid = float(r @ inv_cov @ r)
    if analytic_marginalisation:
        # Load the covariance matrix of comb amplitudes
        covariance = block["amplitudes", "cov"]
        # Again: nested loops; can probably be done more efficiently
        delta_prime = np.zeros((n_tomo * n_comp, int((n_tomo**2+n_tomo)/2) * n_data))
        # Construct delta_prime vector (eq. A.2)
        # t = time.time()
        for mu in range(n_tomo):
            for m in range(n_comp):
                # Determine the vector index for a given z-bin and comb-component
                index = one_dim_index_L_vector(mu, m, n_comp)
                # Loop through redshift bin combinations
                for alpha in range(n_tomo):
                    for beta in range(alpha,n_tomo):
                        # Determine index of redshift bin combination
                        index_zbins = one_dim_index(alpha, beta, n_tomo)
                        x = np.zeros(n_data, 'float64')
                        for i in range(n_comp):
                            factor = 0.
                            # Delta functions (eq. A.2)
                            if alpha == mu:
                                factor += amplitudes[beta,i]
                            if beta == mu:
                                factor += amplitudes[alpha,i]
                            sum_index = one_dim_index(i, m, n_comp)
                            x += mu_combonents[sum_index]*factor
                        # Put everything together
                        delta_prime[index,index_zbins*n_data:(index_zbins+1)*n_data] = (-amplitudes[mu,m]*x)
        # Construct delta_2prime vector (eq. A.4)
        delta_2prime = np.zeros((n_tomo * n_comp, n_tomo * n_comp, int((n_tomo**2+n_tomo)/2) * n_data))
        for mu in range(n_tomo):
            for m in range(n_comp):
                for nu in range(n_tomo):
                    for n in range(n_comp):
                        # Determine the vector index for a given z-bin and comb-component
                        index1 = one_dim_index_L_vector(mu, m, n_comp)
                        index2 = one_dim_index_L_vector(nu, n, n_comp)
                        # Loop through redshift bin combinations
                        for alpha in range(n_tomo):
                            for beta in range(alpha, n_tomo):
                                # Determine index of redshift bin combination
                                index_zbins = one_dim_index(alpha, beta, n_tomo)
                                temp = np.zeros(n_data, 'float64')
                                factor = 0
                                if ((alpha == mu) and (beta == nu)):
                                    factor +=1
                                if ((beta == mu) and (alpha == nu)):
                                    factor +=1
                                if factor > 0:
                                    index_components = one_dim_index(m, n, n_comp)
                                    temp -= (amplitudes[mu,m] * amplitudes[nu,n] * mu_combonents[index_components] * factor)
                                if ((m==n) and (mu==nu)):
                                    s = np.zeros(n_data, 'float64')
                                    for i in range(n_comp):
                                        factor2 = 0
                                        if alpha == mu:
                                            factor2 += (amplitudes[beta,i])
                                        if beta == mu:
                                            factor2 += (amplitudes[alpha,i])
                                        if factor2 > 0:
                                            sum_index = one_dim_index(i, m, n_comp)
                                            s += mu_combonents[index_components] * factor2
                                    temp -= (amplitudes[mu,m] * s)
                                # Put everything together
                                delta_2prime[index1,index2,index_zbins*n_data:(index_zbins+1)*n_data] = temp
        # Compute L' and L'' (eqs. A.1 & A.3)
        L_prime = np.zeros((n_tomo*n_comp))
        L_2prime = np.zeros((n_tomo*n_comp,n_tomo*n_comp))
        for i in range(n_tomo*n_comp):
            L_prime[i] = (delta_prime[i] @ inv_cov @ r) + (r @ inv_cov @ delta_prime[i])
            for j in range(n_tomo*n_comp):
                L_2prime[i,j] = (delta_prime[i] @ inv_cov @ delta_prime[j]) + (delta_prime[j] @ inv_cov @ delta_prime[i]) + (delta_2prime[i,j] @ inv_cov @ r) + (r @ inv_cov @ delta_2prime[i,j])

        l,d,p = linalg.ldl(L_2prime, lower=True)
        l2,d2,p2 = linalg.ldl(L_2prime + 1/2 * np.matmul(L_2prime, np.matmul(covariance, L_2prime)), lower=True)
        y = linalg.solve(l, L_prime, lower=True)
        y2 = linalg.solve(l2, L_prime, lower=True)
        chi2_marg1 = -(np.dot(y, np.dot(np.linalg.inv(d),y)) - np.dot(y2, np.dot(np.linalg.inv(d2),y2)))/2
        # Two options:
        # 1) Tr(ln(1+\Sigma_cal * L''))
        # 2) ln(det(ln(1+\Sigma_cal * L''))
        # Option 2) is significantly faster.
        # chi2_full_marg2 = np.trace(scipy.linalg.logm(np.identity(self.nfitparameters) + np.matmul(self.calibration_matrix, self.L_2prime) / 2.))
        # sign, chi2_full_marg2 = np.linalg.slogdet(np.identity(self.nfitparameters) + np.matmul(self.calibration_matrix, self.L_2prime) / 2.)
        # Use SVD to calculate ln(det) because of numerical stability
        chi2_marg2 = np.sum(np.log(np.linalg.svd(np.identity(n_tomo*n_comp) + np.matmul(covariance, L_2prime)/2, compute_uv=False)))
        chi2 = chi2_fid + chi2_marg1 + chi2_marg2
        if chi2<0:
            print('chi2 < 0 !')
            chi2 = 2e12
        ln_like = -0.5*chi2
        block[names.data_vector, like_name+"_CHI2"] = chi2
        block['comb', 'reconstructed_theory_vector'] =  _mu[sorter]
        block['comb', 'bin1'] =   np.array(bin1, dtype=int)[sorter]
        block['comb', 'bin2'] =   np.array(bin2, dtype=int)[sorter]
        block['comb', 'angbin'] =  np.array(ang, dtype=int)[sorter]
        import os
        import h5py
        out = 'datablock_BS_comb_handle_pp_shear_CCLv2/comb'
        os.makedirs(out, exist_ok=True)
        with h5py.File(os.path.join(out, 'delta_primes.h5'), 'w') as fil:
            fil.create_dataset('delta_prime', data=delta_prime[:, sorter])
            fil.create_dataset('delta_2prime', data=delta_2prime[:, :, sorter])
            fil.close()
    else:
        ln_like = -0.5*chi2_fid
        block[names.data_vector, like_name+"_CHI2"] = chi2_fid
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def clean(config):
    pass
