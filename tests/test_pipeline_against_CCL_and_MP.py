import os
import warnings

import numpy as np

import cosmosis.runtime.config
import cosmosis.runtime.pipeline
import cosmosis.datablock

import pyccl as ccl

pi = np.pi

PIPELINE_FILE = "runs/config/KV450_no_sys.ini"

PATH_TO_MONTEPYTHON = "/home/fkoehlin/soft/montepython_private/"

# this is the root folder to which MontePython will save out the MCMC output
PATH_TO_MP_OUTPUT = "runs/pipeline_tests/montepython/"

# those are the folders where the data vectors are stored and to which Monte
# Python will save out the Cls and theory vectors!
PATH_TO_MP_DATA = {"K1K_COSEBIs": "data/KV450/combined/",
                   "K1K_BandPowers": "modules/scale_cuts/test_files/"}

def set_up_MontePython(likelihood):
    
    import shutil
    
    # copy likelihoods from KCAP to MontePython
    likelihood_dir_mp = os.path.join(PATH_TO_MONTEPYTHON, 'montepython/likelihoods/' + likelihood)
    likelihood_dir_kcap = os.path.join('montepython/', likelihood)
    # easy case: folders don't exist yet:
    #if not os.path.isdir(likelihood_dir_mp):
        # copy likelihood-folder contents from KCAP to MP:
    os.makedirs(likelihood_dir_mp, exist_ok=True)
    shutil.copy(os.path.join(likelihood_dir_kcap, '__init__.py'), likelihood_dir_mp)
    shutil.copy(os.path.join(likelihood_dir_kcap, likelihood  + '.data.benchmark'), os.path.join(likelihood_dir_mp, likelihood + '.data'))
    #else:
    #    warnings.warn("Trying to copy likelihood from KCAP to MP but folder does already exist! \n Check and rename manually then try again!")
    #    exit()

    print('Copied {:} likelihood with benchmark-settings to specified MontePython-folder.'.format(likelihood))
    
    return

def test_no_sys_pipeline(plot=True):
    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "ln_1e10_A_s" : 3.0,
                                               "omega_k" : 0.0,
                                               "w"       : -1.0,
                                               "mnu"     : 0.06},
                  "halo_model_parameters"   : {"A"       : 2.0,}}


    pipeline_ini = cosmosis.runtime.config.Inifile(PIPELINE_FILE)
    values_ini = cosmosis.runtime.config.Inifile(None)
    values_ini.read_dict(param_dict)

    pipeline = cosmosis.runtime.pipeline.LikelihoodPipeline(pipeline_ini, values=values_ini)

    data = pipeline.run_parameters([])

    ccl_cosmo = ccl.Cosmology(Omega_c=data["cosmological_parameters", "omega_c"], 
                              Omega_b=data["cosmological_parameters", "omega_b"],
                              Omega_k=data["cosmological_parameters", "omega_k"],
                              h=data["cosmological_parameters", "h0"], 
                              n_s=data["cosmological_parameters", "n_s"], 
                              A_s=data["cosmological_parameters", "a_s"],
                              Neff=data["cosmological_parameters", "n_eff"], 
                              m_nu=data["cosmological_parameters", "mnu"],
                              w0=data["cosmological_parameters", "w"],
                             )

    print("Loading n(z) and creating CCL tracers.")
    nofz_files = pipeline_ini.get("load_nz", "filepath").split(" ")
    tracers = []
    for nofz_file in nofz_files:
        z, nofz = np.loadtxt(nofz_file, unpack=True)
        z += (z[1]-z[0])/2
        tracers.append(ccl.WeakLensingTracer(cosmo=ccl_cosmo, dndz=(z, nofz)))
    n_tomo_bin = len(tracers)

    print("Comparing P(k) at z=0")
    h = data["cosmological_parameters", "h0"]
    k_lin = data["matter_power_lin", "k_h"]
    k_nl = data["matter_power_nl", "k_h"]
    pk_lin_ccl = ccl.linear_matter_power(ccl_cosmo, k_lin*h, 1.0)*h**3
    pk_nl_ccl = ccl.nonlin_matter_power(ccl_cosmo, k_nl*h, 1.0)*h**3

    frac_diff_pk_lin = data["matter_power_lin", "p_k"][0]/pk_lin_ccl - 1
    frac_diff_pk_nl = data["matter_power_nl", "p_k"][0]/pk_nl_ccl - 1
    print(f"Maximum fractional difference in linear P(k) at z=0 for k<5 h/Mpc: {max(np.abs(frac_diff_pk_lin[k_lin < 5.0]))}")
    print(f"Maximum fractional difference in non-linear P(k) at z=0 for k<5 h/Mpc: {max(np.abs(frac_diff_pk_nl[k_nl < 5.0]))}")

    ell = data["shear_cl", "ell"]
    print(ell.min(), ell.max(), len(ell))
    
    print("Setting up MontePython now.")
    
    import subprocess
    
    statistics_mp = {}
    for idx, likelihood in enumerate(["K1K_COSEBIs", "K1K_BandPowers"]):
        
        set_up_MontePython(likelihood)
        
        # we need to set some paths...
        path_to_param_file = os.path.join("montepython/", "{:}/INPUT/{:}_Benchmark.param".format(likelihood, likelihood[4:]))
        
        path_to_mp_output = os.path.join(PATH_TO_MP_OUTPUT, likelihood) 
        # TODO: unfortunately, MP saves out theory vectors and Cls to folder
        # from which it reads in its data and NOT to the MCMC output folder...
        path_to_mp_input = os.path.join(PATH_TO_MP_DATA[likelihood])
        
        # the call to MontePython:
        cmd = "python {:} run -p {:} -o {:} --conf {:} -N 1".format(
                os.path.join(PATH_TO_MONTEPYTHON, "montepython/MontePython.py"), 
                path_to_param_file, path_to_mp_output, 
                os.path.join(PATH_TO_MONTEPYTHON, "default.conf"))
        subprocess.call(cmd, shell=True)
        
        # we only need the Cls from one MP likelihood:
        if idx == 0:
            print("Loading and re-ordering Cls from MP.")
            
            fname = os.path.join(path_to_mp_input, 'Cls_tot.txt')
            data_mp = np.loadtxt(fname)
            mp_cl_raw = data_mp[:, 1:]
            # bring raw Cls from MP into same sorting order:
            mp_cl_my_sorting = {}
            idx_unique = 0
            for i in range(n_tomo_bin):
                mp_cl_my_sorting[i] = {}
                for j in range(i, n_tomo_bin):
                    #print(f"Bin {i+1}-{j+1}")
                    mp_cl_my_sorting[i][j] = mp_cl_raw[:, idx_unique]
                    idx_unique += 1
            
            mp_cl = {}
            for i in range(n_tomo_bin):
                mp_cl[i] = {}
                for j in range(i+1):
                    print(f"Bin {i+1}-{j+1} = Bin {j+1}-{i+1}")
                    try:
                        mp_cl[i][j] = mp_cl_my_sorting[j][i]
                    except:
                        mp_cl[i][j] = mp_cl_my_sorting[i][j]
        
        fname = os.path.join(path_to_mp_input, "{:}_theory.ascii".format(likelihood[4:]))
        print("Loaded theory vector from MontePython' {:} from: \n {:} \n".format(likelihood, fname))
        statistics_mp[likelihood] = np.loadtxt(fname)
        print(statistics_mp)
        print(statistics_mp["K1K_COSEBIs"])
        
    theta =  data["shear_xi_plus", "theta"]*180/pi
    print("Running CCL for Cls and xis.")
    ccl_cl = {}
    ccl_xi = {}
    for i in range(n_tomo_bin):
        ccl_cl[i] = {}
        ccl_xi[i] = {}
        for j in range(i+1):
            print(f"Bin {i+1}-{j+1}")
            ccl_cl[i][j] = ccl.angular_cl(cosmo=ccl_cosmo, cltracer1=tracers[i], cltracer2=tracers[j], ell=ell)
            
            ell_for_xi = np.concatenate((np.arange(2,100, dtype=float), np.logspace(2,np.log10(ell[-1]), 300)))
            cl_for_xi = ccl.angular_cl(cosmo=ccl_cosmo, cltracer1=tracers[i], cltracer2=tracers[j], ell=ell_for_xi)
            cl_for_xi /= ((ell_for_xi-1)*ell_for_xi*(ell_for_xi+1)*(ell_for_xi+2)/(ell_for_xi+1/2)**4)
            xip = ccl.correlation(cosmo=ccl_cosmo, ell=ell_for_xi, C_ell=cl_for_xi, theta=theta, corr_type="L+")
            xim = ccl.correlation(cosmo=ccl_cosmo, ell=ell_for_xi, C_ell=cl_for_xi, theta=theta, corr_type="L-")
            ccl_xi[i][j] = xip, xim

    ell_range = 10, 10000
    theta_range = 0.05, 10.0
    ell_mask = (ell_range[0] < ell) & (ell < ell_range[1])
    theta_mask = (theta_range[0] < theta) & (theta < theta_range[1])
    for i in range(n_tomo_bin):
        for j in range(i+1):
            frac_diff_cl_KCAP_CCL = data["shear_cl", f"bin_{i+1}_{j+1}"]/ccl_cl[i][j] - 1
            frac_diff_cl_KCAP_MP = data["shear_cl", f"bin_{i+1}_{j+1}"]/mp_cl[i][j] - 1
            print(f"Maximum fractional difference between KCAP and CCL in Cl (ell>10) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_cl_KCAP_CCL[ell_mask]))}")
            print(f"Maximum fractional difference between KCAP and MP in Cl (ell>10) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_cl_KCAP_MP[ell_mask]))}")

            frac_diff_xi_plus = data["shear_xi_plus", f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][0] - 1
            print(f"Maximum fractional difference in xi_p (theta<10 deg) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_xi_plus[theta_mask]))}")

            frac_diff_xi_minus = data["shear_xi_minus", f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][1] - 1
            print(f"Maximum fractional difference in xi_m (theta>0.05 deg) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_xi_minus[theta_mask]))}")

    if plot:
        import matplotlib.pyplot as plt

        print("Plotting P(k)")
 
        fig, ax = plt.subplots(2, 1, sharex=True, sharey=False,
                               figsize=(5, 5))
        fig.subplots_adjust(hspace=0, wspace=0)
 
        ax[0].loglog(k_lin, pk_lin_ccl, ls="--", label=f"CCL linear P(k)")
        ax[0].loglog(k_lin, data["matter_power_lin", "p_k"][0], ls="--", label=f"CosmoSIS linear P(k)")
        ax[0].loglog(k_nl, pk_nl_ccl, label=f"CCL non-linear P(k)")
        ax[0].loglog(k_nl, data["matter_power_nl", "p_k"][0], label=f"CosmoSIS non-linear P(k)")
 
        ax[1].semilogx(k_lin, data["matter_power_lin", "p_k"][0]/pk_lin_ccl-1, ls="--")
        ax[1].semilogx(k_nl, data["matter_power_nl", "p_k"][0]/pk_nl_ccl-1)
 
        ax[0].legend(frameon=False)
        ax[0].set_ylabel("$P(k)$ [$h^{-3}$ Mpc$^{3}$]")
        ax[1].set_ylim(-0.05, 0.05)
        ax[1].set_xlabel("$k$ [$h$ Mpc$^{-1}$]")
        ax[1].set_ylabel("$\Delta P(k)/P(k)$ [$h$ Mpc$^{-1}$]")
         
        fig.suptitle("KCAP vs. CCL, P(k)")
        fig.savefig("KV450_pofk_kcap_vs_ccl.pdf") 
 
        ell_plot_lim = ell_range

        print("Plotting Cls")

        fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                                figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
        fig.subplots_adjust(hspace=0, wspace=0)

        u = ell**2
        for i in range(n_tomo_bin):
            for j in range(n_tomo_bin):
                if j > i:
                    ax[i][j].axis("off")
                else:
                    ax[i][j].loglog(ell, u*ccl_cl[i][j], label=f"CCL bin {i+1}-{j+1}")
                    ax[i][j].loglog(ell, u*data["shear_cl", f"bin_{i+1}_{j+1}"], label=f"CosmoSIS")
                    ax[i][j].loglog(ell, u*mp_cl[i][j], label=f"MP")

                    ax[i][j].legend(fontsize="small", frameon=False)
                    ax[i][j].set_xlim(*ell_plot_lim)

        for p in ax[-1]:
            p.set_xlabel(r"$\ell$")
        for p in ax:
            p[0].set_ylabel(r"$\ell^2\ C_\ell$")

        fig.suptitle("KCAP vs. CCL vs. MP, Cls")
        fig.savefig("KV450_Cl_kcap_vs_ccl_vs_mp.pdf")

        print("Plotting Cl fractional differences.")
        fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                            figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
        fig.subplots_adjust(hspace=0, wspace=0)
        for i in range(n_tomo_bin):
            for j in range(n_tomo_bin):
                if j > i:
                    ax[i][j].axis("off")
                else:
                    label = f"bin {i+1}-{j+1}"
                    ax[i][j].text(0.75, 0.80, label, horizontalalignment='center', transform=ax[i][j].transAxes)
                    ax[i][j].semilogx(ell, data["shear_cl", f"bin_{i+1}_{j+1}"]/ccl_cl[i][j]-1, label=f"KCAP vs. CCL")
                    ax[i][j].semilogx(ell, data["shear_cl", f"bin_{i+1}_{j+1}"]/mp_cl[i][j]-1, label="KCAP vs. MP")
                    #ax[i][j].legend(fontsize="small", frameon=False)
                
                    ax[i][j].axhline(0., ls='-', color='gray')
                    for val in [0.01, 0.05]:
                        ax[i][j].axhline(val, ls=':', color='gray')
                        ax[i][j].axhline(-val, ls=':', color='gray')
                    
                ax[i][j].set_xlim(*ell_plot_lim)
                ax[i][j].set_ylim(-0.1, 0.1)
        
        ax[1][1].legend(fontsize="small", frameon=False, bbox_to_anchor=(.95, 1.05))
        
        for p in ax[-1]:
            p.set_xlabel(r"$\ell$")
        for p in ax:
            p[0].set_ylabel(r"$|\Delta C_\ell|/C_\ell$")

        fig.suptitle("KCAP vs. CCL vs. MP, Cls")
        fig.savefig("KV450_Cl_kcap_vs_ccl_vs_mp_frac_diff.pdf") 

        print("Plotting xi fractional differences.")
        theta_plot_lim = theta_range
        fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                                figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
        fig.subplots_adjust(hspace=0, wspace=0)
        for i in range(n_tomo_bin):
            for j in range(n_tomo_bin):
                if j > i:
                    ax[i][j].axis("off")
                else:
                    ax[i][j].semilogx(theta, data["shear_xi_plus", f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][0]-1, label=f"xip bin {i+1}-{j+1}")
                    ax[i][j].semilogx(theta, data["shear_xi_minus", f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][1]-1, label=f"xim bin {i+1}-{j+1}")
                    ax[i][j].legend(fontsize="small", frameon=False)
        
                ax[i][j].set_xlim(*theta_plot_lim)
                ax[i][j].set_ylim(-0.1, 0.1)

        for p in ax[-1]:
            p.set_xlabel(r"$\theta$ [deg]")
        for p in ax:
            p[0].set_ylabel(r"$|\Delta \xi|/\xi$")

        fig.suptitle("KCAP vs. CCL, xis")
        fig.savefig("KV450_xi_kcap_vs_ccl_frac_diff.pdf") 

        print("Plotting COSEBIs fractional differences.")
        # TODO: this range should be taken from CosmoSIS keywords!!!
        n_modes = 5
        n_tomo_corrs = n_tomo_bin * (n_tomo_bin + 1) // 2
        modes = np.arange(1, n_modes + 1)
        cosebis_mp = statistics_mp["K1K_COSEBIs"].reshape(n_tomo_corrs, n_modes)
        fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                               figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
        fig.subplots_adjust(hspace=0, wspace=0)
        idx_corr = 0
        for i in range(n_tomo_bin):
            for j in range(n_tomo_bin):
                if j > i:
                    ax[i][j].axis("off")
                else:
                    # TODO: add also output from CosmoSIS!
                    #ax[i][j].semilogx(modes, data["cosebis", f"bin_{i+1}_{j+1}"] / cosebis_mp[idx_corr, :] -1, label=f"COSEBIs bin {i+1}-{j+1}")
                    ax[i][j].plot(modes, cosebis_mp[idx_corr, :] / cosebis_mp[idx_corr, :] -1, label=f"COSEBIs bin {i+1}-{j+1}")
                    ax[i][j].legend(fontsize="small", frameon=False)
                    idx_corr += 1
                
                ax[i][j].set_xlim([0., 6.])
                ax[i][j].set_ylim([-0.1, 0.1])
                    
        for p in ax[-1]:
            p.set_xlabel(r"$n$")
        for p in ax:
            p[0].set_ylabel(r"$|\Delta E_n|/E_n$")
            
        fig.suptitle("KCAP vs. MP, COSEBIs")
        fig.savefig("KV450_cosebis_kcap_vs_mp_frac_diff.pdf") 

        print("Plotting BandPowers fractional differences.")
        # TODO: this range should be taken from CosmoSIS keywords!!!
        n_ell_bins = 8
        n_tomo_corrs = n_tomo_bin * (n_tomo_bin + 1) // 2
        plot_ell = np.logspace(np.log10(100.), np.log10(1500.), n_ell_bins)
        bandpowers_mp = statistics_mp["K1K_BandPowers"].reshape(n_tomo_corrs, n_ell_bins)
        fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                               figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
        fig.subplots_adjust(hspace=0, wspace=0)
        idx_corr = 0
        for i in range(n_tomo_bin):
            for j in range(n_tomo_bin):
                if j > i:
                    ax[i][j].axis("off")
                else:
                    # TODO: add also output from CosmoSIS!
                    #ax[i][j].semilogx(modes, data["bandpower", f"bin_{i+1}_{j+1}"] / bandpowers_mp[idx_corr, :] -1, label=f"BandPowers bin {i+1}-{j+1}")
                    ax[i][j].semilogx(plot_ell, bandpowers_mp[idx_corr, :] / bandpowers_mp[idx_corr, :] -1, label=f"BandPowers bin {i+1}-{j+1}")
                    ax[i][j].legend(fontsize="small", frameon=False)
                    idx_corr += 1
                
                ax[i][j].set_xlim([100., 1500.])
                ax[i][j].set_ylim([-0.1, 0.1])
                    
        for p in ax[-1]:
            p.set_xlabel(r"$\ell$")
        for p in ax:
            p[0].set_ylabel(r"$|\Delta PeeE|/PeeE $")
            
        fig.suptitle("KCAP vs. MP, BandPowers")
        fig.savefig("MOCK_bandpowers_kcap_vs_mp_frac_diff.pdf")

        plt.show()

if __name__ == "__main__":
    if not os.path.isdir("cosmosis-standard-library"):
        warnings.warn("cosmosis-standard-library not in cwd. Are you running this test from the kcap root directory?")
    test_no_sys_pipeline()