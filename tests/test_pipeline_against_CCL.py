import os
import warnings

import numpy as np

import cosmosis.runtime.config
import cosmosis.runtime.pipeline
import cosmosis.datablock

import pyccl as ccl

pi = np.pi

PIPELINE_FILE = "runs/config/KV450_no_sys.ini"

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
            frac_diff_cl = data["shear_cl", f"bin_{i+1}_{j+1}"]/ccl_cl[i][j] - 1
            print(f"Maximum fractional difference in Cl (ell>10) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_cl[ell_mask]))}")

            frac_diff_xi_plus = data["shear_xi_plus", f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][0] - 1
            print(f"Maximum fractional difference in xi_p (theta<10 deg) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_xi_plus[theta_mask]))}")

            frac_diff_xi_minus = data["shear_xi_minus", f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][1] - 1
            print(f"Maximum fractional difference in xi_m (theta>0.05 deg) for bin {i+1}-{j+1}: {max(np.abs(frac_diff_xi_plus[theta_mask]))}")


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
        
        fig.suptitle("KCAP vs CCL, P(k)")
        # fig.savefig("KV450_pofk_kcap_vs_ccl.pdf") 


        print("Plotting Cls")
        ell_plot_lim = ell_range

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

                    ax[i][j].legend(fontsize="small", frameon=False)
                    ax[i][j].set_xlim(*ell_plot_lim)

        for p in ax[-1]:
            p.set_xlabel(r"$\ell$")
        for p in ax:
            p[0].set_ylabel(r"$\ell^2\ C_\ell$")

        fig.suptitle("KCAP vs CCL, Cls")
        #fig.savefig("KV450_Cl_kcap_vs_ccl.pdf")

        print("Plotting Cl fractional differences.")
        fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                            figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
        fig.subplots_adjust(hspace=0, wspace=0)
        for i in range(n_tomo_bin):
            for j in range(n_tomo_bin):
                if j > i:
                    ax[i][j].axis("off")
                else:
                    ax[i][j].semilogx(ell, data["shear_cl", f"bin_{i+1}_{j+1}"]/ccl_cl[i][j]-1, label=f"bin {i+1}-{j+1}")
                    ax[i][j].legend(fontsize="small", frameon=False)
        
                ax[i][j].set_xlim(*ell_plot_lim)
                ax[i][j].set_ylim(-0.1, 0.1)

        for p in ax[-1]:
            p.set_xlabel(r"$\ell$")
        for p in ax:
            p[0].set_ylabel(r"$|\Delta C_\ell|/C_\ell$")

        fig.suptitle("KCAP vs CCL, Cls")
        #fig.savefig("KV450_Cl_kcap_vs_ccl_frac_diff.pdf") 

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

        fig.suptitle("KCAP vs CCL, xis")
        #fig.savefig("KV450_xi_kcap_vs_ccl_frac_diff.pdf") 
        
        plt.show()

    

if __name__ == "__main__":
    if not os.path.isdir("cosmosis-standard-library"):
        warnings.warn("cosmosis-standard-library not in cwd. Are you running this test from the kcap root directory?")
    test_no_sys_pipeline()