import os
import glob
import re
import sys

import numpy as np
pi = np.pi

import pyccl as ccl

def load_cosmosis_params(output_path, section="cosmological_parameters"):
    params = {}
    with open(os.path.join(output_path, section, "values.txt"), "r") as f:
        for l in f.readlines():
            m = re.match("^([0-9A-Za-z_]+) = (.+)\w?$", l)
            if m is None:
                print(f"Ill formed entry: {l[:-1]}.")
                continue
            key, value = m.groups()
            if re.match("^[0-9\-]+$", value) is not None:
                value = int(value)
            elif re.match("^[0-9eE\-\.]+$", value) is not None:
                value = float(value)
            elif value.lower() == "true":
                value = True
            elif value.lower() == "false":
                value = False
            params[key] = value
            
    return params

def load_cosmosis_2pt(output_path, x_name="ell", y_name="shear", suffix="cl"):
    cosmosis_2pt = {}

    name = y_name + "_" + suffix
    for path in glob.glob(os.path.join(output_path, f"{name}*")):
        section = os.path.split(path)[1]
        cosmosis_2pt[section] = {}

        values = load_cosmosis_params(os.path.split(path)[0], section=section)
        cosmosis_2pt[section] = {**values}
        
        if os.path.isfile(os.path.join(path, x_name+".txt")):
            x_file = True
            x = np.loadtxt(os.path.join(path, x_name+".txt"))
        else:
            x_file = False
            
        for bin_path in glob.glob(os.path.join(path, "*_*_*")):
            key = os.path.split(bin_path)[1]
            key = key[:key.rfind(".txt")]
            data = np.loadtxt(bin_path).T
            if x_file:
                y = data
            else:
                x = data[0]
                y = data[1]
            
            cosmosis_2pt[section][key] = y
            
        cosmosis_2pt[section][x_name] = x
        
    return cosmosis_2pt

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    cosmo_params = load_cosmosis_params("../../runs/output/KV450/output_no_sys/")
    
    cosmosis_pk = {"matter_power_lin": {"p_k" : np.loadtxt("../../runs/output/KV450/output_no_sys/matter_power_lin/p_k.txt"),
                                        "k_h" : np.loadtxt("../../runs/output/KV450/output_no_sys/matter_power_lin/k_h.txt"),
                                        },
                   "matter_power_nl":  {"p_k" : np.loadtxt("../../runs/output/KV450/output_no_sys/matter_power_nl/p_k.txt"),
                                        "k_h" : np.loadtxt("../../runs/output/KV450/output_no_sys/matter_power_nl/k_h.txt"),
                                        },
                  }
    k_lin_cosmosis = cosmosis_pk["matter_power_lin"]["k_h"]
    k_nl_cosmosis = cosmosis_pk["matter_power_nl"]["k_h"]

    cosmosis_cl = load_cosmosis_2pt("../../runs/output/KV450/output_no_sys/",
                                    suffix="cl")
    ell_cosmosis = cosmosis_cl["shear_cl"]["ell"]

    cosmosis_xi = load_cosmosis_2pt("../../runs/output/KV450/output_no_sys/",
                                    suffix="xi_plus", x_name="theta")
    theta_xi_plus_cosmosis = cosmosis_xi["shear_xi_plus"]["theta"]
    
    cosmosis_xi.update(load_cosmosis_2pt("../../runs/output/KV450/output_no_sys/",
                                    suffix="xi_minus", x_name="theta"))
    theta_xi_minus_cosmosis = cosmosis_xi["shear_xi_minus"]["theta"]

    n_tomo_bin = int(cosmosis_cl["shear_cl"]["nbin_a"])

    print("Cosmological parameters")
    print(f'Omega_c:    {cosmo_params["omega_c"]}')
    print(f'Omega_b:    {cosmo_params["omega_b"]}')
    print(f'h:          {cosmo_params["h0"]}')
    print(f'n_s:        {cosmo_params["n_s"]}')
    print(f'sigma_8:    {cosmo_params["sigma_8"]}')
    print()
    print(f'n_bin:      {n_tomo_bin}')

    ccl_cosmo = ccl.Cosmology(Omega_c=cosmo_params["omega_c"], Omega_b=cosmo_params["omega_b"],
                              h=cosmo_params["h0"], n_s=cosmo_params["n_s"], A_s=cosmo_params["a_s"],
                              Neff=cosmo_params["nnu"], m_nu=cosmo_params["mnu"],
                                )
    
    k_lin = k_lin_cosmosis
    k_nl = k_nl_cosmosis
    pk_lin = ccl.linear_matter_power(ccl_cosmo, k_lin*cosmo_params["h0"], 1.0)*cosmo_params["h0"]**3
    pk_nl = ccl.nonlin_matter_power(ccl_cosmo, k_nl*cosmo_params["h0"], 1.0)*cosmo_params["h0"]**3

    print("Plotting P(k)")

    fig, ax = plt.subplots(2, 1, sharex=True, sharey=False,
                           figsize=(5, 5))
    fig.subplots_adjust(hspace=0, wspace=0)

    ax[0].loglog(k_lin, pk_lin, ls="--", label=f"CCL linear P(k)")
    ax[0].loglog(k_lin, cosmosis_pk["matter_power_lin"]["p_k"][0], ls="--", label=f"CosmoSIS linear P(k)")
    ax[0].loglog(k_nl, pk_nl, label=f"CCL non-linear P(k)")
    ax[0].loglog(k_nl, cosmosis_pk["matter_power_nl"]["p_k"][0], label=f"CosmoSIS non-linear P(k)")

    ax[1].semilogx(k_lin, cosmosis_pk["matter_power_lin"]["p_k"][0]/pk_lin-1, ls="--")
    ax[1].semilogx(k_nl, cosmosis_pk["matter_power_nl"]["p_k"][0]/pk_nl-1)

    ax[0].legend(frameon=False)
    ax[0].set_ylabel("$P(k)$ [$h^{-3}$ Mpc$^{3}$]")
    ax[1].set_ylim(-0.05, 0.05)
    ax[1].set_xlabel("$k$ [$h$ Mpc$^{-1}$]")
    ax[1].set_ylabel("$\Delta P(k)/P(k)$ [$h$ Mpc$^{-1}$]")
    
    fig.suptitle("KCAP vs CCL, P(k)")
    fig.savefig("KV450_pofk_kcap_vs_ccl.pdf") 

    nofz_files = ["../../data/KV450/nofz/Nz_DIR_z0.1t0.3.asc",
                  "../../data/KV450/nofz/Nz_DIR_z0.3t0.5.asc",
                  "../../data/KV450/nofz/Nz_DIR_z0.5t0.7.asc",
                  "../../data/KV450/nofz/Nz_DIR_z0.7t0.9.asc",
                  "../../data/KV450/nofz/Nz_DIR_z0.9t1.2.asc",]

    print("Loading n(z) and creating CCL tracers.")
    tracers = []
    for nofz_file in nofz_files:
        z, nofz = np.loadtxt(nofz_file, unpack=True)
        z += (z[1]-z[0])/2
        tracers.append(ccl.WeakLensingTracer(cosmo=ccl_cosmo, dndz=(z, nofz)))

    # ell_min = 10
    # ell_max = 1e4
    # n_ell = 100
    # ell = np.logspace(np.log10(ell_min), np.log10(ell_max), n_ell, endpoint=True)
    ell = ell_cosmosis
    theta =  theta_xi_plus_cosmosis*180/pi
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

    print("Plotting Cls")
    ell_plot_lim = 10, 5000

    fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                           figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
    fig.subplots_adjust(hspace=0, wspace=0)

    u_cosmosis = ell_cosmosis**2

    u = ell**2
    for i in range(n_tomo_bin):
        for j in range(n_tomo_bin):
            if j > i:
                ax[i][j].axis("off")
            else:
                ax[i][j].loglog(ell, u*ccl_cl[i][j], label=f"CCL bin {i+1}-{j+1}")
                ax[i][j].loglog(ell_cosmosis, u_cosmosis*cosmosis_cl["shear_cl"][f"bin_{i+1}_{j+1}"], label=f"CosmoSIS")

                ax[i][j].legend(fontsize="small", frameon=False)
                ax[i][j].set_xlim(*ell_plot_lim)

    for p in ax[-1]:
        p.set_xlabel(r"$\ell$")
    for p in ax:
        p[0].set_ylabel(r"$\ell^2\ C_\ell$")

    fig.suptitle("KCAP vs CCL, Cls")
    fig.savefig("KV450_Cl_kcap_vs_ccl.pdf")

    print("Plotting Cl fractional differences.")
    fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                           figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
    fig.subplots_adjust(hspace=0, wspace=0)
    for i in range(n_tomo_bin):
        for j in range(n_tomo_bin):
            if j > i:
                ax[i][j].axis("off")
            else:
                ax[i][j].semilogx(ell, cosmosis_cl["shear_cl"][f"bin_{i+1}_{j+1}"]/ccl_cl[i][j]-1, label=f"bin {i+1}-{j+1}")
                ax[i][j].legend(fontsize="small", frameon=False)
    
            ax[i][j].set_xlim(*ell_plot_lim)
            ax[i][j].set_ylim(-0.1, 0.1)

    for p in ax[-1]:
        p.set_xlabel(r"$\ell$")
    for p in ax:
        p[0].set_ylabel(r"$|\Delta C_\ell|/C_\ell$")

    fig.suptitle("KCAP vs CCL, Cls")
    fig.savefig("KV450_Cl_kcap_vs_ccl_frac_diff.pdf") 


    print("Plotting xi fractional differences.")
    theta_plot_lim = 1/60, 5.0
    fig, ax = plt.subplots(n_tomo_bin, n_tomo_bin, sharex=True, sharey=True,
                           figsize=(2*n_tomo_bin, 1.5*n_tomo_bin))
    fig.subplots_adjust(hspace=0, wspace=0)
    for i in range(n_tomo_bin):
        for j in range(n_tomo_bin):
            if j > i:
                ax[i][j].axis("off")
            else:
                ax[i][j].semilogx(theta, cosmosis_xi["shear_xi_plus"][f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][0]-1, label=f"xip bin {i+1}-{j+1}")
                ax[i][j].semilogx(theta, cosmosis_xi["shear_xi_minus"][f"bin_{i+1}_{j+1}"]/ccl_xi[i][j][1]-1, label=f"xim bin {i+1}-{j+1}")
                ax[i][j].legend(fontsize="small", frameon=False)
    
            ax[i][j].set_xlim(*theta_plot_lim)
            ax[i][j].set_ylim(-0.1, 0.1)

    for p in ax[-1]:
        p.set_xlabel(r"$\theta$ [deg]")
    for p in ax:
        p[0].set_ylabel(r"$|\Delta \xi|/\xi$")

    fig.suptitle("KCAP vs CCL, xis")
    fig.savefig("KV450_xi_kcap_vs_ccl_frac_diff.pdf") 

    

    


