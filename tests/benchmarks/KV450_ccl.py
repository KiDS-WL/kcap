import os
import glob

import numpy as np
pi = np.pi

import pyccl as ccl

def load_cosmosis_params(output_path, section="cosmological_parameters"):
    with open(os.path.join(output_path, section, "values.txt"), "r") as f:
        params = {line.split("=")[0].strip() : float(line.split("=")[1].strip()) for line in f.readlines()}
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

    cosmo_params = load_cosmosis_params("../../examples/output_kv450_no_sys/kv450_no_sys_output/")

    cosmosis_cl = load_cosmosis_2pt("../../examples/output_kv450_no_sys/kv450_no_sys_output/",
                                    suffix="cl")
    ell_cosmosis = cosmosis_cl["shear_cl"]["ell"]

    cosmosis_xi = load_cosmosis_2pt("../../examples/output_kv450_no_sys/kv450_no_sys_output/",
                                    suffix="xi", x_name="theta")
    theta_cosmosis = cosmosis_xi["shear_xi"]["theta"]

    ccl_cosmo = ccl.Cosmology(Omega_c=cosmo_params["omega_c"], Omega_b=cosmo_params["omega_b"],
                              h=cosmo_params["h0"], n_s=cosmo_params["n_s"], sigma8=cosmo_params["sigma_8"],
                              Neff=2.046, m_nu=0.06, mnu_type="sum")

    n_tomo_bin = int(cosmosis_cl["shear_cl"]["nbin_a"])

    print("Cosmological parameters")
    print(f'Omega_c:    {cosmo_params["omega_c"]}')
    print(f'Omega_b:    {cosmo_params["omega_b"]}')
    print(f'h:          {cosmo_params["h0"]}')
    print(f'n_s:        {cosmo_params["n_s"]}')
    print(f'sigma_8:    {cosmo_params["sigma_8"]}')
    print()
    print(f'n_bin:      {n_tomo_bin}')

    nofz_files = ["../../examples/kv450_data/nofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB01t03_blindA_Nz.asc",
                  "../../examples/kv450_data/nofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB03t05_blindA_Nz.asc",
                  "../../examples/kv450_data/nofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB05t07_blindA_Nz.asc",
                  "../../examples/kv450_data/nofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB07t09_blindA_Nz.asc",
                  "../../examples/kv450_data/nofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB09t12_blindA_Nz.asc"]

    print("Loading n(z) and creating CCL tracers.")
    tracers = []
    for nofz_file in nofz_files:
        z, nofz = np.loadtxt(nofz_file, unpack=True)
        z += (z[1]-z[0])/2
        tracers.append(ccl.ClTracerLensing(cosmo=ccl_cosmo, z=z, n=nofz, has_intrinsic_alignment=False))

    # ell_min = 10
    # ell_max = 1e4
    # n_ell = 100
    # ell = np.logspace(np.log10(ell_min), np.log10(ell_max), n_ell, endpoint=True)
    ell = ell_cosmosis
    theta =  theta_cosmosis*180/pi
    print("Running CCL for Cls and xis.")
    ccl_cl = {}
    ccl_xi = {}
    for i in range(n_tomo_bin):
        ccl_cl[i] = {}
        ccl_xi[i] = {}
        for j in range(i+1):
            print(f"Bin {i+1}-{j+1}")
            ccl_cl[i][j] = ccl.angular_cl(cosmo=ccl_cosmo, cltracer1=tracers[i], cltracer2=tracers[j], ell=ell)

            xip = ccl.correlation(cosmo=ccl_cosmo, ell=ell, C_ell=ccl_cl[i][j], theta=theta, corr_type="L+")
            xim = ccl.correlation(cosmo=ccl_cosmo, ell=ell, C_ell=ccl_cl[i][j], theta=theta, corr_type="L-")
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
                ax[i][j].semilogx(theta, cosmosis_xi["shear_xi"][f"xiplus_{i+1}_{j+1}"]/ccl_xi[i][j][0]-1, label=f"xip bin {i+1}-{j+1}")
                ax[i][j].semilogx(theta, cosmosis_xi["shear_xi"][f"ximinus_{i+1}_{j+1}"]/ccl_xi[i][j][1]-1, label=f"xim bin {i+1}-{j+1}")
                ax[i][j].legend(fontsize="small", frameon=False)
    
            ax[i][j].set_xlim(*theta_plot_lim)
            ax[i][j].set_ylim(-0.1, 0.1)

    for p in ax[-1]:
        p.set_xlabel(r"$\theta$ [deg]")
    for p in ax:
        p[0].set_ylabel(r"$|\Delta \xi|/\xi$")

    fig.suptitle("KCAP vs CCL, xis")
    fig.savefig("KV450_xi_kcap_vs_ccl_frac_diff.pdf") 

    

    


