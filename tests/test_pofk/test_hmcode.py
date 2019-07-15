import numpy as np
import scipy.interpolate

def pofk_interpolator(pofk, k, z=None):
    if z is None:
        intp = scipy.interpolate.InterpolatedUnivariateSpline(np.log(k), np.log(pofk))
        return lambda k: np.exp(intp(np.log(k))).squeeze()
    else:
        intp = scipy.interpolate.RectBivariateSpline(z, np.log(k), np.log(pofk))
        return lambda k, z: np.exp(intp(z, np.log(k), grid=True)).squeeze()

def test_hmcode(plot=False):
    import camb
    cosmological_parameters = {"ombh2" : 0.2228089E-01,
                            "omch2" : 0.1181838E+00,
                            "h0"    : 0.679565,
                            "A_s"   : 2.170229e-9,
                            "n_s"   : 0.9679701,
                            "tau"    : 0.0739,
                            "omnuh2" : 0.00065,
                            "omega_nu" : 0.00139699,
                            "massive_nu" : 1,
                            "massless_nu" : 2.046,
                            }

    pars = camb.set_params(H0=cosmological_parameters["h0"]*100, 
                        ombh2=cosmological_parameters["ombh2"], 
                        omch2=cosmological_parameters["omch2"],
                        omnuh2=cosmological_parameters["omnuh2"],
                        num_nu_massive=cosmological_parameters["massive_nu"],
                        num_nu_massless=cosmological_parameters["massless_nu"],
                        As=cosmological_parameters["A_s"], 
                        ns=cosmological_parameters["n_s"],
                        r=0,
                        max_eta_k=1000,
                        lmax=2500,
                        )

    pars.set_matter_power(redshifts=np.linspace(0, 4, 100), kmax=100.0)
    pars.NonLinearModel.set_params(halofit_version='mead')

    results = camb.get_results(pars)

    kh_camb, z_camb, pk_camb_mead = results.get_nonlinear_matter_power_spectrum(params=pars)
    kh_lin_camb, z_lin_camb, pk_lin_camb = results.get_linear_matter_power_spectrum(params=pars)

    p_k_camb = pofk_interpolator(pk_camb_mead, kh_camb, z_camb)
    p_k_lin_camb = pofk_interpolator(pk_lin_camb, kh_lin_camb, z_lin_camb)

    if plot:
        import matplotlib
        import matplotlib.pyplot as plt

        z = np.loadtxt("cosmosis_output/matter_power_nl/z.txt")
        k = np.loadtxt("cosmosis_output/matter_power_nl/k_h.txt")
        p_k_nonlin = np.loadtxt("cosmosis_output/matter_power_nl/p_k.txt")

        z_lin = np.loadtxt("cosmosis_output/matter_power_lin/z.txt")
        k_lin = np.loadtxt("cosmosis_output/matter_power_lin/k_h.txt")
        p_k_lin = np.loadtxt("cosmosis_output/matter_power_lin/p_k.txt")

        cmap = plt.get_cmap("magma_r")

        fig, ax = plt.subplots(2, 1, figsize=(5,6))
        fig.subplots_adjust(hspace=0.5, left=0.2, right=0.95)


        cb_ax = matplotlib.colorbar.make_axes(ax)
        norm = matplotlib.colors.Normalize(vmin=z[0], vmax=z[-1])
        cb1 = matplotlib.colorbar.ColorbarBase(cb_ax[0], cmap=cmap,
                                        norm=norm, **cb_ax[1])
        cb1.set_label('z')

        _ = [ax[0].semilogx(k_lin, p_k_lin[i]/p_k_lin_camb(k_lin, z_lin[i]) - 1, c=cmap(i/len(z_lin))) for i in range(len(z_lin))]
        ax[0].set_title("Linear power spectrum")
        ax[0].set_xlabel("k [h/Mpc]")
        ax[0].set_ylabel("CosmoSIS/CAMB-1")

        _ = [ax[1].semilogx(k, p_k_nonlin[i]/p_k_camb(k ,z[i]) - 1, c=cmap(i/len(z))) for i in range(len(z))]
        ax[1].set_title("Non-linear power spectrum")
        ax[1].set_xlabel("k [h/Mpc]")
        ax[1].set_ylabel("CosmoSIS/CAMB-1")

    plt.show()

def test_hmcode_no_camb(plot=False):
    z = np.loadtxt("cosmosis_output/matter_power_nl/z.txt")
    k = np.loadtxt("cosmosis_output/matter_power_nl/k_h.txt")
    p_k_nonlin = np.loadtxt("cosmosis_output/matter_power_nl/p_k.txt")

    z_camb = np.loadtxt("camb_output/z.txt")
    k_camb = np.loadtxt("camb_output/k_h.txt")
    p_k_camb = np.loadtxt("camb_output/pk_mead.txt")

    p_k_camb = pofk_interpolator(p_k_camb, k_camb, z_camb)
    assert np.allclose(p_k_nonlin, p_k_camb(k, z), rtol=1e-2)

    if plot:
        import matplotlib.pyplot as plt

        cmap = plt.get_cmap("magma_r")
        _ = [plt.semilogx(k, p_k_nonlin[i]/p_k_camb(k,z[i]) - 1, c=cmap(i/len(z))) for i in range(len(z))]
        plt.show()
    



if __name__ == "__main__":
    #test_hmcode(plot=True)
    test_hmcode_no_camb()