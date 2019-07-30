import os

import numpy as np
import scipy.interpolate

import cosmosis.runtime.module

import sys
sys.path.append("../kcap")
from cosmosis_utils import dict_to_datablock

KCAP_MODULE_PATH = "../utils"

def create_pipeline(config, verbose=True):
    modules = []
    for module_name, c in config.items():
        filename = c.pop("file")
        module = cosmosis.runtime.module.Module(module_name=module_name,
                                                file_path=filename)
        module.setup(dict_to_datablock({module_name : c}))
        modules.append(module)

    def pipeline(block):
        for m in modules:
            if verbose: print(f"Running {m.name}")
            status = m.execute(block)
            if status != 0:
                raise RuntimeError(f"Module {m.name} failed at execute.")
    
    return pipeline

def pofk_interpolator(pofk, k, z=None):
    if z is None:
        intp = scipy.interpolate.InterpolatedUnivariateSpline(np.log(k), np.log(pofk))
        return lambda k: np.exp(intp(np.log(k))).squeeze()
    else:
        intp = scipy.interpolate.RectBivariateSpline(z, np.log(k), np.log(pofk))
        return lambda k, z: np.exp(intp(z, np.log(k), grid=True)).squeeze()

def test_meadcb(plot=False):
    config_camb = {     "camb"     : {"file" : os.path.join(os.path.join(os.environ["CSL_PATH"], 
                                                                    "boltzmann/pycamb/camb_interface.py")),
                                    "do_reionization"    : False,
                                    "mode"               : "transfer",
                                    "halofit_version"    : "mead",
                                    "nonlinear"          : "pk",
                                    "kmax"               : 20.0,}}

    config_meadcb = {   "camb"        : {"file" : os.path.join(os.path.join(os.environ["CSL_PATH"], 
                                                                      "boltzmann/pycamb/camb_interface.py")),
                                         "do_reionization" : False,
                                         "mode"            : "transfer",
                                         "nonlinear"       : "none",
                                         "kmax"               : 20.0,},
                        "consistency" : {"file" : os.path.join(os.path.join(os.environ["CSL_PATH"], 
                                                                      "utility/consistency/consistency_interface.py"))},
                        "meadcb"      : {"file" : os.path.join(os.path.join(os.environ["CSL_PATH"], 
                                                                      "structure/meadcb/mead_interface.so")),
                                         "one_baryon_parameter" : False,
                                         "feedback"             : False}}


    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "omega_k" : 0.0,
                                               "w"       : -0.9,
                                               "mnu"     : 0.06}}
    
    camb_pipeline = create_pipeline(config_camb)

    meadcb_pipeline = create_pipeline(config_meadcb)

    block_camb = dict_to_datablock(param_dict)
    camb_pipeline(block_camb)
    
    block_meadcb = dict_to_datablock(param_dict)
    meadcb_pipeline(block_meadcb)

    assert np.allclose(block_camb["matter_power_nl", "z"], block_meadcb["matter_power_nl", "z"])
    assert np.allclose(block_camb["matter_power_nl", "k_h"], block_meadcb["matter_power_nl", "k_h"])

    k = block_camb["matter_power_nl", "k_h"]
    z = block_camb["matter_power_nl", "z"]

    assert np.allclose(block_camb["matter_power_nl", "p_k"][z < 2.0], block_meadcb["matter_power_nl", "p_k"][z < 2.0], rtol=3e-3)
    assert np.allclose(block_camb["matter_power_nl", "p_k"][z >= 2.0], block_meadcb["matter_power_nl", "p_k"][z >= 2.0], rtol=1e-2)
    print("WARNING: Non-linear powerspectrum tolerance is set to 3e-3 for z < 2 and 1e-2 for z > 2.")

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.colorbar

        cmap = plt.get_cmap("magma_r")
        fig, ax = plt.subplots(1, 1, )
        fig.subplots_adjust(left=0.2, right=0.95)

        cb_ax = matplotlib.colorbar.make_axes(ax)
        norm = matplotlib.colors.Normalize(vmin=z[0], vmax=z[-1])
        cb1 = matplotlib.colorbar.ColorbarBase(cb_ax[0], cmap=cmap,
                                        norm=norm, **cb_ax[1])
        cb1.set_label('z')

        _ = [ax.semilogx(k, block_meadcb["matter_power_nl", "p_k"][i]/block_camb["matter_power_nl", "p_k"][i] - 1, 
                         c=cmap(i/len(z))) for i in range(len(z))]
        ax.set_title("Non-linear power spectrum")
        ax.set_xlabel("k [h/Mpc]")
        ax.set_ylabel("meadcb/CAMB-1")

        plt.show()

def test_camb(plot=False):
    

    config_camb = {   "hmcode_parameters" : {"file" : os.path.join(KCAP_MODULE_PATH, "one_parameter_hmcode.py"),
                                             "a_0"  : 1.0,
                                             "a_1"  : -0.5}, 
                      "camb"     : {"file" : os.path.join(os.path.join(os.environ["CSL_PATH"], 
                                                                    "boltzmann/pycamb/camb_interface.py")),
                                    "do_reionization"    : False,
                                    "mode"               : "transfer",
                                    "halofit_version"    : "mead",
                                    "nonlinear"          : "pk",
                                    "zmax"               : 4.0,
                                    "nz"                 : 100,
                                    "kmax"               : 20.0,
                                    # "DoLateRadTruncation" : False,
                                    # "AccuracyBoost"       : 2.0
                                    }}

    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "omega_k" : 0.0,
                                               "w"       : -0.9,
                                               "mnu"     : 0.06},
                  "halo_model_parameters"   : {"A"       : 2.0,}}
    # Get cosmosis results
    camb_pipeline = create_pipeline(config_camb)
    block_camb = dict_to_datablock(param_dict)
    camb_pipeline(block_camb)


    # Get CAMB results
    import camb
    # Setting WantCls matters. If not matched, causes 2e-3 differences in the 
    # linear power spectrum.
    pars = camb.CAMBparams(WantCls=False)
    pars.set_cosmology(H0=param_dict["cosmological_parameters"]["h0"]*100, 
                       ombh2=param_dict["cosmological_parameters"]["ombh2"], 
                       omch2=param_dict["cosmological_parameters"]["omch2"],
                       omk=param_dict["cosmological_parameters"]["omega_k"],
                       mnu=param_dict["cosmological_parameters"]["mnu"],
                        )
    pars.InitPower.set_params(As=param_dict["cosmological_parameters"]["A_s"], 
                              ns=param_dict["cosmological_parameters"]["n_s"])

    pars.set_matter_power(redshifts=np.linspace(0, 4, 100)[::-1], kmax=20.0)
    pars.set_dark_energy(w=param_dict["cosmological_parameters"]["w"])

    pars.NonLinearModel.set_params(halofit_version='mead', 
                                   HMCode_A_baryon=param_dict["halo_model_parameters"]["A"], 
                                   HMCode_eta_baryon=1.0 - 0.5*param_dict["halo_model_parameters"]["A"])

    results = camb.get_results(pars)

    kh_camb, z_camb, pk_camb_mead = results.get_nonlinear_matter_power_spectrum()
    kh_lin_camb, z_lin_camb, pk_lin_camb = results.get_linear_matter_power_spectrum()

    p_k_nonlin_camb = pofk_interpolator(pk_camb_mead, kh_camb, z_camb)
    p_k_lin_camb = pofk_interpolator(pk_lin_camb, kh_lin_camb, z_lin_camb)

    k = block_camb["matter_power_nl", "k_h"]
    z = block_camb["matter_power_nl", "z"]
    p_k_nonlin = block_camb["matter_power_nl", "p_k"]

    k_lin = block_camb["matter_power_lin", "k_h"]
    z_lin = block_camb["matter_power_lin", "z"]
    p_k_lin = block_camb["matter_power_lin", "p_k"]

    assert(np.allclose(p_k_lin, p_k_lin_camb(k_lin, z_lin)))
    assert(np.allclose(p_k_nonlin, p_k_nonlin_camb(k, z)))

    if plot:
        import matplotlib
        import matplotlib.pyplot as plt

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

        _ = [ax[1].semilogx(k, p_k_nonlin[i]/p_k_nonlin_camb(k ,z[i]) - 1, c=cmap(i/len(z))) for i in range(len(z))]
        _ = [ax[1].semilogx(k, p_k_nonlin[i]/p_k_nonlin_camb(k ,z[i]) - 1, ls="--", c=cmap(i/len(z))) for i in range(len(z))]

        ax[1].set_title("Non-linear power spectrum")
        ax[1].set_xlabel("k [h/Mpc]")
        ax[1].set_ylabel("CosmoSIS/CAMB-1")

        plt.show()


if __name__ == "__main__":
    test_meadcb(plot=False)
    test_camb(plot=False)