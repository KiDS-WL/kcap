import os
import warnings

import numpy as np
import cosmosis.runtime.module

CAMB_INTERFACE = os.path.join(os.path.dirname(__file__), "camb_interface.py")
CONSISTENCY_INTERFACE = "../../utility/consistency/consistency_interface.py"

def dict_to_datablock(d={}):
    b = cosmosis.datablock.DataBlock()
    for section in d.keys():
        for name, value in d[section].items():
            b[section, name] = value

    return b

def test_camb_massless_neutrinos():
    import camb

    Neff = 2.5

    p = camb.CAMBparams()
    p.set_cosmology(ombh2=0.022, omch2=0.11, H0=70.0, nnu=Neff, 
                    num_massive_neutrinos=0,
                    mnu=0.0)

    # Fails with CAMB 1.0.9 but works >1.0.10
    assert np.isclose(Neff, p.N_eff)

def test_setup():
    config = {        "camb"     : {"file"               : CAMB_INTERFACE,
                                    "mode"               : "background",
                                    "zmin"               : 0.0,
                                    "zmid"               : 1.0,
                                    "zmax"               : 1.0,}}

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=config["camb"]["file"])
    np.testing.assert_raises(ValueError, camb_module.setup, dict_to_datablock(config))

def test_parameters():
    config = {        "camb"     : {"file"               : CAMB_INTERFACE,
                                    "do_reionization"    : False,
                                    "mode"               : "background",}}

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=config["camb"]["file"])
    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "w"       : -1.0,
                                               "mnu"     : 0.0,
                                               "num_massive_neutrinos" : 0,
                                               "n_eff"   : 3.046}}

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=config["camb"]["file"])
    camb_module.setup(dict_to_datablock(config))

    block_camb = dict_to_datablock(param_dict)
    camb_module.execute(block_camb)


def test_neutrinos():
    config_camb = {   "camb"     : {"file"               : CAMB_INTERFACE,
                                    "do_reionization"    : False,
                                    "mode"               : "transfer",
                                    "halofit_version"    : "mead",
                                    "nonlinear"          : "pk",
                                    "kmax"               : 20.0,
                                    "zmax_background"    : 1000.0}}

    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "w"       : -1.0,
                                               "mnu"     : 0.0,
                                               "num_massive_neutrinos" : 0,
                                               "n_eff"   : 3.046}}

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=config_camb["camb"]["file"])
    camb_module.setup(dict_to_datablock(config_camb))

    block_camb = dict_to_datablock(param_dict)
    camb_module.execute(block_camb)

    # for k in block_camb.keys('cosmological_parameters'):
    #     print(k[1], block_camb[k[0], k[1]])

    assert np.isclose(block_camb["cosmological_parameters", "massive_nu"], 0)
    assert np.isclose(block_camb["cosmological_parameters", "omega_nu"], 0)
    assert np.isclose(param_dict["cosmological_parameters"]["n_eff"], block_camb["cosmological_parameters", "n_eff"])


    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "w"       : -1.0,
                                               "mnu"     : 0.06,
                                               "massive_nu" : 2,
                                               "n_eff"   : 3.046}}

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=config_camb["camb"]["file"])
    camb_module.setup(dict_to_datablock(config_camb))

    block_camb = dict_to_datablock(param_dict)
    camb_module.execute(block_camb)

    # for k in block_camb.keys('cosmological_parameters'):
    #     print(k[1], block_camb[k[0], k[1]])

    assert np.isclose(block_camb["cosmological_parameters", "num_nu_massive"], 2)
    assert block_camb["cosmological_parameters", "omega_nu"] > 0


def test_consistency():
    config = {     "consistency" : {"file"               : CONSISTENCY_INTERFACE,},
        
                      "camb"     : {"file"               : CAMB_INTERFACE,
                                    "do_reionization"    : False,
                                    "mode"               : "background",
                                    #"halofit_version"    : "mead",
                                    #"nonlinear"          : "pk",
                                    #"kmax"               : 20.0,
                                    "zmax_background"    : 1000.0}}

    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "w"       : -1.0,
                                               "mnu"     : 0.0,
                                               "omega_k" : 0.1,
                                               "massless_nu" : 2,
                                               "num_massive_neutrinos" : 0,
                                               "n_eff"   : 3.046}}

    consistency_module = cosmosis.runtime.module.Module(module_name="consistency", 
                                                 file_path=config["consistency"]["file"])
    consistency_module.setup(dict_to_datablock(config))

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=config["camb"]["file"])
    camb_module.setup(dict_to_datablock(config))

    block = dict_to_datablock(param_dict)
    consistency_module.execute(block)
    camb_module.execute(block)

    block_no_consistency = dict_to_datablock(param_dict)
    camb_module.execute(block_no_consistency)

    for k in block_no_consistency.keys('cosmological_parameters'):
        print(k[1], block_no_consistency[k[0], k[1]])

if __name__ == "__main__":
    test_setup()
    test_parameters()
    test_neutrinos()
    test_camb_massless_neutrinos()
    test_consistency()