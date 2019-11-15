import os

import numpy as np

import cosmosis.runtime.module

import sys
sys.path.append("../kcap")
from cosmosis_utils import dict_to_datablock

KCAP_MODULE_PATH = "../utils"

CSL_PATH = os.environ["CSL_PATH"] if "CSL_PATH" in os.environ else "../cosmosis-standard-library"

def test_camb_background_module():
    config_background = {"camb_background" : {"w_name" : "w_prime",
                                              "feedback" : 4,
                                              "zmin_background" : 0.0,
                                              "zmax_background" : 10.0,
                                              "nz_background"   : 1000,}}
    config_camb = {      "camb" : {"do_reionization" : False,
                                   "zmin_background" : 0.0,
                                   "zmax_background" : 10.0,
                                   "nz_background"   : 1000,}}

    config_camb_old = { "camb_old" : {"mode" : "background",
                                       "feedback" : 0,
                                       "accuracy_boost" : 2.0,
                                       "background_zmin" : 0.0,
                                       "background_zmax" : 10.0,
                                       "background_nz"   : 1000,}}

    param_dict = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "omega_k" : 0.0,
                                               "w"       : -0.9,
                                               "w_prime" : -0.9,
                                               "mnu"     : 0.06}}
    
    camb_background_module = cosmosis.runtime.module.Module(module_name="camb_background", 
                                                            file_path=os.path.join(KCAP_MODULE_PATH, "camb_background.py"))
    camb_background_module.setup(dict_to_datablock(config_background))

    camb_module = cosmosis.runtime.module.Module(module_name="camb", 
                                                 file_path=os.path.join(CSL_PATH, 
                                                                        "boltzmann/pycamb/camb_interface.py"))
    camb_module.setup(dict_to_datablock(config_camb))
    
    camb_old_module = cosmosis.runtime.module.Module(module_name="camb_old", 
                                                     file_path=os.path.join(CSL_PATH, 
                                                                            "boltzmann/camb/camb.so"))
    camb_old_module.setup(dict_to_datablock(config_camb_old))

    consistency_module = cosmosis.runtime.module.Module(module_name="consistency", 
                                                        file_path=os.path.join(CSL_PATH, 
                                                                               "utility/consistency/consistency_interface.py"))
    consistency_module.setup(dict_to_datablock({}))

    block_background = dict_to_datablock(param_dict)
    camb_background_module.execute(block_background)
    
    block_camb = dict_to_datablock(param_dict)
    camb_module.execute(block_camb)

    assert np.allclose(block_background["distances", "z"], block_camb["distances", "z"])
    assert np.allclose(block_background["distances", "D_A"], block_camb["distances", "D_A"])

    param_dict_nu0 = {"cosmological_parameters" : {"omch2"   : 0.1,
                                               "ombh2"   : 0.022,
                                               "h0"      : 0.7,
                                               "n_s"     : 0.96,
                                               "A_s"     : 2.1e-9,
                                               "omega_k" : 0.0,
                                               "w"       : -0.9,
                                               "w_prime" : -0.9,
                                               "tau"     : 0.08,
                                               "mnu"     : 0.0,
                                               "num_massive_neutrinos" : 0,
                                               "massive_nu" : 0,
                                               "omega_nu"  : 0.0}}

    block_camb_old = dict_to_datablock(param_dict_nu0)
    consistency_module.execute(block_camb_old)
    camb_old_module.execute(block_camb_old)

    block_background_nu0 = dict_to_datablock(param_dict_nu0)
    camb_background_module.execute(block_background_nu0)

    print("WARNING: Angular diameter distance tolerance is set to 5e-5. The agreement should be much better than this!")
    assert np.allclose(block_background_nu0["distances", "z"], block_camb_old["distances", "z"])
    # print(block_camb_old["distances", "D_A"]/block_background_nu0["distances", "D_A"] - 1)
    assert np.allclose(block_background_nu0["distances", "D_A"], block_camb_old["distances", "D_A"], rtol=5e-5)    



if __name__ == "__main__":
    test_camb_background_module()