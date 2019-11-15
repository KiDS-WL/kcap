import os
import warnings

import numpy as np
import cosmosis.runtime.module

import sys
sys.path.append("./kcap")
from cosmosis_utils import dict_to_datablock

CORRELATED_PRIORS_INTERFACE = "./utils/correlated_priors.py"
PARAMETER_COVARIANCE_FILE = "tests/parameter_covariance.txt"

def test_parameters(plot=True):
    corr = np.array([[1.0, 0.8, 0.5], [0.8, 1.0, 0.8], [0.5, 0.8, 1.0]])
    std = np.diag([2.0, 1.5, 0.8])
    cov = std @ corr @ std
    np.savetxt(PARAMETER_COVARIANCE_FILE, cov)

    config = {"correlated_priors" : {"file"                     : CORRELATED_PRIORS_INTERFACE,
                                     "uncorrelated_parameters"  : "cosmological_parameters/p_1 cosmological_parameters/p_2 other_parameters/p_1",
                                     "output_parameters"        : "cosmological_parameters/ombh2 cosmological_parameters/omch2 cosmological_parameters/a_s",
                                     "covariance"               : PARAMETER_COVARIANCE_FILE,}}

    module = cosmosis.runtime.module.Module(module_name="correlated_priors", 
                                            file_path=config["correlated_priors"]["file"])
    param_dict = {"cosmological_parameters" : {"p_1"     : 0.1,
                                               "p_2"     : -0.3},
                  "other_parameters"        : {"p_1"     : 0.5}}

    module.setup(dict_to_datablock(config))

    block = dict_to_datablock(param_dict)
    module.execute(block)

    p = np.array([param_dict["cosmological_parameters"]["p_1"], param_dict["cosmological_parameters"]["p_2"], param_dict["other_parameters"]["p_1"]])
    p = np.linalg.cholesky(cov) @ p

    p_cosmosis = np.array([block["cosmological_parameters", "ombh2"], block["cosmological_parameters", "omch2"], block["cosmological_parameters", "a_s"]])

    np.testing.assert_almost_equal(p, p_cosmosis)

    if plot:
        import matplotlib.pyplot as plt

        p = np.random.normal(size=(1000, 3))

        config = {"correlated_priors" : {"file"                     : CORRELATED_PRIORS_INTERFACE,
                                         "uncorrelated_parameters"  : "other_parameters/p_1 other_parameters/p_2 other_parameters/p_3",
                                         "output_parameters"        : "cosmological_parameters/ombh2 cosmological_parameters/omch2 cosmological_parameters/a_s",
                                         "covariance"               : PARAMETER_COVARIANCE_FILE,}}

        module = cosmosis.runtime.module.Module(module_name="correlated_priors", 
                                                file_path=config["correlated_priors"]["file"])
        module.setup(dict_to_datablock(config))

        p_out = []
        for p_in in p:
            param_dict = {"other_parameters"        : {"p_1" : p_in[0], "p_2" : p_in[1], "p_3" : p_in[2]}}
            block = dict_to_datablock(param_dict)
            module.execute(block)

            p_out.append(np.array([block["cosmological_parameters", "ombh2"], block["cosmological_parameters", "omch2"], block["cosmological_parameters", "a_s"]]))

        p_out = np.array(p_out)
        plt.scatter(p_out[:,0], p_out[:,1], s=1)
        plt.show()




if __name__ == "__main__":
    test_parameters()

