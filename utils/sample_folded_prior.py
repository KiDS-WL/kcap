import numpy as np

from cosmosis.datablock import option_section, names

def setup(options):
    param_name = options[option_section, "name"]
    param_sec, param_name = param_name.split("/")

    fold = options.get_double(option_section, "fold", default=0.0)
    
    return param_sec, param_name, fold

def execute(block, config):
    param_sec, param_name, fold = config

    proxy_param = block[param_sec, param_name + "_folded"]

    block[param_sec, param_name] = fold + np.abs(proxy_param-fold)
    return 0

def cleanup(config):
    pass
