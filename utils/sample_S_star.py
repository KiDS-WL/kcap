import numpy as np

from cosmosis.datablock import option_section

def setup(options):
    alpha = options.get_double(option_section, "alpha", default=0.5)
    use_ln_As = options.get_bool(option_section, "use_ln_A_s", default=False)
    return alpha, use_ln_As

def execute(block, config):
    alpha, use_ln_As = config
    if not use_ln_As:
        S_star = block["cosmological_parameters", "S_star"]
        Omega_m = block["cosmological_parameters", "omega_m"]
        A_s = S_star*2e-9/(Omega_m/0.3)**alpha
    else:
        S_ln_star = block["cosmological_parameters", "S_ln_star"]
        Omega_m = block["cosmological_parameters", "omega_m"]
        ln_A_s_1e10 = S_ln_star/(Omega_m/0.3)**alpha
        A_s = np.exp(ln_A_s_1e10)*1e-10    
    block["cosmological_parameters", "A_s"] = A_s
    return 0

def clean(config):
    pass
