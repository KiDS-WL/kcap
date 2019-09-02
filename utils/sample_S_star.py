import numpy as np

from cosmosis.datablock import option_section

def setup(options):
    alpha = options.get_double(option_section, "alpha", default=0.5)
    use_ln_As = options.get_bool(option_section, "use_ln_A_s", default=True)
    return alpha, use_ln_As

def execute(block, config):
    alpha, use_ln_As = config
    if not use_ln_As:
        S_star = block["cosmological_parameters", "S_star"]
        omch2 = block["cosmological_parameters", "omch2"]
        A_s = S_star*3e-9/(omch2/0.1)**alpha
    else:
        S_star = block["cosmological_parameters", "S_star"]
        omch2 = block["cosmological_parameters", "omch2"]
        ln_As = S_star*3.0/(omch2/0.1)**alpha
        A_s = np.exp(ln_As)*1e-10    
    block["cosmological_parameters", "A_s"] = A_s
    return 0

def clean(config):
    pass
