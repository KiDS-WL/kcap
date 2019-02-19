from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    try:
        parameter_names = options.get_string(option_section, "parameters")
        parameter_names = [p.strip().lower() for p in parameter_names.split(" ")]
        config = parameter_names
    except errors.BlockNameNotFound:
        parameter_names = []

    config = {}
    config["parameters"] = parameter_names
    for p in parameter_names:
        if p not in ["s8", "s_star"]:
            raise ValueError(f"Derived parameter {p} not supported.")
        if p == "s8" or p == "s_star":
            config["alpha"] = options.get_double(option_section, "alpha", default=0.5)

    return config

def execute(block, config):
    parameter_names = config
    for parameter_name in config["parameters"]:
        if parameter_name == "s8":
            alpha = config["alpha"]
            Omega_m = block["cosmological_parameters", "omega_m"]
            sigma_8 = block["cosmological_parameters", "sigma_8"]
            S8 = sigma_8*(Omega_m/0.3)**alpha
            block["cosmological_parameters", "S8"] = S8
        
        if parameter_name == "s_star":
            alpha = config["alpha"]
            Omega_m = block["cosmological_parameters", "omega_m"]
            A_s = block["cosmological_parameters", "A_s"]
            S_star = A_s/2e-9*(Omega_m/0.3)**alpha
            block["cosmological_parameters", "S_star"] = S_star

    return 0

def cleanup(config):
    pass
