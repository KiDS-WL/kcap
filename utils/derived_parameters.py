from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    try:
        parameter_names = options.get_string(option_section, "parameters")
        parameter_names = [p.strip() for p in parameter_names.split(" ")]
        config = parameter_names
    except errors.BlockNameNotFound:
        parameter_names = []

    for p in parameter_names:
        if p not in ["S8"]:
            raise ValueError(f"Derived parameter {p} not supported.")

    config = parameter_names
    return config

def execute(block, config):
    parameter_names = config
    for parameter_name in parameter_names:
        if parameter_name == "S8":
            Omega_m = block["cosmological_parameters", "omega_m"]
            sigma_8 = block["cosmological_parameters", "sigma_8"]
            S8 = sigma_8*(Omega_m/0.3)**0.5

            block["cosmological_parameters", "S8"] = S8

    return 0

def cleanup(config):
    pass
