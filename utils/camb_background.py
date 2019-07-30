import camb

import numpy as np
from cosmosis.datablock import option_section, names

def get_optional_params(block, section, names):
    params = {}    
    for name in names:
        cosmosis_name = name
        output_name = name
        if isinstance(name, (list, tuple)):
            if len(name) == 2 and isinstance(name[1], str):
                # Output name specified
                output_name = name[1]
                cosmosis_name = name[0]
        if block.has_value(section, cosmosis_name):
            params[output_name] = block[section, cosmosis_name]
    return params

def setup(options):
    config = {}

    config["w_name"] = options.get_string(option_section, "w_name", default="w")
    config["wa_name"] = options.get_string(option_section, "wa_name", default="wa")

    config["cosmology_settings"] = get_optional_params(options, option_section, ["neutrino_hierarchy" ,"theta_H0_range"])

    config['zmin_background'] = options.get_double(option_section, 'zmin_background', default=0.0)
    config['zmax_background'] = options.get_double(option_section, 'zmax_background', default=10.0)
    config['nz_background'] = options.get_int(option_section, 'nz_background', default=1000)

    return config

def execute(block, config):
    # Get optional parameters from datablock.
    cosmology_params = get_optional_params(block, names.cosmological_parameters, ["TCMB", "mnu", "nnu", "standard_neutrino_neff",
                                                          "num_massive_neutrinos",])

    dark_energy_params = get_optional_params(block, names.cosmological_parameters, [(config["w_name"], "w"), 
                                                            (config["wa_name"], "wa")])

    # Set h if provided, otherwise look for theta_mc
    if block.has_value(names.cosmological_parameters, "hubble"):
        cosmology_params["H0"] = block[names.cosmological_parameters, "hubble"]
    elif block.has_value(names.cosmological_parameters, "h0"):
        cosmology_params["H0"] = block[names.cosmological_parameters, "h0"]*100
    else:
        cosmology_params["cosmomc_theta"] = block[names.cosmological_parameters, "cosmomc_theta"]/100

    p = camb.CAMBparams()
    p.set_cosmology(ombh2 = block[names.cosmological_parameters, 'ombh2'],
                    omch2 = block[names.cosmological_parameters, 'omch2'],
                    omk = block[names.cosmological_parameters, 'omega_k'],
                    **config["cosmology_settings"],
                    **cosmology_params)

    p.set_dark_energy(**dark_energy_params)
    
    r = camb.get_background(p, no_thermo=True)

    z_background = np.linspace(config["zmin_background"], config["zmax_background"], config["nz_background"])

    block[names.distances, "z"] = z_background
    block[names.distances, "a"] = 1/(z_background+1)
    block[names.distances, "D_A"] = r.angular_diameter_distance(z_background)
    block[names.distances, "D_C"] = r.comoving_radial_distance(z_background)
    block[names.distances, "D_M"] = r.comoving_radial_distance(z_background)
    d_L = r.luminosity_distance(z_background)
    block[names.distances, "D_L"] = d_L
    block[names.distances, "MU"] = np.insert(5*np.log10(d_L[1:])+25, 0, np.nan)

    block[names.distances, "H"] = r.h_of_z(z_background)

    return 0

def cleanup(config):
    pass
