from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
cmb = section_names.cmb_cl
matter_powspec = section_names.matter_power_lin


def setup(options):
    alpha = options.get_double(option_section, "alpha", default=0.5)
    return alpha


def execute(block, config):
    alpha = config
    # Get parameters from sampler and CAMB output
    S8_input = block[cosmo, "S8_input"]
    omch2 = block[cosmo, "omch2"]
    ombh2 = block[cosmo, "ombh2"]
    h = block[cosmo, "h0"]
    Omega_m=(omch2+ombh2)/h/h
    # Omega_m=block[cosmo, "omega_m"]
    sigma8_input  = S8_input/(Omega_m/0.3)**alpha
    # print('Omega_m=',Omega_m)
    # print('S8_input=',S8_input)
    # print('sigma8_input=',sigma8_input)
    block[cosmo, "sigma8_input"] = sigma8_input
    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
