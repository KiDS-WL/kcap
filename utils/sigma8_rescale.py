from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
cmb = section_names.cmb_cl
matter_powspec_lin = section_names.matter_power_lin
# matter_powspec_nl = section_names.matter_power_nl


def setup(options):
    return 0


def execute(block, config):

    # Get parameters from sampler and CAMB output
    sigma8_input = block[cosmo, 'sigma8_input']
    sigma8_camb = block[cosmo, 'sigma_8']

    A_s = block[cosmo, 'A_s']
    P_k_lin = block[matter_powspec_lin, 'P_k']


    zmin = block[matter_powspec_lin, 'z'].min()
    if zmin != 0.0:
        raise ValueError(
            "You need to set zmin=0 in CAMB to use the sigma8_rescale module.")

    # Calculate rescale factor
    r = (sigma8_input**2) / (sigma8_camb**2)

    # Rescale CMB Cl and matter power outputs
    A_s *= r
    P_k_lin *= r

    # Save back to block
    block[cosmo, 'A_s'] = A_s
    block[matter_powspec_lin, 'P_k'] = P_k_lin


    block[cosmo, 'sigma_8'] = sigma8_input

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
