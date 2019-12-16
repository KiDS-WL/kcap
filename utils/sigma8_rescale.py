from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
cmb = section_names.cmb_cl
matter_powspec_lin = section_names.matter_power_lin
growth_params = section_names.growth_parameters


def setup(options):
    print(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ")
    print(" sigma8_rescale only deals with linear power spectra. Nonlinear power is not changed! ")
    print(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    return 0




def execute(block, config):

    # Get parameters from sampler and CAMB output
    sigma8_input = block[cosmo, 'sigma8_input']
    sigma8_camb = block[cosmo, 'sigma_8']

    zmin = block[matter_powspec_lin, 'z'].min()
    if zmin != 0.0:
        raise ValueError(
            "You need to set zmin=0 in CAMB to use the sigma8_rescale module.")

    # Calculate rescale factor
    r = (sigma8_input**2) / (sigma8_camb**2)

    # rescale parameters first
    A_s = block[cosmo, 'A_s']
    A_s *= r
    block[cosmo, 'A_s'] = A_s
    block[cosmo, 'sigma_8'] = sigma8_input

    if(block.has_value(growth_params, "sigma_8")):
        block[growth_params, "sigma_8"] = sigma8_input

    if(block.has_value(growth_params, "fsigma_8")):
        fsigma8 = block[growth_params, "fsigma_8"]
        fsigma8 *= sigma8_input/sigma8_camb
        block[growth_params, "fsigma_8"] = fsigma8

    # Check if power spectra exists
    # Then Rescale CMB Cl and matter power outputs
    if(block.has_value(cmb, 'TT')):
        TT = block[cmb, 'TT']
        TT *= r
        block[cmb, 'TT'] = TT

    if(block.has_value(cmb, 'EE')):
        EE = block[cmb, 'EE']
        EE *= r
        block[cmb, 'EE'] = EE

    if(block.has_value(cmb, 'BB')):
        BB = block[cmb, 'BB']
        BB *= r
        block[cmb, 'BB'] = BB

    if(block.has_value(cmb, 'TE')):
        TE = block[cmb, 'TE']
        TE *= r
        block[cmb, 'TE'] = TE

    if(block.has_value(cmb, 'PP')):
        PP = block[cmb, 'PP']
        PP *= r
        block[cmb, 'PP'] = PP

    if(block.has_value(cmb, 'PT')):
        PT = block[cmb, 'PT']
        PT *= r
        block[cmb, 'PT'] = PT

    if(block.has_value(cmb, 'PE')):
        PE = block[cmb, 'PE']
        PE *= r
        block[cmb, 'PE'] = PE

    if(block.has_value(cmb, 'P_k')):
        P_k_lin = block[matter_powspec_lin, 'P_k']
        P_k_lin *= r
        block[matter_powspec_lin, 'P_k'] = P_k_lin

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
