#coding: utf-8
#import cl_to_xi_full
from __future__ import print_function
from builtins import range
import numpy as np
from cosmosis.datablock import option_section, names as section_names
from cl_to_xi import save_xi_00_02, save_xi_22, arcmin_to_radians, SpectrumInterp
from legendre import get_legfactors_00, get_legfactors_02, precomp_GpGm


def setup(options):
    if options.has_value(option_section, "theta"):
        theta = options[option_section, 'theta']
        if np.isscalar(theta):
            theta = np.array([theta])
        theta = arcmin_to_radians(theta)
    else:
        n_theta = options[option_section, "n_theta"]
        theta_min = options[option_section, "theta_min"]
        theta_max = options[option_section, "theta_max"]
        theta_min = arcmin_to_radians(theta_min)
        theta_max = arcmin_to_radians(theta_max)
        theta = np.logspace(np.log10(theta_min), np.log10(theta_max), n_theta)

    corr_type = options.get_int(option_section, 'corr_type')
    ell_max = options.get_int(option_section, "ell_max")

    cl_section = options.get_string(option_section, "input_section_name", "")
    output_section = options.get_string(
        option_section, "output_section_name", "")
    # setup precompute functions and I/O sections
    if corr_type == 0:
        precomp_func = precomp_GpGm
        cl_to_xi_func = save_xi_22
        if not cl_section:
            cl_section = "shear_cl"
        if not output_section:
            output_section = "shear_xi"
    elif corr_type == 1:
        precomp_func = get_legfactors_00
        cl_to_xi_func = save_xi_00_02
        if not cl_section:
            cl_section = "galaxy_cl"
        if not output_section:
            output_section = "galaxy_xi"
    elif corr_type == 2:
        precomp_func = get_legfactors_02
        cl_to_xi_func = save_xi_00_02
        if not cl_section:
            cl_section = "galaxy_shear_cl"
        if not output_section:
            output_section = "galaxy_shear_xi"
    else:
        print("corr_type should be 0 (for spin 2 autocorrelations e.g. xi+/-(theta)),")
        print("1 (for scalar autocorrelations e.g. w(theta) or 2")
        print("for spin 0 x spin 2 correlations e.g. gamma_t(theta)")
        raise ValueError()

    legfacs = precomp_func(np.arange(ell_max + 1), theta)
    return theta, ell_max, legfacs, cl_to_xi_func, cl_section, output_section


def execute(block, config):

    thetas, ell_max, legfacs, cl_to_xi_func, cl_section, output_section = config
    n_theta = len(thetas)

    ell = block[cl_section, "ell"]

    nbina, nbinb = block[cl_section, 'nbin_a'], block[cl_section, 'nbin_b']

    block[output_section, "nbin_a"] = nbina
    block[output_section, "nbin_b"] = nbinb
    block[output_section, "theta"] = thetas
    #block.put_metadata(output_section, "theta", "unit", "radians")

    for i in range(1, nbina + 1):
        for j in range(1, nbinb + 1):
            name = 'bin_%d_%d' % (i, j)
            if block.has_value(cl_section, name):
                c_ell = block[cl_section, name]
            else:
                continue
            cl_interp = SpectrumInterp(ell, c_ell)
            cl_to_xi_func(block, output_section, i, j,
                          cl_interp, thetas, legfacs)
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness.  The joy of python.
    return 0
