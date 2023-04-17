# Copyright (c) 2021 Tilman Troester

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.




# # Example setup:
#
# [bin_xipm]
# file = %(REANALYSIS_PATH)s/bin_xi_pm.py
# twopoint_file = %(2PT_FILE)s
#
# # 9 log-spaced bins between 0.5 and 300.0
# angular_bin_edges = 0.5 1.01777898 2.0717481 4.21716333 8.58428037 17.4738002 35.56893304 72.40262468 147.37973879 300.0
#
# # How are the datablock bins counted: from 1 or 0?
# bin_offset = 1
#
# # Use measured n_pair to weigh the correlation function? The curly brackets
# # {bin_a} and {bin_b} will get replaced with the bin numbers.
# n_pair_file_template = %(KIDS1000_CAT_TO_OBS_DIR)s/data/kids/xipm/XI_K1000_ALL_BLIND_C_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_nBins_5_Bin{bin_b}_Bin{bin_a}.ascii
#
# # Specify the columns to use. These are the default values
# n_pair_theta_col = 0
# n_pair_weight_col = 10
#
# # Use a hack to reproduce the exact binning of Asgari et al. 2021?
# asgari2019_compatability_mode = F
#
# # Specify the name of the input/output sections.
# xi_plus_output_section = shear_xi_plus
# xi_minus_output_section = shear_xi_minus

import numpy as np
import scipy.integrate
import scipy.interpolate

import astropy.io.fits

from cosmosis.datablock import option_section


def setup(options):
    twopoint_file = options[option_section, "twopoint_file"]
    # Load data from the twopoint file
    with astropy.io.fits.open(twopoint_file) as hdu:
        n_z_bin = hdu["xim"].header["N_ZBIN_1"]
        n_angular_bin_2pt = hdu["xim"].header["N_ANG"]
        theta_binned_2pt = hdu["xim"].data["ANG"][:n_angular_bin_2pt]
        ang_index = hdu["xim"].data.names.index("ANG")
        angle_unit = hdu["xim"].header.get('TUNIT{}'.format(ang_index + 1))
        if angle_unit.lower() == "arcmin":
            theta_binned_2pt *= np.pi/180/60

    angular_bin_edges = options[option_section, "angular_bin_edges"]
    angular_bin_edges *= np.pi/180/60
    n_angular_bin = len(angular_bin_edges) - 1

    if n_angular_bin != n_angular_bin_2pt:
        raise ValueError("Number of bins in 2pt file and number of bin edges "
                         "do not match.")

    # Bin counting from bin_0_0 or bin_1_1?
    bin_offset = options.get_int(option_section, "bin_offset", 1)

    # Use n_pair weighting? If so, the file name template should be of the form
    # ...{bin_a}...{bin_b}...
    n_pair_file_template = options.get_string(option_section,
                                              "n_pair_file_template", "")
    if n_pair_file_template != "":
        asgari2021_compatability_mode = options.get_bool(
                                            option_section,
                                            "asgari2019_compatability_mode",
                                            False)
        theta_col = options.get_int(option_section, "n_pair_theta_col", 0)
        weight_col = options.get_int(option_section, "n_pair_weight_col", 10)

        bins = [(i+bin_offset, j+bin_offset) for i in range(n_z_bin)
                for j in range(i+1)]
        if asgari2021_compatability_mode:
            bins_asgari2021_ordering = [(i+bin_offset, j+bin_offset)
                                        for i in range(n_z_bin)
                                        for j in range(i, n_z_bin)]
        weights = {}
        for i, (bin_a, bin_b) in enumerate(bins):
            if asgari2021_compatability_mode:
                bin_w_a, bin_w_b = bins_asgari2021_ordering[i]
                n_pair_file = n_pair_file_template.format(bin_a=bin_w_b,
                                                          bin_b=bin_w_a)
            else:
                n_pair_file = n_pair_file_template.format(bin_a=bin_a,
                                                          bin_b=bin_b)
            theta, n_pair = np.loadtxt(n_pair_file,
                                       usecols=[theta_col, weight_col],
                                       unpack=True)
            # Assume theta is in arcmin
            theta = theta/(180*60)*np.pi
            weights[(bin_a, bin_b)] = theta, n_pair

        interpolate = "weights"
    else:
        weights = None
        interpolate = "log"

    constant_c_offset = options.get_bool(option_section, "add_c_term", False)

    xi_plus_input_section = options.get_string(option_section,
                                               "xi_plus_input_section",
                                               "shear_xi_plus")
    xi_minus_input_section = options.get_string(option_section,
                                                "xi_minus_input_section",
                                                "shear_xi_minus")
    xi_plus_output_section = options.get_string(option_section,
                                                "xi_plus_output_section",
                                                "shear_xi_plus")
    xi_minus_output_section = options.get_string(option_section,
                                                 "xi_minus_output_section",
                                                 "shear_xi_minus")

    return (angular_bin_edges, n_angular_bin,
            xi_plus_input_section, xi_plus_output_section,
            xi_minus_input_section, xi_minus_output_section,
            theta_binned_2pt, bin_offset, weights, interpolate,
            constant_c_offset)


def execute(block, config):
    (angular_bin_edges, n_angular_bin,
     xi_plus_input_section, xi_plus_output_section,
     xi_minus_input_section, xi_minus_output_section,
     theta_binned_2pt, bin_offset, weights, interpolate,
     constant_c_offset) = config

    for kind, xi_input_section, xi_output_section in zip(
            ["plus", "minus"],
            [xi_plus_input_section, xi_minus_input_section],
            [xi_plus_output_section, xi_minus_output_section]):
        theta = block[xi_input_section, "theta"]
        n_z_bin = block[xi_input_section, "nbin_a"]
        bins = [(i+bin_offset, j+bin_offset) for i in range(n_z_bin)
                for j in range(i+1)]

        for (bin_a, bin_b) in bins:
            bin_name = f"bin_{bin_a}_{bin_b}"
            xi = block[xi_input_section, bin_name]

            # When using quadrature, interpolate and integrate in log(theta)
            # or theta.
            if interpolate == "log":
                w = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=np.log(theta), y=theta)
                y = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=np.log(theta), y=xi)
            elif interpolate == "lin":
                w = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=theta, y=theta)
                y = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=theta, y=xi)
            elif interpolate == "weights":
                x, w = weights[(bin_a, bin_b)]
                y = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=np.log(theta), y=xi)(np.log(x))

            xi_binned = np.zeros(n_angular_bin)
            for i, (l, u) in enumerate(zip(angular_bin_edges[:-1],
                                           angular_bin_edges[1:])):
                if interpolate == "log":
                    normalisation = scipy.integrate.quad(
                                        lambda log_x: np.exp(log_x)*w(log_x),
                                        a=np.log(l), b=np.log(u))[0]
                    xi_binned[i] = scipy.integrate.quad(
                                        lambda log_x: (np.exp(log_x)
                                                       * w(log_x)*y(log_x)),
                                        a=np.log(l), b=np.log(u))[0]
                    xi_binned[i] /= normalisation
                elif interpolate == "lin":
                    normalisation = scipy.integrate.quad(
                                        lambda x: w(x), a=l, b=u)[0]
                    xi_binned[i] = scipy.integrate.quad(
                                        lambda x: w(x)*y(x),
                                        a=l, b=u)[0]
                    xi_binned[i] /= normalisation
                elif interpolate == "weights":
                    mask = np.zeros(w.size, dtype=bool)
                    # Find closest theta value to bin boundary
                    l_idx = np.argmin(np.abs(l-x))
                    u_idx = np.argmin(np.abs(u-x))
                    mask[l_idx:u_idx+1] = True
                    normalisation = w[mask].sum()
                    xi_binned[i] = (w[mask]*y[mask]).sum()/normalisation
                else:
                    mask = (l <= theta) & (theta < u)
                    x = theta[mask]
                    log_x = np.log(x)
                    w = theta[mask]
                    y = xi[mask]
                    normalisation = np.trapz(x=log_x, y=x*w)
                    xi_binned[i] = np.trapz(x=log_x, y=x*w*y)
                    xi_binned[i] /= normalisation

            # Add constant c term to xi_plus
            if kind == "plus" and constant_c_offset:
                xi_binned += block["shear_c_bias", "delta_c"]
            # Overwrite the xi_pm sections
            block[xi_output_section, bin_name] = xi_binned
            # For KiDS-1000 scale_cuts module compatability
            block[xi_output_section, f"theta_bin_{bin_a}_{bin_b}"] = \
                theta_binned_2pt/np.pi*180*60
        block[xi_output_section, "theta"] = theta_binned_2pt

    return 0


def cleanup(config):
    pass
