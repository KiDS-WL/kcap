import os

import numpy as np
import scipy.interpolate
import scipy.integrate

from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors
import tpst

pi = np.pi


def log_interpolate(x_arr, y_arr, x):
    log_x_arr = np.log(x_arr)
    if np.any(y_arr <= 0):
        intp = scipy.interpolate.InterpolatedUnivariateSpline(
                    log_x_arr, y_arr, ext=2)
    else:
        log_intp = scipy.interpolate.InterpolatedUnivariateSpline(
                        log_x_arr, np.log(y_arr), ext=2)

        def intp(log_x):
            return np.exp(log_intp(log_x))

    def low_extrap(log_x):
        dx = log_x_arr[1] - log_x_arr[0]
        dy = np.log(y_arr[1]) - np.log(y_arr[0])
        return y_arr[0] * np.exp(dy/dx * (log_x - log_x_arr[0]))

    def high_extrap(log_x):
        dx = log_x_arr[-1] - log_x_arr[-2]
        dy = np.log(y_arr[-1]) - np.log(y_arr[-2])
        return y_arr[-1] * np.exp(dy/dx * (log_x - log_x_arr[-1]))

    log_x = np.log(x)
    y = np.piecewise(log_x,
                     [log_x < log_x_arr[0],
                      (log_x >= log_x_arr[0]) & (log_x <= log_x_arr[-1])],
                     [low_extrap, intp, high_extrap])
    return y


def create_bandpower_window_function(n_ell_bin=8,
                                     ell_bin_min=100,
                                     ell_bin_max=1500,
                                     ell_window_min=0.1,
                                     ell_window_max=10000,
                                     n_ell_window=1000,
                                     theta_range=(0.5, 300.0),
                                     apodise_theta=0.5):
    ell_bin_edges = np.geomspace(ell_bin_min, ell_bin_max, n_ell_bin+1)

    N_BP = np.log(ell_bin_edges[1:]/ell_bin_edges[:-1])
    ell = np.geomspace(ell_window_min, ell_window_max, n_ell_window)

    theta_range_bandpower = (theta_range[0]/60/180*pi,
                             theta_range[1]/60/180*pi)
    theta_range_apo_bandpower = \
        (theta_range[0]/60/180*pi*np.exp(-apodise_theta/2),
         theta_range[1]/60/180*pi*np.exp(apodise_theta/2))

    W_EE = np.zeros((n_ell_bin, ell.size))

    for i, (l, u) in enumerate(zip(ell_bin_edges[:-1], ell_bin_edges[1:])):
        print(f"Bin {i+1}: {l:.1f}-{u:.1f}")
        try:
            import tqdm
            iter = tqdm.tqdm(range(ell.size))
        except ImportError:
            iter = range(ell.size)
        for j in iter:
            # if ell[j] < 2:
            #     continue
            W_EE[i, j] = tpst.flat_sky.compute_weights(
                                    ell=ell[j], ell_bin=(l, u),
                                    theta_range=theta_range_bandpower,
                                    apodise_theta=(apodise_theta,
                                                   *theta_range_apo_bandpower),
                                    estimator="EE",
                                    integrator="quad")

    return ell, ell/(2*N_BP[:, None])*W_EE, ell_bin_edges


def create_cosebis_window_function(n_mode=5,
                                   ell_window_min=1,
                                   ell_window_max=50000,
                                   n_ell_window=50000,
                                   theta_range=(0.5, 300.0),):
    ell = np.geomspace(ell_window_min, ell_window_max, n_ell_window)
    theta_range_cosebis = (theta_range[0]/60/180*pi, theta_range[1]/60/180*pi)

    cosebis_modes = np.array(list(range(1, n_mode+1)))
    W_log_cosebis = np.zeros((len(cosebis_modes), ell.size))

    print("Computing roots")
    roots, normalisations = tpst.cosebis.roots_normalization(
                                n_max=max(cosebis_modes),
                                z_max=np.log(theta_range_cosebis[1]
                                             / theta_range_cosebis[0]))

    for i, n in enumerate(cosebis_modes):
        print(f"Mode {n}")
        try:
            import tqdm
            iter = tqdm.tqdm(range(ell.size))
        except ImportError:
            iter = range(ell.size)
        for j in iter:
            # if ell[j] < 1:
            #     continue
            W_log_cosebis[i, j] = tpst.flat_sky.compute_weights(
                                    ell=ell[j], ell_bin=n,
                                    theta_range=theta_range_cosebis,
                                    estimator="log-cosebis",
                                    log_cosebis_roots=roots,
                                    log_cosebis_normalisations=normalisations)

    return ell, ell*W_log_cosebis/(2*pi), cosebis_modes


def setup(options):
    extra_info = {}

    do_extrapolate_low = options.get_bool(
        option_section, "do_extrapolate_low", True)
    do_extrapolate_high = options.get_bool(
        option_section, "do_extrapolate_high", True)

    extra_info["do_extrapolate_low"] = do_extrapolate_low
    if do_extrapolate_low:
        extra_info["ell_min_extrapolate"] = options.get_double(
                                                option_section,
                                                "ell_min_extrapolate", -1.0)
        if extra_info["ell_min_extrapolate"] < 0:
            extra_info["ell_min_extrapolate"] = None

    extra_info["do_extrapolate_high"] = do_extrapolate_high
    if do_extrapolate_high:
        extra_info["ell_max_extrapolate"] = options.get_double(
                                                option_section,
                                                "ell_max_extrapolate", -1.0)
        if extra_info["ell_max_extrapolate"] < 0:
            extra_info["ell_max_extrapolate"] = None

    extra_info["integration_mode"] = options.get_string(option_section,
                                                        "integration_mode",
                                                        "simps")

    window_function_file = options.get_string(option_section,
                                              "window_function_file")
    if os.path.isfile(window_function_file):
        d = np.load(window_function_file)
        ell = d["ell"]
        window_function = d["EE"]
        stat_name = d["statistic"]
        try:
            requested_stat = options[option_section, "statistic"]
            if requested_stat != stat_name:
                raise ValueError(f"Requested statistic {requested_stat} does not "
                                 f"agree with statistic of the provided "
                                 f"window function {stat_name}.")
        except errors.BlockError:
            pass
        if stat_name == "bandpower":
            ell_bin_edges = d["bin_edges"]
        elif stat_name == "cosebis":
            cosebis_modes = d["modes"]
    else:
        print(f"Window function file {window_function_file} does not exist.")

        stat_name = options[option_section, "statistic"]
        if stat_name == "bandpower":
            ell, window_function, ell_bin_edges = \
                create_bandpower_window_function(
                                     n_ell_bin=options[option_section, "n_ell"],
                                     ell_bin_min=options[option_section, "ell_min"],
                                     ell_bin_max=options[option_section, "ell_max"],
                                     n_ell_window=3000,
                                     ell_window_max=30000)
            os.makedirs(os.path.dirname(window_function_file), exist_ok=True)
            np.savez(window_function_file, ell=ell, EE=window_function,
                     bin_edges=ell_bin_edges, statistic="bandpower", n_ell_bin=options[option_section, "n_ell"],
                                     ell_bin_min=options[option_section, "ell_min"],
                                     ell_bin_max=options[option_section, "ell_max"])
        elif stat_name == "cosebis":
            ell, window_function, cosebis_modes = \
                create_cosebis_window_function(n_mode=options[option_section, "n_mode"],
                                               n_ell_window=50000)
            os.makedirs(os.path.dirname(window_function_file), exist_ok=True)
            np.savez(window_function_file, ell=ell, EE=window_function,
                     modes=cosebis_modes, statistic="cosebis", n_mode=options[option_section, "n_mode"])

    print(f"Window function for {stat_name}: "
          f"ell_min={ell[0]}, ell_max={ell[-1]}, n_ell={len(ell)}")
    if stat_name == "bandpower":
        extra_info["l_min_vec"] = ell_bin_edges[:-1]
        extra_info["l_max_vec"] = ell_bin_edges[1:]
        extra_info["ell"] = 10**(0.5 * (np.log10(ell_bin_edges[:-1]) + np.log10(ell_bin_edges[1:])))
    elif stat_name == "cosebis":
        extra_info["cosebis_n"] = cosebis_modes

    input_section_name = options.get_string(option_section,
                                            "input_section_name",
                                            "shear_cl")
    output_section_name = options.get_string(option_section,
                                             "output_section_name",
                                             "")

    if output_section_name == "":
        if stat_name == "bandpower":
            output_section_name = "bandpower_shear_e"
        elif stat_name == "cosebis":
            output_section_name = "cosebis"

    return (ell, window_function,
            input_section_name, output_section_name, extra_info)


def execute(block, config):
    (ell, window_function,
     input_section_name, output_section_name, extra_info) = config

    ell_data_block = block[input_section_name, "ell"]
    n_z_bin_a = block[input_section_name, "nbin_a"]
    n_z_bin_b = block[input_section_name, "nbin_b"]
    bin_indices = []
    for i in range(n_z_bin_a):
        jmax = i+1 if block.get_bool(input_section_name, 'is_auto') else n_z_bin_b
        for j in range(jmax):
            if block.get_bool(input_section_name, 'auto_only'):
                if j!=i:
                    continue
            bin_indices.append((i, j))

    m = np.ones(ell.size, dtype=bool)
    if extra_info["do_extrapolate_low"]:
        if extra_info["ell_min_extrapolate"] is not None:
            m &= (ell >= extra_info["ell_min_extrapolate"])
    else:
        m &= (ell >= ell_data_block.min())
    if extra_info["do_extrapolate_high"]:
        if extra_info["ell_max_extrapolate"] is not None:
            m &= (ell <= extra_info["ell_max_extrapolate"])
    else:
        m &= (ell <= ell_data_block.max())

    for idx_1, idx_2 in bin_indices:
        bin_tag = f"bin_{idx_1+1}_{idx_2+1}"
        Cl_data_block = block[input_section_name, bin_tag]
        Cl_EE = np.zeros(len(ell))
        Cl_EE[m] = log_interpolate(ell_data_block, Cl_data_block, ell[m])

        if extra_info["integration_mode"] == "simps":
            EE = scipy.integrate.simps(window_function * Cl_EE, ell, axis=1)
        elif extra_info["integration_mode"] == "quad":
            EE = []
            for i in range(window_function.shape[0]):
                intp = scipy.interpolate.InterpolatedUnivariateSpline(
                    np.log(ell), ell*window_function[i] * Cl_EE)
                EE.append(scipy.integrate.quad(
                                intp,
                                a=np.log(ell[0]), b=np.log(ell[-1]),
                                limit=2000)[0])
            EE = np.array(EE)

        block[output_section_name, bin_tag] = EE

    for name, value in extra_info.items():
        if value is None:
            continue
        block[output_section_name, name] = value
        
    block[output_section_name, 'auto_only'] = block.get_bool(input_section_name, 'auto_only')
    block[output_section_name, 'is_auto'] = block.get_bool(input_section_name, 'is_auto')
    if block.has_value(input_section_name, 'nbin'):
        block[output_section_name, 'nbin'] = block[input_section_name, 'nbin']
    block[output_section_name, 'nbin_a'] = block[input_section_name, 'nbin_a']
    block[output_section_name, 'nbin_b'] = block[input_section_name, 'nbin_b']
    block[output_section_name, 'sample_a'] = block[input_section_name, 'sample_a']
    block[output_section_name, 'sample_b'] = block[input_section_name, 'sample_b']
    block[output_section_name, 'sep_name'] = block[input_section_name, 'sep_name']
    block[output_section_name, 'save_name'] = block[input_section_name, 'save_name']

    return 0


def cleanup(block):
    pass
