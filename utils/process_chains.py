import re
import warnings
import os
import configparser

import numpy as np

try:
    import getdist
except ImportError:
    warnings.warn("Could not import getdist. Loading chains is not supported.")

def compute_sample_range_and_coverage(samples, weights, idx):
    """Compute range of values of samples[idx] and their coverage."""
    l, u = samples[idx].min(), samples[idx].max()
    coverage_1d = np.sum(weights[(l <= samples) & (samples <= u)])
    coverage_nd = np.sum(weights[idx])
    return (l,u), coverage_1d, coverage_nd

def find_HPD_CI(samples, weights=None, coverage_1d_threshold=0.68, 
            sort_idx=None, log_posterior=None, verbose=False):
    """Find the highest posterior density credible interval.
    
    Finds the HPD CI of samples, using the posterior masses given my weights. 
    The result is interpolated between the two samples that lie on either side
    of the target coverage. For example, if the covarage of the n samples with
    the highest posterior density is 67% and the covarage of the n+1 samples 
    is 72%, the result is interpolated to the target coverage (e.g., 68%).
    
    Arguments
    ---------
    samples : numpy.array
        Samples of the parameter to compute the CI for.
    weights : numpy.array
        Posterior mass of the samples used to compute the coverage. If not 
        provided is assumed to be unifrom 1/n_sample.
    coverage_1d_threshold : float
        Target 1D coverage of the CI. (Default 0.68)
    sort_idx : numpy.array
        Array of indices that sort the samples in descending posterior value.
        If not provided, this will be computed internally from log_posterior.
    log_posterior : numpy.array
        Array of the log posterior values of the samples.
    verbose : bool
        Print extra outputs.
        
    Returns
    -------
    (l, u) : tuple
        CI of samples.
    coverage_1d : float
        Achieved 1D coverage.
    coverage_nd : float
        N-dimensional coverage of the samples in CI.
    n_sample : int
        Number of samples in CI.
    """
    if weights is None:
        weights = np.ones_like(samples)/samples.size
    if sort_idx is None:
        if log_posterior is None:
            raise ValueError("If sorting indicies are not provided, the log posterior must be given.")
        sort_idx = np.argsort(log_posterior)[::-1]        
        
    for i in range(3, len(sort_idx)):
        (l_this, u_this), coverage_1d_this, coverage_nd_this = compute_sample_range_and_coverage(samples, weights, sort_idx[:i])

        if coverage_1d_this > coverage_1d_threshold:
            (l_last, u_last), coverage_1d_last, coverage_nd_last = compute_sample_range_and_coverage(samples, weights, sort_idx[:i-1])

            # Interpolate between this sample and the last at the target coverage
            t = (coverage_1d_this - coverage_1d_threshold)/(coverage_1d_this-coverage_1d_last)
            l = l_this*(1-t) + l_last*t
            u = u_this*(1-t) + u_last*t
            coverage_1d = coverage_1d_this*(1-t) + coverage_1d_last*t
            coverage_nd = coverage_nd_this*(1-t) + coverage_nd_last*t
            
            if verbose: print(f"1D coverage change in last sample: {coverage_1d_last:.2f}->{coverage_1d_this:.2f}")
            return (l,u), coverage_1d, coverage_nd, i
    else:
        warnings.warn(f"Could not match 1D covarage of {coverage_1d_threshold}. Got coverage of {coverage_1d_this:.2f} for CI ({l_this:.3f}, {u_this:.3f}).")
        return (l_this, u_this), coverage_1d_this, coverage_nd_this, i

def plot_CI(plot, chains, params, MAP=None, CI=None):
    """Plot CI bands on a getdist corner plot.
    
    Arguments
    ---------
    plot : getdist.plots.GetDistPlotter
        Plotter to add bands to.
    chains : list
        List of chains (getdist.MCSamples) that are in the plot.
    params : list
        List of parameters.
    MAP : dict
        Dictionary of MAP. Needs to have structure 
        {chain.name_tag : {param_name : MAP}}.
    CI : dict
        Dictionary of CI. Needs to have structure 
        {chain.name_tag : {param_name : (l, u)}}.
    """
    for chain_idx, chain in enumerate(chains[::-1]):
        if MAP is not None:
            markers = {p : None if p not in MAP[chain.name_tag] else MAP[chain.name_tag][p] for p in params}
        if CI is not None:
            HPD_range = {p : None if p not in CI[chain.name_tag] else CI[chain.name_tag][p] for p in params}
        
        for i in range(plot.subplots.shape[0]):
            for j in range(i+1):
                ax = plot.subplots[i,j]
                param_x, param_y = params[j], params[i]
                if MAP is not None and markers[param_x] is not None:
                    plot.add_x_marker(markers[param_x], ax=ax, c=plot.settings.solid_colors[chain_idx], ls="-")
                if MAP is not None and j < i and markers[param_y] is not None:
                    plot.add_y_marker(markers[param_y], ax=ax, c=plot.settings.solid_colors[chain_idx], ls="-")
                
                if CI is not None and HPD_range[param_x] is not None and i == j:
                    ax.axvspan(*HPD_range[param_x], alpha=0.5, facecolor=plot.settings.solid_colors[chain_idx])


parameter_dictionary = {
        "omega_m" :                  {"cosmosis" :    "cosmological_parameters--omega_m",
                                      "montepython" : "Omega_m",
                                      "cosmomc" :     "omegam",
                                      "latex" :       "\\Omega_m"},
        "omega_c h^2" :              {"cosmosis" :    "cosmological_parameters--omch2",
                                      "montepython" : "omega_cdm",
                                      "cosmomc" :     "omegach2",
                                      "latex" :       "h^2\\Omega_c"},
        "omega_b h^2" :              {"cosmosis" :    "cosmological_parameters--ombh2", 
                                      "montepython" : "omega_b",
                                      "cosmomc" :     "omegabh2",
                                      "latex" :       "h^2\\Omega_b"},
        "h" :                        {"cosmosis" :    "cosmological_parameters--h0", 
                                      "montepython" : "h",
                                      "cosmomc" :     "h",
                                      "latex" :       "h"},
        "w" :                        {"cosmosis" :    "cosmological_parameters--w", 
                                      "montepython" : "w",
                                      "cosmomc" :     "w",
                                      "latex" :       "w"},
        "w_background" :             {"cosmosis" :    "cosmological_parameters--w_background", 
                                      "montepython" : "w",
                                      "cosmomc" :     "w_b",
                                      "latex" :       "w_b"},
        "n_s"   :                    {"cosmosis" :    "cosmological_parameters--n_s",
                                      "montepython" : "n_s",
                                      "cosmomc" :     "ns",
                                      "latex" :       "n_s"},
        "A_s"   :                    {"cosmosis" :    "cosmological_parameters--a_s",
                                      "montepython" : "A_s",
                                      "cosmomc" :     "as",
                                      "latex" :       "A_s"},
        "log A_s"   :                {"cosmosis" :    "cosmological_parameters--ln_1e10_a_s",
                                      "montepython" : "ln10^{10}A_s",
                                      "cosmomc" :     "logA",
                                      "latex" :       "{\\rm ln}(10^{10} A_s)"},
        "sigma_8" :                  {"cosmosis" :    "cosmological_parameters--sigma_8",
                                      "montepython" : "sigma8",
                                      "cosmomc" :     "sigma8",
                                      "latex" :       "\\sigma_8"},
        "S_8" :                      {"cosmosis" :    "cosmological_parameters--s_8",
                                      "montepython" : "S8",
                                      "cosmomc" :     "s8",
                                      "latex" :       "S_8"},
        "m_nu" :                     {"cosmosis" :    "cosmological_parameters--mnu",
                                      "cosmomc" :     "mnu",
                                      "latex" :       "M_\\nu"},
        "omega_k" :                  {"cosmosis" :    "cosmological_parameters--omega_k",
                                      "cosmomc" :     "omegak",
                                      "latex" :       "\\Omega_k"},
        "omega_Lambda" :             {"cosmosis" :    "cosmological_parameters--omega_lambda",
                                      "cosmomc" :     "omegal",
                                      "latex" :       "\\Omega_\\Lambda"},
        "omega_nu" :                 {"cosmosis" :    "cosmological_parameters--omega_nu",
                                      "cosmomc" :     "omeganu",
                                      "latex" :       "\\Omega_\\nu"},
        "b_1 sigma_8 S_8 lowz" :     {"cosmosis" :    "cosmological_parameters--bsigma8S8_bin_1",
                                      "cosmomc" :     "b1l_sigma8_s8",
                                      "latex" :       "b_1^{\\rm lowz}\\sigma_8 S_8"},
        "b_1 sigma_8 S_8 highz" :    {"cosmosis" :    "cosmological_parameters--bsigma8S8_bin_2",
                                      "cosmomc" :     "b1h_sigma8_s8",
                                      "latex" :       "b_1^{\\rm highz}\\sigma_8 S_8"},
        "cosmomc_theta" :            {"cosmosis" :    "cosmological_parameters--cosmomc_theta",
                                      "cosmomc" :     "theta",
                                      "latex" :       "100\\theta_{\\rm MC}"},
        "tau" :                      {"cosmosis" :    "cosmological_parameters--tau",
                                      "cosmomc" :     "tau",
                                      "latex" :       "\\tau"},
        "b_1 lowz" :                 {"cosmosis" :    "bias_parameters--b1_bin_1",
                                      "cosmomc" :     "b1l",
                                      "latex" :       "b_1^{\\rm lowz}"},
        "b_2 lowz" :                 {"cosmosis" :    "bias_parameters--b2_bin_1",
                                      "cosmomc" :     "b2l",
                                      "latex" :       "b_2^{\\rm lowz}"},
        "gamma_3 lowz" :             {"cosmosis" :    "bias_parameters--gamma3_bin_1",
                                      "cosmomc" :     "gamma3l",
                                      "latex" :       "\\gamma_3^{\\rm lowz}"},
        "a_vir lowz" :               {"cosmosis" :    "bias_parameters--a_vir_bin_1",
                                      "cosmomc" :     "a_virl",
                                      "latex" :       "a_{\\rm vir}^{\\rm lowz}"},
        "b_1 highz" :                {"cosmosis" :    "bias_parameters--b1_bin_2",
                                      "cosmomc" :     "b1h",
                                      "latex" :       "b_1^{\\rm highz}"},
        "b_2 highz" :                {"cosmosis" :    "bias_parameters--b2_bin_2",
                                      "cosmomc" :     "b2h",
                                      "latex" :       "b_2^{\\rm highz}"},
        "gamma_3 highz" :            {"cosmosis" :    "bias_parameters--gamma3_bin_2",
                                      "cosmomc" :     "gamma3h",
                                      "latex" :       "\\gamma_3^{\\rm highz}"},
        "a_vir highz" :              {"cosmosis" :    "bias_parameters--a_vir_bin_2",
                                      "cosmomc" :     "a_virh",
                                      "latex" :       "a_{\\rm vir}^{\\rm highz}"},
        "A_IA" :                     {"cosmosis" :    "intrinsic_alignment_parameters--a",
                                      "montepython" : "A_IA",
                                      "cosmomc" :     "a_ia",
                                      "latex" :       "A_{\\rm IA}"},
        "A_baryon" :                 {"cosmosis" :    "halo_model_parameters--a",
                                      "montepython" : "c_min",
                                      "cosmomc" :     "a_baryon",
                                      "latex" :       "A_{\\rm baryon}"},
        "eta_baryon" :               {"cosmosis" :    "halo_model_parameters--eta0",
                                      "montepython" : "eta_0",
                                      "cosmomc" :     "eta_baryon",
                                      "latex" :       "\\eta_{\\rm baryon}"},
        "delta_c" :                  {"cosmosis" :    "shear_c_bias--delta_c",
                                      "montepython" : "dc",
                                      "cosmomc" :     "delta_c",
                                      "latex" :       "\\delta c"},
        "A_c" :                      {"cosmosis" :    "shear_c_bias--a_c",
                                      "montepython" : "Ac",
                                      "cosmomc" :     "a_c",
                                      "latex" :       "A_{\\rm c}"},
        "delta z_1" :                {"cosmosis" :    "nofz_shifts--bias_1",
                                      "montepython" : "D_z1",
                                      "cosmomc" :     "delta_z1",
                                      "latex" :       "\\delta z_1"},
        "delta z_2" :                {"cosmosis" :    "nofz_shifts--bias_2",
                                      "montepython" : "D_z2",
                                      "cosmomc" :     "delta_z2",
                                      "latex" :       "\\delta z_2"},
        "delta z_3" :                {"cosmosis" :    "nofz_shifts--bias_3",
                                      "montepython" : "D_z3",
                                      "cosmomc" :     "delta_z3",
                                      "latex" :       "\\delta z_3"},
        "delta z_4" :                {"cosmosis" :    "nofz_shifts--bias_4",
                                      "montepython" : "D_z4",
                                      "cosmomc" :     "delta_z4",
                                      "latex" :       "\\delta z_4"},
        "delta z_5" :                {"cosmosis" :    "nofz_shifts--bias_5",
                                      "montepython" : "D_z5",
                                      "cosmomc" :     "delta_z5",
                                      "latex" :       "\\delta z_5"},
        "A_Planck" :                 {"cosmosis" :    "planck--a_planck",
                                      "cosmomc" :     "calPlanck",
                                      "latex" :       "A_{\\rm Planck}"},
        "F_AP lowz" :                {"cosmosis" :    "lss_parameters--f_ap_bin_1",
                                      "cosmomc" :     "FAP1",
                                      "latex" :       "F_{\\rm AP lowz}"},
        "rs_DV lowz" :               {"cosmosis" :    "lss_parameters--rs_dv_bin_1",
                                      "cosmomc" :     "rsDv1",
                                      "latex" :       "r_{\\rm d}/D_V {\\rm lowz}"},
        "fsigma_8 lowz" :            {"cosmosis" :    "lss_parameters--fsigma_8_bin_1",
                                      "cosmomc" :     "fsigma8z1",
                                      "latex" :       "f\\sigma_8 {\\rm highz}"},
        "F_AP highz" :               {"cosmosis" :    "lss_parameters--f_ap_bin_2",
                                      "cosmomc" :     "FAP2",
                                      "latex" :       "F_{\\rm AP highz}"},
        "rs_DV highz" :              {"cosmosis" :    "lss_parameters--rs_dv_bin_2",
                                      "cosmomc" :     "rsDv2",
                                      "latex" :       "r_{\\rm d}/D_V {\\rm lowz}"},
        "fsigma_8 highz" :           {"cosmosis" :    "lss_parameters--fsigma_8_bin_2",
                                      "cosmomc" :     "fsigma8z2",
                                      "latex" :       "f\\sigma_8 {\\rm highz}"},
        "lnlike" :                   {"cosmosis" :    "like",
                                      "cosmomc" :     "loglike",
                                      "latex" :       "\\log \\mathcal{L}"},
        "lnpost" :                   {"cosmosis" :    "post",
                                      "cosmomc" :     "logpost",
                                      "latex" :       "\\log \\mathcal{P}"},
        "lnprior" :                  {"cosmosis" :    "prior",
                                      "cosmomc" :     "logprior",
                                      "latex" :       "\\log \\pi"},
        "weight" :                   {"cosmosis" :    "weight",
                                      "cosmomc" :     "weight",
                                      "latex" :       "w"},
           
                              }

def parameter_mapping(chain_file,
                      parameters,
                      statistics=None, chain_format="cosmosis", 
                      parameter_name_map=None):
    chain_format = chain_format.lower()
    parameter_names = parameter_dictionary
    
    with open(chain_file, "r") as f:
        params = f.readline()[1:]
    
    if chain_format == "cosmosis":
        params = [p.strip().lower() for p in params.split("\t")]
    elif chain_format == "montepython":
        params = [p.strip() for p in params.split(",")]
    else:
        raise ValueError(f"Chain format {chain_format} not supported.")

    param_idx = []
    output_parameter_names = {mapping : [] for mapping in parameter_name_map}
    for p in parameters:
        if p in parameter_names and parameter_names[p][chain_format] in params:
            param_idx.append(params.index(parameter_names[p][chain_format]))
            if parameter_name_map is not None:
                for mapping in parameter_name_map:
                    output_parameter_names[mapping].append(parameter_names[p][mapping])
        else:
            raise ValueError(f"Parameter {p} not in chain.")
            
    stats_idx = []
    statistics = statistics or []
    for s in statistics:
        s = s.lower()
        if s in params:
            stats_idx.append(params.index(s))
        else:
            raise ValueError(f"Statistic {s} not in chain.")
    
    if parameter_name_map is not None:
        return (param_idx, stats_idx) + tuple(output_parameter_names.values())
    else:
        return param_idx, stats_idx

def load_nested_sampling_file(filename):
    info = {}
    with open(filename, "r") as f:
        for s in f.readlines()[-3:]:
            m = re.match("^#([a-z_]+)=([0-9\.\-]+)", s)
            if m is None:
                raise ValueError(f"Can't read property from line {s}")
            key, value = m.groups()
            info[key] = value

    n_sample = int(info["nsample"])
    chain = np.loadtxt(filename)
    print(f"Using {n_sample} samples out of {chain.shape[0]} in the chain.")
    return chain[-n_sample:], float(info["log_z"]), float(info["log_z_error"])


class MissingParameterError(KeyError):
    pass

def get_sampled_params_and_ranges(values_file, only_varied=True):
    config = configparser.ConfigParser()
    with open(values_file, "r") as f:
        config.read_file(f)

    params = []
    for sec in config.sections():
        for key in config[sec].keys():
            val = config[sec][key]
            val = [float(s) for s in config[sec][key].split()]
            if len(val) == 3:
                params.append((sec, key, [val[0], val[2]]))
            
    return params

def cosmosis_to_cosmomc_param_names(param_names):
    cosmomc_params = []
    for param in param_names:
        for m in parameter_dictionary.values():
            if m["cosmosis"] == "--".join(param).lower():
                cosmomc_params.append(m["cosmomc"])
                break
        else:
            cosmomc_params.append("--".join(param).lower())
        
    return cosmomc_params

def load_chain(chain_file, parameters=None, run_name=None, 
               chain_format="cosmosis", parameter_map="cosmomc", strict_mapping=False, 
               values=None, burn_in=0.3):
    with open(chain_file, "r") as f:
        params = f.readline()[1:]
    
    if chain_format == "cosmosis":
        chain_params = [p.strip().lower() for p in params.split("\t")]
    elif chain_format == "montepython":
        chain_params = [p.strip() for p in params.split(",")]
    else:
        raise ValueError(f"Chain format {chain_format} not supported.")
        
    parameter_names = []
    parameter_names_latex = []
    for p in chain_params:
        for mapping in parameter_dictionary.values():
            if mapping[chain_format] == p:
                parameter_names.append(mapping[parameter_map])
                parameter_names_latex.append(mapping["latex"])
                break
        else:
            if strict_mapping:
                raise MissingParameterError(f"Parameter {p} in chain does not have mapping to {parameter_map} format.")
            else:
                warnings.warn(f"Parameter {p} in chain does not have mapping to {parameter_map} format.")
                parameter_names.append(p)
                parameter_names_latex.append(p)

    column_idx = list(range(len(parameter_names)))
    
    stat_column_idx = {}
    
    # remove statistics columns
    for s in ["weight", "lnlike"]:
        s_map = parameter_dictionary[s][parameter_map]
        if s_map in parameter_names:
            stat_column_idx[s] = column_idx[parameter_names.index(s_map)]
            
    # remove statistic columns from parameter list
    for s, i in stat_column_idx.items():
        del parameter_names[i]
        del parameter_names_latex[i]
        del column_idx[i]
    
    try:
        chain, log_Z, log_Z_err = load_nested_sampling_file(chain_file)
        nested_sampling = True
    except:
        chain = np.loadtxt(chain_file)
        n_sample = chain.shape[0]
        if isinstance(burn_in, float) and burn_in > 0.0 and burn_in < 1.0:
            chain = chain[int(n_sample*burn_in):]
        elif isinstance(burn_in, int):
            chain = chain[burn_in:]
        else:
            raise ValueError(f"Invalid burn_in value: {burn_in}")
        print(f"Using {chain.shape[0]} samples out of {n_sample} in the chain.")
        nested_sampling = False
        
    if values:
        if chain_format != "cosmosis" or parameter_map != "cosmomc":
            raise ValueError("Loading value files is only supported for cosmosis chains and cosmomc parameters.")
        sampled_params_ranges = get_sampled_params_and_ranges(values)
        cosmomc_names = cosmosis_to_cosmomc_param_names([(s,k) for s,k,_ in sampled_params_ranges])
        ranges = {cosmomc_names[i] : sampled_params_ranges[i][2] for i in range(len(cosmomc_names)) if cosmomc_names[i] in parameter_names}
    else:
        ranges = {}
    
    run_name = run_name or os.path.split(chain_file)[1]
    samples = getdist.MCSamples(name_tag=run_name,
                                samples=chain[:,column_idx],
                                weights=chain[:,stat_column_idx["weight"]] if "weight" in stat_column_idx else None,
                                loglikes=chain[:,stat_column_idx["lnlike"]] if "lnlike" in stat_column_idx else None,
                                names=parameter_names,
                                labels=parameter_names_latex,
                                sampler="nested" if nested_sampling else None,
                                ranges=ranges)
    
    if nested_sampling:
        samples.log_Z = log_Z
        samples.log_Z_err = log_Z_err
    else:
        samples.chain_offsets = [0, chain.shape[0]]
    
    return samples