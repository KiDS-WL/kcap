import warnings

import numpy as np

import scipy.optimize
import scipy.stats


def find_CI(method, samples, weights=None, coverage=0.683,
            logpost=None, logpost_sort_idx=None,
            return_point_estimate=False, return_coverage=False, 
            return_extras=False, options=None):
    """Compute credible intervals and point estimates from samples.

    Arguments
    ---------
    method : str
        Method to compute CI. Options are "PJ-HPD", "tail CI", "std", and "HPD".
        PJ-HPD:  Compute the CI from the joint posterior HPD region such that
                 the projected range of the HPDR has coverage ``coverage``.
                 See Joachimi et al. 2020.
                 The point estimate is the joint posterior MAP.
        tail CI: This is the usual quantile CI. I.e., for CI (l,u) and
                 coverage c, P(x<l) = (1-c)/2 and P(x>u) = 1-(1-c)/2.
                 The point estimate is the median.
        std:     Compute the CI as (mean - n_sigma*std, mean + n_sigma*std).
                 ``n_sigma`` is the number of standard devations that cover
                 ``coverage`` in a normal distribution.
                 The point estimate is the mean.
        HPD:     Compute the HPDI of the samples.
                 The point estimate is the MAP.
    samples : array
        Samples to use.
    weights : array, optional
        Sample weights.
    coverage : float, optional
        Target coverage. This gets converted into sigmas. Default: 0.683.
    logpost : array, optional
        Array of the log posterior values of the samples. Required for method ``PJ-HPD``.
    logpost_sort_idx : array, optional
        Array of indices that sort the samples in descending posterior value.
        If method is ``PJ-HPD`` and it is not provided, this will be computed
        internally from logpost.
    return_point_estimate : bool, optional
        Whether to return the point_estimate.
    return_coverage : bool, optional
        Whether to return the actual coverage of the CI.
    options : dict, optional
        Additional options passed to the CI methods.


    Returns
    -------
    (l, u) : tuple
        Credible interval of the samples.
    p : float
        Point estimate. Only returned if return_point_estimate is true.
    coverage : float
        The achieved coverage of the returned CI.
    """
    options = options or {}

    extras = None
    if method.lower() == "pj-hpd" or method.lower() == "projected joint hpd":
        if logpost is None and logpost_sort_idx is None:
            raise ValueError("For method PJ-HPD, either logpost or logpost_sort_idx need to be specified.")
            
        CI, MAP, alpha, n_sample = find_projected_joint_HPDI(samples, weights, 
                                                             coverage_1d_threshold=coverage,
                                                             sort_idx=logpost_sort_idx, log_posterior=logpost, 
                                                             return_map=True, return_coverage_1d=True, return_n_sample=True,
                                                             **options)
        point_estimate = MAP
        extras = n_sample
        
    elif method.lower() == "hpd" or method.lower() == "m-hpd":
        CI, marg_MAP, alpha, no_constraints = find_marginal_HPDI(samples, weights, coverage=coverage, 
                                                                 return_map=True, return_coverage=True,
                                                                 check_prior_edges=True,
                                                                 **options)
        point_estimate = marg_MAP
        extras = no_constraints
        
    elif method.lower() == "tail ci" or method.lower() == "quantile ci":
        CI, marg_median, alpha = find_quantile_CI(samples, weights, coverage=coverage, 
                                                  return_median=True, return_coverage=True)
        point_estimate = marg_median

    elif method.lower() == "std":
        CI, marg_mean, alpha = find_std_CI(samples, weights, coverage=coverage, 
                                           return_mean=True, return_coverage=True)
        point_estimate = marg_mean
        
    else:
        raise NotImplementedError(f"Method {method} not supported.")
        
    result = [CI]
    if return_point_estimate:
        result += [point_estimate]
    if return_coverage:
        result += [alpha]
    if return_extras:
        result += [extras]

    if len(result) == 1: 
        # Only CI
        return result[0]
    else:
        return tuple(result)


def find_std_CI(samples, weights=None, coverage=0.683, 
                return_mean=False, return_coverage=False):
    """Find mean ± std of samples and compute the coverage of that interval.

    Arguments
    ---------
    samples : array
        Samples to use.
    weights : array, optional
        Sample weights.
    coverage : float, optional
        Target coverage. This gets converted into sigmas. Default: 0.683.
    return_mean : bool, optional
        Whether to return the mean.
    return_coverage : bool, optional
        Whether to return the actual coverage of the CI.

    Returns
    -------
    (l, u) : tuple
        The CI given by (mean - n_sigma*std, mean + n_sigma*std).
    mean : float
        Mean of the samples. Only returned if return_mean is true.
    coverage : float
        Coverage of the CI.
    """
    if samples.ndim != 1:
        raise ValueError("Sample array not 1d.")

    std = np.sqrt(np.cov(samples, aweights=weights, ddof=1))
    mean = np.average(samples, weights=weights)

    n_sigma = -scipy.stats.norm.ppf((1-coverage)/2)
    l = mean - n_sigma*std
    u = mean + n_sigma*std

    result = [(l, u)]
    if return_mean:
        result += [mean]
    if return_coverage:
        cdf = create_interpolated_cdf(samples, weights, threshold=samples.min())
        alpha = cdf(u) - cdf(l)
        result += [alpha]

    if len(result) == 1:
        # Only CI
        return result[0]
    else:
        return tuple(result)


def find_quantile_CI(samples, weights=None, coverage=0.683, 
                     return_median=False, return_coverage=False):
    """Find quantile/tail credible interval.

    Arguments
    ---------
    samples : array
        Samples to use.
    weights : array, optional
        Sample weights. If not specifed is set to 1/len(samples).
    coverage : float
        Target coverage of the returned interval.
    return_median : bool, optional
        Whether to return the median.
    return_coverage : bool, optional
        Whether to return the actual coverage of the CI.

    Returns
    -------
    (l, u) : tuple
        The quantile/tail credible interval.
    median : float
        Median of the samples. Only returned if return_median is true.
    coverage : float
        The achieved coverage of the returned CI.
    """
    if coverage <= 0.0 or coverage >= 1.0:
        raise ValueError("Requested coverage outside of range (0,1)")

    if weights is None:
        weights = np.ones_like(samples)/samples.size

    samples_min = samples.min().squeeze()
    samples_max = samples.max().squeeze()

    cdf = create_interpolated_cdf(samples, weights, threshold=samples.min())

    # Get lower interval bound
    c = (1-coverage)/2
    res = scipy.optimize.root_scalar(lambda x: cdf(x)-c, bracket=(samples_min, samples_max))
    if not res.converged:
        raise RuntimeError(f"Failed to find lower bound for tail coverage {c}.")
    l = res.root

    # Get upper interval bound
    c = 1-(1-coverage)/2
    res = scipy.optimize.root_scalar(lambda x: cdf(x)-c, bracket=(samples_min, samples_max))
    if not res.converged:
        raise RuntimeError(f"Failed to find upper bound for tail coverage {c}.")
    u = res.root

    alpha = cdf(u) - cdf(l)

    result = [(l, u)]
    if return_median:
        res = scipy.optimize.root_scalar(lambda x: cdf(x)-0.5, bracket=(samples_min, samples_max))
        if not res.converged:
            raise RuntimeError(f"Failed to find median.")
        median = res.root
        result += [median] 
    if return_coverage:
        result += [alpha]

    if len(result) == 1:
        # Only CI
        return result[0]
    else:
        return tuple(result)



def find_marginal_HPDI(samples, weights=None, coverage=0.683,
                       kde_bandwidth="silverman",
                       return_map=False, return_coverage=False,
                       check_prior_edges=False, prior_edge_threshold=0.1):
    """Find the highest posterior density interval of a marginal posterior.

    Arguments
    ---------
    samples : array
        Array of the samples.
    weights : array, optional
        Sample weights.
    coverage : float
        Coverage of the returned CI. Default 0.683.
    kde_bandwidth : str or float
        Bandwidth used for KDE. Same as scipy.stats.gaussian_kde.
    return_map : bool, optional
        Whether to return the marginal MAP.
    return_coverage : bool, optional
        Whether to return the actual coverage of the CI.

    Returns
    -------
    (l, u) : tuple
        HPD credible interval.
    map : float
        MAP of the samples. Only returned if return_map is true.
    alpha : float
        Achieved coverage of the CI. Only returned if return_coverage is true.

    Based on Chieh-An Lin's code.
    """
    if coverage <= 0.0 or coverage >= 1.0:
        raise ValueError("Requested coverage outside of range (0,1)")

    if samples.ndim != 1:
        raise ValueError("Samples need to be 1D marginals.")
    if weights is not None and weights.shape != samples.shape:
        raise ValueError(f"Samples and weights have mismatching shapes: {samples.shape} vs {weights.shape}.")

    samples_min = samples.min().squeeze()
    samples_max = samples.max().squeeze()

    kde = scipy.stats.gaussian_kde(samples, bw_method=kde_bandwidth, weights=weights)

    # Find maximum of marginal posterior
    res = scipy.optimize.minimize_scalar(lambda x: -kde(x), bracket=(samples_min, samples_max))
    if not res.success:
        raise RuntimeError(f"Failed to find maximum marginal posterior: {res}")

    marg_MAP = res.x.squeeze()
    marg_MAP_value = -res.fun.squeeze()

    marg_samples_min = kde(samples_min)
    marg_samples_max = kde(samples_max)

    # Find CI
    def get_coverage(level, return_interval=False):        
        # Find locations of marginal posterior where P(x) = level
        if level < marg_samples_min:
            l = samples_min
        else:
            res = scipy.optimize.root_scalar(lambda x: kde(x)-level, 
                                             bracket=(samples_min, marg_MAP))
            if not res.converged:
                raise RuntimeError(f"Failed to find lower bound for level {level}.")
            l = res.root

        if level < marg_samples_max:
            u = samples_max
        else:
            res = scipy.optimize.root_scalar(lambda x: kde(x)-level, 
                                             bracket=(marg_MAP, samples_max))
            if not res.converged:
                raise RuntimeError(f"Failed to find upper bound for level {level}.")
            u = res.root

        # Get coverage between l and u
        alpha = kde.integrate_box_1d(l, u)

        if return_interval:
            return alpha, (l, u)
        else:
            return alpha

    res = scipy.optimize.root_scalar(lambda x: get_coverage(x)-coverage, 
                                     bracket=(0, marg_MAP_value))
    if not res.converged:
        raise RuntimeError(f"Failed to find coverage.")

    level = res.root
    alpha, CI = get_coverage(level, return_interval=True)

    if check_prior_edges:
        if marg_samples_min/marg_MAP_value > prior_edge_threshold or marg_samples_max/marg_MAP_value > prior_edge_threshold:
            prior_edge_hit = True
        else:
            prior_edge_hit = False

    result = [CI]
    if return_map:
        result += [marg_MAP]
    if return_coverage:
        result += [alpha]
    if check_prior_edges:
        result += [prior_edge_hit]

    if len(result) == 1:
        # Only CI
        return result[0]
    else:
        return tuple(result)


def compute_sample_range_and_coverage(samples, weights, idx):
    """Compute range of values of samples[idx] and their coverage."""
    l, u = samples[idx].min(), samples[idx].max()
    coverage_1d = np.sum(weights[(l <= samples) & (samples <= u)])
    coverage_nd = np.sum(weights[idx])
    return (l,u), coverage_1d, coverage_nd

def create_interpolated_cdf(samples, weights, 
                            threshold=None, reverse=False,
                            start_at_zero=True):
    unique_samples, unique_inverse_idx = np.unique(samples, return_inverse=True)
    if len(unique_samples) != len(samples):
        unique_weights = np.zeros_like(unique_samples)
        np.add.at(unique_weights, unique_inverse_idx, weights)
        samples = unique_samples
        weights = unique_weights

    if threshold is None:
        selection = np.ones_like(samples, dtype=bool)
    else:
        selection = samples <= threshold if reverse else samples >= threshold
    if reverse:
        sample_sort_idx = np.argsort(samples[selection])
        cdf = np.cumsum(weights[selection][sample_sort_idx])
        cdf = cdf[-1] - cdf
        if not start_at_zero:
            cdf += weights[selection][sample_sort_idx[-1]]
    else:
        sample_sort_idx = np.argsort(samples[selection])
        cdf = np.cumsum(weights[selection][sample_sort_idx])
        if start_at_zero:
            cdf -= cdf[0]
    
    if len(cdf) < 2:
        # Only one sample
        return lambda x: np.zeros_like(x)

    cdf_func = scipy.interpolate.InterpolatedUnivariateSpline(x=samples[selection][sample_sort_idx],
                                                              y=cdf, k=1, ext=3)
    return cdf_func

def weighted_median(a, weights):
    if len(a) == 1:
        return a
    if np.all(np.isclose(a, a[0])):
        return a[0]

    if weights.sum()/weights.max()-1 < 1e-2:
        return a[np.argmax(weights)]

    sort_idx = np.argsort(a)
    
    l_cum = create_interpolated_cdf(a, weights, start_at_zero=False)
    r_cum = create_interpolated_cdf(a, weights, reverse=True, start_at_zero=False)
    
    return scipy.optimize.root_scalar(lambda x: l_cum(x) - r_cum(x), bracket=(a.min(),a.max())).root


def find_projected_joint_HPDI(samples, weights=None, coverage_1d_threshold=0.683, 
                              MAP=None,
                              sort_idx=None, log_posterior=None, 
                              method="interpolate", twosided=True, verbose=False, strict=False,
                              return_map=False, return_coverage_1d=False,
                              return_coverage_nd=False, return_n_sample=False):
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
        Target 1D coverage of the CI. (Default 0.683)
    MAP : float, optional
        If the exact MAP is known, it can be specified here. If not provided, 
        the sample with the highest posterior value is used. 
    sort_idx : numpy.array
        Array of indices that sort the samples in descending posterior value.
        If not provided, this will be computed internally from log_posterior.
    log_posterior : numpy.array
        Array of the log posterior values of the samples.
    method : str
        How to interplate the MHPD CI for finite number of samples. Options are
        'interpolate', 'expand', 'expand symmetric', 'expand minimal'.
        Default 'interpolate'
    twosided: bool
        Whether to keep search samples until both the lower and upper bounds
        extend beyond the threshold coverage. Default: True.
    verbose : bool
        Print extra outputs.
    strict : bool
        Whether to raise an exception of no CI can be found. Default: False.
    return_map : bool, optional
        Whether to return the MAP of the samples.
    return_coverage_1d : bool, optional
        Whether to return the actual, 1D coverage of the CI.
    return_coverage_nd : bool, optional
        Whether to return the nd coverage of the HPD used to compute the CI.
    return_n_sample : bool, optional
        Whether to return the number of (joint posterior) samples used to compute the CI.

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
    if coverage_1d_threshold <= 0.0 or coverage_1d_threshold >= 1.0:
        raise ValueError("Requested coverage outside of range (0,1)")

    if weights is None:
        weights = np.ones_like(samples)/samples.size
    if sort_idx is None:
        if log_posterior is None:
            raise ValueError("If sorting indicies are not provided, the log posterior must be given.")
        sort_idx = np.argsort(log_posterior)[::-1]

    if MAP is not None:
        samples = np.insert(samples, 0, MAP)
        weights = np.insert(weights, 0, 0.0)
        sort_idx = np.insert(sort_idx + 1, 0, 0)

    if samples[sort_idx[0]] <= samples.min():
        samples[sort_idx[0]] = samples.min()
        twosided = False
        warnings.warn("The starting point is at or lower than the sample range. Only using one sided interpolation.")
    elif samples[sort_idx[0]] >= samples.max():
        samples[sort_idx[0]] = samples.max()
        twosided = False
        warnings.warn("The starting point is at or higher than the sample range. Only using one sided interpolation.")


    l_outer = l_inner = u_inner = u_outer = samples[sort_idx[0]]
    coverage_1d_inner = coverage_nd_inner = 0
    n_inner = 0
    for i in range(2, len(sort_idx)):
        (l_this, u_this), coverage_1d_this, coverage_nd_this = compute_sample_range_and_coverage(samples, weights, sort_idx[:i])

        if coverage_1d_this < coverage_1d_threshold:
            # Still under the target coverage
            l_inner = l_outer = l_this
            u_inner = u_outer = u_this
            coverage_1d_inner = coverage_1d_this
            coverage_nd_inner = coverage_nd_this
            n_inner = i
        else:
            # Over the target coverage. Only update outer limits
            if l_this < l_inner:
                l_outer = l_this
            if u_this > u_inner:
                u_outer = u_this

        if method == "interpolate":
            if (twosided and (l_inner != l_outer and u_inner != u_outer)) \
            or (not twosided and (l_inner != l_outer or u_inner != u_outer)):
                # Got two estimates for both sides

                l_cdf_func = create_interpolated_cdf(samples, weights, l_inner, reverse=True)
                u_cdf_func = create_interpolated_cdf(samples, weights, u_inner, reverse=False)

                def coverage(t, l_cdf, u_cdf):
                    l = l_outer*(1-t) + l_inner*t
                    u = u_outer*(1-t) + u_inner*t
                    return coverage_1d_inner + l_cdf(l) + u_cdf(u)

                res = scipy.optimize.root_scalar(f=lambda t, l_cdf=l_cdf_func, u_cdf=u_cdf_func: coverage(t, l_cdf, u_cdf)-coverage_1d_threshold,
                                                 bracket=(0,1))
                t = res.root
                l = l_outer*(1-t) + l_inner*t
                u = u_outer*(1-t) + u_inner*t
                coverage_1d = np.sum(weights[(l <= samples) & (samples <= u)])
                break

        elif method == "expand symmetric":
            if(l_inner != l_outer or u_inner != u_outer):
                delta_init = max(l_inner-l_outer, u_outer-u_inner)

                l_cdf_func = create_interpolated_cdf(samples, weights, l_inner, reverse=True)
                u_cdf_func = create_interpolated_cdf(samples, weights, u_inner, reverse=False)

                def coverage(d, l_cdf, u_cdf):
                    l = l_inner - d
                    u = u_inner + d
                    return coverage_1d_inner + l_cdf(l) + u_cdf(u)

                res = scipy.optimize.root_scalar(f=lambda t, l_cdf=l_cdf_func, u_cdf=u_cdf_func: coverage(t, l_cdf, u_cdf)-coverage_1d_threshold,
                                                 bracket=(0,delta_init))
                d = res.root
                l = l_inner - d
                u = u_inner + d
                coverage_1d = np.sum(weights[(l <= samples) & (samples <= u)])
                break

        elif method == "expand":
            if(l_inner != l_outer or u_inner != u_outer):
                # Only need two estimates for one side.

                delta_coverage = (coverage_1d_threshold-coverage_1d_inner)/2
                sample_sort_idx = np.argsort(samples[samples < l_inner])[::-1]
                l_idx = np.searchsorted(np.cumsum(weights[samples < l_inner][sample_sort_idx]), delta_coverage)

                if l_idx == samples[samples < l_inner].size:
                    l_idx -= 1
                l = samples[samples < l_inner][sample_sort_idx][l_idx]

                sample_sort_idx = np.argsort(samples[samples > u_inner])
                u_idx = np.searchsorted(np.cumsum(weights[samples > u_inner][sample_sort_idx]), delta_coverage)
                if u_idx == samples[samples > u_inner].size:
                    u_idx -= 1
                u = samples[samples > u_inner][sample_sort_idx][u_idx]

                coverage_1d = np.sum(weights[(l <= samples) & (samples <= u)])
                break

        elif method == "expand minimal":
            if(l_inner != l_outer and u_inner != u_outer):
                # Initial guess
                l_delta = l_inner-l_outer
                u_delta = u_outer-u_inner

                l_cdf_func = create_interpolated_cdf(samples, weights, l_inner, reverse=True)
                u_cdf_func = create_interpolated_cdf(samples, weights, u_inner, reverse=False)

                def constraint(x, l_cdf, u_cdf):
                    l, u = x
                    coverage = coverage_1d_inner + l_cdf(l) + u_cdf(u)
                    return coverage - coverage_1d_threshold

                def constraint_jac(x, l_cdf, u_cdf):
                    l, u = x
                    return l_cdf(l, nu=1), u_cdf(u, nu=1)

                def fun(x):
                    l, u = x
                    return u-l

                def fun_jac(x):
                    l, u = x
                    return -1.0, 1.0

                res = scipy.optimize.minimize(fun=fun, x0=[l_inner-l_delta/2, u_inner+u_delta/2],
                                              jac=fun_jac,
                                              bounds=[(l_outer, l_inner), (u_inner, u_outer)],
                                              constraints=[{"type" : "eq", "fun" : constraint, "jac" : constraint_jac, "args" : (l_cdf_func, u_cdf_func)},
                                                          ],
                                              options={"ftol" : 0.05},
                                              method="SLSQP")
                if not res.success:
                    if res.message != "Positive directional derivative for linesearch":
                        print(res)

                l, u = res.x
                coverage_1d = np.sum(weights[(l <= samples) & (samples <= u)])
                break
        else:
            raise ValueError(f"Method {method} not supported.")      
    else:
        if strict:
            raise RuntimeError(f"Could not match 1D coverage of {coverage_1d_threshold}. Got coverage of {coverage_1d_this:.2f} for CI ({l_this:.3f}, {u_this:.3f}).")
        else:
            warnings.warn(f"Could not match 1D coverage of {coverage_1d_threshold}. Got coverage of {coverage_1d_this:.2f} for CI ({l_this:.3f}, {u_this:.3f}).")
            l, u = l_this, u_this 
            coverage_1d = coverage_1d_this

    result = [(l, u)]
    if return_map:
        if MAP is not None:
            result += [samples[sort_idx[1]]]
        else:
            result += [samples[sort_idx[0]]]
    if return_coverage_1d:
        result += [coverage_1d]
    if return_coverage_nd:
        result += [coverage_nd_inner]
    if return_n_sample:
        result += [n_inner]

    if len(result) == 1:
        # Only CI
        return result[0]
    else:
        return tuple(result)
