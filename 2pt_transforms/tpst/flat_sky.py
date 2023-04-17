import numpy as np

import scipy.integrate
import scipy.special

from . import cosebis

pi = np.pi

def integrate_bessel_functions(func, a, b, bessel_frequencies=[], bessel_orders=[], func_args=(), mode="split"):
    # Split integral at roots of fastest varying Bessel function
    f = max(bessel_frequencies)
    nu = bessel_orders[np.argmax(bessel_frequencies)]
    n_root = int(b/(pi)*f)+1
    
    if mode.lower() == "split":
        roots = scipy.special.jn_zeros(nu, n_root)/f
        roots = roots[:np.searchsorted(roots,b)]
        # Remove every 2nd root to integrate over a full period.
        roots = roots[::2]

        if len(roots) > 0:
            intervals = [(roots[i], roots[i+1]) for i in range(len(roots)-1)]
            intervals.insert(0, (a, roots[0]))
            intervals.append((roots[-1], b))
        else:
            intervals = [(a,b)]
    
        integral = 0
        for interval in intervals:
            integral += scipy.integrate.quad(func, *interval, args=func_args)[0]
    else:
        integral = scipy.integrate.quad(func, a, b, args=func_args, limit=max(50, 2*n_root))[0]

    return integral


def compute_weights(ell, ell_bin, theta_range, apodise_theta=None, apodise_ell=False, estimator="gg", 
                    log_cosebis_roots=None, log_cosebis_normalisations=None,
                    integrator="quad"):
    theta_l, theta_u = theta_range
    if estimator.lower() in ("cosebis", "log-cosebis"):
        n = int(ell_bin)
    else:
        ell_l, ell_u = ell_bin
    
    J_0 = scipy.special.j0
    J_1 = scipy.special.j1
    J_2 = lambda x: scipy.special.jv(2, x)
    J_4 = lambda x: scipy.special.jv(4, x)
    
    if apodise_theta is not None:
        Delta_x, theta_l_apo, theta_u_apo = apodise_theta
        # x_l and x_u should be set independently from theta_l and theta_u
        x_l, x_u = np.log(theta_l), np.log(theta_u)
        theta_l = theta_l_apo
        theta_u = theta_u_apo
        def T(theta):
            x = np.log(theta)
            return np.piecewise(x, condlist=(x < x_l-Delta_x/2, 
                                             (x_l-Delta_x/2 <= x) & (x < x_l+Delta_x/2),
                                             (x_l+Delta_x/2 <= x) & (x < x_u-Delta_x/2),
                                             (x_u-Delta_x/2 <= x) & (x < x_u+Delta_x/2),
                                             (x_u+Delta_x/2 <= x)),
                                   funclist=(0, 
                                             lambda x: np.cos(pi/2*(x-(x_l+Delta_x/2))/Delta_x)**2, 
                                             1, 
                                             lambda x: np.cos(pi/2*(x-(x_u-Delta_x/2))/Delta_x)**2,
                                             0))
    else:
        T = lambda theta: 1
        
    if apodise_ell:
        raise NotImplementedError("Apodising in ell not implemented yet.")
    else:
        g_plus = lambda theta: 1/theta*(ell_u*J_1(ell_u*theta) - ell_l*J_1(ell_l*theta))
        
        G = lambda x: (x-8/x)*J_1(x) - 8*J_2(x)
        g_minus = lambda theta: 1/theta**2*(G(ell_u*theta) - G(ell_l*theta))
        
        h = lambda theta: -1/theta*(ell_u*J_1(ell_u*theta) - ell_l*J_1(ell_l*theta) \
                                    + 2*J_0(ell_u*theta)/theta - 2*J_0(ell_l*theta)/theta)
    
    if estimator.lower() == "gg":
        integral_info = [{"func" : lambda theta: theta*T(theta)*g_plus(theta)*J_0(ell*theta),
                          "bessel_frequencies" : [ell_u, ell_l, ell],
                          "bessel_orders" : [1, 1, 0]}]
    elif estimator.lower() in ("ee", "bb", "eb", "be"):
        integral_info = [{"func" : lambda theta: theta*T(theta)*g_plus(theta)*J_0(ell*theta),
                          "bessel_frequencies" : [ell_u, ell_l, ell_u, ell_l, ell],
                          "bessel_orders" : [1, 1, 2, 2, 0]},
                         {"func" : lambda theta: theta*T(theta)*g_minus(theta)*J_4(ell*theta),
                          "bessel_frequencies" : [ell_u, ell_l, ell_u, ell_l, ell],
                          "bessel_orders" : [1, 1, 2, 2, 4]},]
        if estimator.lower() in ("eb", "be"):
            integral_info[1]["sign"] = -1
    elif estimator.lower() == "ge":
        integral_info = [{"func" : lambda theta: theta*T(theta)*h(theta)*J_2(ell*theta),
                          "bessel_frequencies" : [ell_u, ell_l, ell_u, ell_l, ell],
                          "bessel_orders" : [1, 1, 0, 0, 2]}]
    elif estimator.lower() in ("cosebis", "log-cosebis"):
        if estimator.lower() == "cosebis":
            theta_mean = (theta_l + theta_u)/2
            Delta_theta = theta_u - theta_l
            B = Delta_theta/theta_mean/2
            if n == 1:
                X_1 = np.sqrt(8/5*(25 + 5*B**2 + 6*B**4))
                t_plus = lambda x: 1/X_1*(3*B**2 - 5 - 6*B*x + 3*(5-B**2)*x**2)
            elif n == 2:
                X_2 = np.sqrt(8*(25 + 5*B**2 + 6*B**4)*(175 + 35*B**2 + 45*B**4 + B**6))
                t_plus = lambda x: 1/X_2*(B**3*(25 + 3*B**2) - 15*(35 + 9*B**2 + 8*B**4)*x \
                                          -15*B**3*(3 + B**2)*x**2 + 35*(25 + 5*B**2 + 6*B**4)*x**3)
            else:
                t_plus = lambda x: np.sqrt((2*n+3)/2)*scipy.special.eval_legendre(n+1, x)
            T_plus = lambda theta: t_plus(2*(theta-theta_mean)/Delta_theta)
        else:
            #log-COSEBIS
            if log_cosebis_roots is None or log_cosebis_normalisations is None:
                # Compute the roots and normalisations
                log_cosebis_roots, log_cosebis_normalisations = cosebis.roots_normalization(n_max=n, z_max=np.log(theta_u/theta_l))

            r = log_cosebis_roots[n]
            N = log_cosebis_normalisations[n]
            t_plus = lambda z: N*np.prod(z - r)
            T_plus = lambda theta: t_plus(np.log(theta/theta_l))
        
        integral_info = [{"func" : lambda theta: theta*T_plus(theta)*J_0(ell*theta),
                          "bessel_frequencies" : [ell],
                          "bessel_orders" : [0]}]
    else:
        raise ValueError(f"Mode {estimator} not supported.")
        
    return sum([info.get("sign", 1)*integrate_bessel_functions(info["func"], theta_l, theta_u, 
                                                               bessel_frequencies=info["bessel_frequencies"], 
                                                               bessel_orders=info["bessel_orders"],
                                                               mode=integrator) for info in integral_info])
    
        