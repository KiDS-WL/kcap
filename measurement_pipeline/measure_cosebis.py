import numpy as np
import math, os
from scipy.interpolate import interp1d
from scipy import pi,sqrt,exp
from scipy.special.orthogonal import p_roots
from numpy.polynomial.legendre import legcompanion, legval, legder
import numpy.linalg as la
from scipy.integrate import quad
# from decimal import *
# import sympy as smp
# import mpmath as mp

# getcontext().prec = 28
# mp.mp.dps=28


#  improve this
def leggauss(deg):
    """
    Gauss-Legendre quadrature.
    Computes the sample points and weights for Gauss-Legendre quadrature.
    These sample points and weights will correctly integrate polynomials of
    degree :math:`2*deg - 1` or less over the interval :math:`[-1, 1]` with
    the weight function :math:`f(x) = 1`.
    Parameters
    ----------
    deg : int
        Number of sample points and weights. It must be >= 1.
    Returns
    -------
    x : ndarray
        1-D ndarray containing the sample points.
    y : ndarray
        1-D ndarray containing the weights.
    Notes
    -----
    .. versionadded:: 1.7.0
    The results have only been tested up to degree 100, higher degrees may
    be problematic. The weights are determined by using the fact that
    .. math:: w_k = c / (L'_n(x_k) * L_{n-1}(x_k))
    where :math:`c` is a constant independent of :math:`k` and :math:`x_k`
    is the k'th root of :math:`L_n`, and then scaling the results to get
    the right value when integrating 1.
    """
    ideg = int(deg)
    if ideg != deg or ideg < 1:
        raise ValueError("deg must be a non-negative integer")

    # first approximation of roots. We use the fact that the companion
    # matrix is symmetric in this case in order to obtain better zeros.
    c = np.array([0]*deg + [1])
    m = legcompanion(c)
    x = la.eigvalsh(m)

    # improve roots by one application of Newton
    dy = legval(x, c)
    df = legval(x, legder(c))
    x -= dy/df

    # # improve roots again
    # dy = legval(x, c)
    # df = legval(x, legder(c))
    # x -= dy/df

    # compute the weights. We scale the factor to avoid possible numerical
    # overflow.
    fm = legval(x, c[1:])
    fm /= np.abs(fm).max()
    df /= np.abs(df).max()
    w = 1/(fm * df)

    # for Legendre we can also symmetrize
    w = (w + w[::-1])/2
    x = (x - x[::-1])/2

    # scale w to get the right value
    w *= 2. / w.sum()

    return x, w

#  improve this
def leggauss_improve(deg):
    """
    Gauss-Legendre quadrature.
    Computes the sample points and weights for Gauss-Legendre quadrature.
    These sample points and weights will correctly integrate polynomials of
    degree :math:`2*deg - 1` or less over the interval :math:`[-1, 1]` with
    the weight function :math:`f(x) = 1`.
    Parameters
    ----------
    deg : int
        Number of sample points and weights. It must be >= 1.
    Returns
    -------
    x : ndarray
        1-D ndarray containing the sample points.
    y : ndarray
        1-D ndarray containing the weights.
    Notes
    -----
    .. versionadded:: 1.7.0
    The results have only been tested up to degree 100, higher degrees may
    be problematic. The weights are determined by using the fact that
    .. math:: w_k = c / (L'_n(x_k) * L_{n-1}(x_k))
    where :math:`c` is a constant independent of :math:`k` and :math:`x_k`
    is the k'th root of :math:`L_n`, and then scaling the results to get
    the right value when integrating 1.
    """
    ideg = int(deg)
    if ideg != deg or ideg < 1:
        raise ValueError("deg must be a non-negative integer")

    # first approximation of roots. We use the fact that the companion
    # matrix is symmetric in this case in order to obtain better zeros.
    c = np.array([0]*deg + [1])
    m = legcompanion(c)
    x = la.eigvalsh(m)

    # # improve roots by one application of Newton
    dy = legval(x, c)
    df = legval(x, legder(c))
    x -= dy/df

    # # improve roots again
    # dy = legval(x, c)
    # df = legval(x, legder(c))
    # x -= dy/df

    # compute the weights. We scale the factor to avoid possible numerical
    # overflow.
    fm = legval(x, c[1:])
    fm /= np.abs(fm).max()
    df /= np.abs(df).max()
    w = 1/(fm * df)

    # for Legendre we can also symmetrize
    w = (w + w[::-1])/2
    x = (x - x[::-1])/2

    # scale w to get the right value
    w *= 2. / w.sum()

    return x, w, dy/df

# calculates T_plus logarithmic functions for COSEBIs
def tplus(tmin,tmax,n,norm,root,ntheta=10000):
    theta=np.logspace(np.log10(tmin),np.log10(tmax),ntheta)
    # 
    tplus=np.zeros((ntheta,2))
    tplus[:,0]=theta
    z=np.log(theta/tmin)
    result=1.
    for r in range(n+1):
        result*=(z-root[r])
# 
    result*=norm
    tplus[:,1]=result
    return tplus

# integrant for T_minus
def tminus_integ(y,z,tplus_func):
    # if y<0:
        # return 0
    return 4.*tplus_func(y)*(np.exp(2.*(y-z))-3.*np.exp(4.*(y-z)))


def gauss1(f,n):
    [x,w] = p_roots(n+1)
    G=sum(w*f(x))
    return G


def gauss(f,n,a,b):
    [x,w] = p_roots(n+1)
    G=0.5*(b-a)*sum(w*f(0.5*(b-a)*x+0.5*(b+a)))
    return G

# T_minus
def tminus(tmin,tmax,n,norm,root,tp,ntheta=10000,nG=20):
    tplus_func=interp1d(np.log(tp[:,0]/tmin),tp[:,1])
    theta=np.logspace(np.log10(tmin),np.log10(tmax),ntheta)
    # 
    tminus=np.zeros((ntheta,2))
    tminus[:,0]=theta
    z=np.log(theta/tmin)
    tminus[:,1]=tplus_func(z)
    [x,w] = p_roots(nG+1)
    integ_limits=np.insert(root/tmin,0,0)
    for iz in range(len(z)):
        result=0.
        good_integ=(integ_limits<=z[iz])
        integ_limits_good=integ_limits[good_integ]
        for il in range(1,len(integ_limits_good)):
            delta_limit=integ_limits_good[il]-integ_limits_good[il-1]
            y_in=0.5*delta_limit*x+0.5*(integ_limits_good[il]+integ_limits_good[il-1])
            y=y_in[y_in>=0.]
            result+=delta_limit*0.5*sum(w[y_in>=0.]*tminus_integ(y,z[iz],tplus_func))
        # print(il)
        delta_limit=z[iz]-integ_limits_good[-1]
        y_in=x*(delta_limit*0.5)+(z[iz]+integ_limits_good[-1])*0.5
        y=y_in[y_in>=0.]
        result+=delta_limit*0.5*sum(w[y_in>=0.]*tminus_integ(y,z[iz],tplus_func))
        tminus[iz,1]+=result
    return tminus

def tminus_quad(tmin,tmax,n,norm,root,tp,ntheta=10000):
    tplus_func=interp1d(np.log(tp[:,0]/tmin),tp[:,1])
    theta=np.logspace(np.log10(tmin),np.log10(tmax),ntheta)
    # 
    tminus=np.zeros((ntheta,2))
    tminus[:,0]=theta
    z=np.log(theta/tmin)
    tminus[:,1]=tplus_func(z)
    integ_limits=np.insert(root,0,0)
    # print(integ_limits)
    for iz in range(len(z)):
        good_integ=(integ_limits<=z[iz])
        integ_limits_good=integ_limits[good_integ]
        for il in range(1,len(integ_limits_good)):
            # print(il,z[iz],len(integ_limits_good))
            result=quad(tminus_integ,integ_limits[il-1] , integ_limits[il], args=(z[iz],tplus_func))
            tminus[iz,1]+=result[0]
        result=quad(tminus_integ,integ_limits[len(integ_limits_good)-1] ,z[iz], args=(z[iz],tplus_func))
        tminus[iz,1]+=result[0]
    return tminus


def tminus_longdouble(tmin,tmax,n,norm,root,tp,ntheta=10000,nG=20):
    tplus_func=interp1d(np.log(np.longdouble(tp[:,0])/np.longdouble(tmin)),np.longdouble(tp[:,1]))
    theta=np.logspace(np.log10(np.longdouble(tmin)),np.log10(np.longdouble(tmax)),ntheta)
    # 
    tminus=np.zeros((ntheta,2))
    tminus[:,0]=np.longdouble(theta)
    z=np.log(theta/np.longdouble(tmin))
    tminus[:,1]=tplus_func(z)
    [x,w] = p_roots(nG+1)
    integ_limits=np.insert(np.log(root/tmin),0,0)
    for iz in range(len(z)):
        result=0.
        good_integ=(integ_limits<=z[iz])
        integ_limits_good=integ_limits[good_integ]
        for il in range(1,len(integ_limits_good)):
            delta_limit=integ_limits_good[il]-integ_limits_good[il-1]
            y_in=0.5*delta_limit*x+0.5*(integ_limits_good[il]+integ_limits_good[il-1])
            y=y_in[y_in>=0.]
            result+=delta_limit*0.5*sum(w[y_in>=0.]*tminus_integ(y,z[iz],tplus_func))
        # print(il)
        delta_limit=z[iz]-integ_limits_good[-1]
        y_in=x*(delta_limit*0.5)+(z[iz]+integ_limits_good[-1])*0.5
        y=y_in[y_in>=0.]
        result+=delta_limit*0.5*sum(w[y_in>=0.]*tminus_integ(y,z[iz],tplus_func))
        tminus[iz,1]+=result
    return tminus

def integ_xi(xi_func,theta_edges, ntheta):
    ix=np.linspace(0,ntheta-1,ntheta)
    xip_integrated=np.zeros(len(theta_edges)-1)
    for tbin in range(len(theta_edges)-1):
        theta_in_range=np.exp(np.log(theta_edges[tbin])+(np.log(theta_edges[tbin+1])-np.log(theta_edges[tbin]))/(ntheta)*(ix+0.5))
        xip_integrated[tbin]=sum(xi_func(theta_in_range)*theta_in_range)/sum(theta_in_range)
    return xip_integrated


def integ_xi_noisy(xi_func,theta_edges, ntheta):
    ix=np.linspace(0,ntheta-1,ntheta)
    xip_integrated=np.zeros(len(theta_edges)-1)
    for tbin in range(len(theta_edges)-1):
        theta_in_range=np.exp(np.log(theta_edges[tbin])+(np.log(theta_edges[tbin+1])-np.log(theta_edges[tbin]))/(ntheta)*(ix+0.5))
        xip_integrated[tbin]=sum(xi_func(theta_in_range)*theta_in_range)/sum(theta_in_range)
    return xip_integrated
