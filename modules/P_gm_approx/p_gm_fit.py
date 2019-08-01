import numpy as np

# evaluate a 2nd-order polynomial in 2 dimensions
def polyval_2d_order2(x, y, coeff):
    X, Y = np.meshgrid(x, y)
    return coeff[0] + coeff[1]*X + coeff[2]*Y + coeff[3]*X**2 + coeff[4]*X*Y + coeff[5]*Y**2


# main routine for fit functions
def prefac_nonlin(k, b2, g2, g3, omch2, h, ns, coeffs):
    logk = np.log(k)
    nonlin_term = np.zeros(shape=(3, len(k)))

    for j, ident in enumerate(['b2','g2','g3']):
        coeff = coeffs[ident]
        sum_poly_coeff = 0.0
        for i in range(3):
            poly_coeff_k = polyval_2d_order2(ns, omch2/h, coeff[i])
            sum_poly_coeff = sum_poly_coeff + poly_coeff_k * logk**(2-i)
        nonlin_term[j] = np.exp(sum_poly_coeff)  
    
    return np.where(k<1.e-2, 0.0, b2*nonlin_term[0] - g2*nonlin_term[1] - g3*nonlin_term[2])  # g2/3 F terms negative

def P_gm_approx(k, pofk_lin, pofk_nonlin, omch2, h, ns, b1, b2, g2, g3, coeffs):
    f = prefac_nonlin(k, b2, g2, g3, omch2, h, ns, coeffs)
    return b1*pofk_nonlin + f*pofk_lin**2