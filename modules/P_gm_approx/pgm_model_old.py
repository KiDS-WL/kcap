# prototype code for fast approximate P_gm model

### settings
modelpath = "/share/splinter/joachimi/kids/kids1000/models/"
fitpath = "/share/splinter/joachimi/kids/kids1000/models/fitresults/"

omch2=0.1181838E+00
h=0.679565
ns=0.9679701
b1=1.8 #1.8871010
b2=-0.11710970
g2=0.2 #(-2.)/7.*(b1-1.)
g3=1.0 #0.79485430


### declarations
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate




def pofk_interpolator(pofk, k, z=None):
    if z is None:
        intp = scipy.interpolate.InterpolatedUnivariateSpline(np.log(k), np.log(pofk))
        return lambda k: np.exp(intp(np.log(k))).squeeze()
    else:
        intp = scipy.interpolate.RectBivariateSpline(z, np.log(k), np.log(pofk))
        return lambda k, z: np.exp(intp(z, np.log(k), grid=True)).squeeze()


# evaluate a 2nd-order polynomial in 2 dimensions
def polyval_2d_order2(x,y,coeff):
    X, Y = np.meshgrid(x, y)
    return coeff[0]+coeff[1]*X+coeff[2]*Y+coeff[3]*np.square(X)+coeff[4]*X*Y+coeff[5]*np.square(Y)


# main routine for fit functions
def prefac_nonlin(k,b2,g2,g3,omch2,h,ns):
    logk=np.log(k)
    nonlin_term=np.zeros(shape=(3,len(k)))

    for j,ident in enumerate(['b2','g2','g3']):
        coeff = np.loadtxt(fitpath+"/parameter_gridfit_"+ident+".dat")
        sum_poly_coeff=0.0
        for i in range(3):
            poly_coeff_k = polyval_2d_order2(ns,omch2/h,coeff[i])
            sum_poly_coeff = sum_poly_coeff + poly_coeff_k * np.power(logk,2-i)
        nonlin_term[j]=np.exp(sum_poly_coeff)  
    return np.where(k<1.e-2,0.0,b2*nonlin_term[0]-g2*nonlin_term[1]-g3*nonlin_term[2])  # g2/3 F terms negative



### main

## baseline power spectra (some only for comparison)
zeff = np.loadtxt(modelpath+"matter_matter_power_spectrum_pt/z.txt")
k = np.loadtxt(modelpath+"matter_matter_power_spectrum_pt/k_h.txt")
pk_mm = np.loadtxt(modelpath+"matter_matter_power_spectrum_pt/p_k.txt")
pk_gm = np.loadtxt(modelpath+"galaxy_matter_power_spectrum_pt/p_k.txt")
pk_gg = np.loadtxt(modelpath+"galaxy_galaxy_power_spectrum_pt/p_k.txt")

pk_lin_raw = np.loadtxt(modelpath+"matter_power_lin/p_k.txt")
pk_nl_raw = np.loadtxt(modelpath+"matter_power_nl/p_k.txt")
z_lin_raw = np.loadtxt(modelpath+"matter_power_lin/z.txt")
pk_lin = pofk_interpolator(pk_lin_raw,k,z_lin_raw)(k,zeff)
pk_nl = pofk_interpolator(pk_nl_raw,k,z_lin_raw)(k,zeff)



## get approximate power spectrum
# use pk_mm (PT matter PS) in first term to compare against PT results; use pk_nl (halofit) for KiDS modelling
pk_gm_approx_pt = b1*pk_mm + prefac_nonlin(k,b2,g2,g3,omch2,h,ns)*np.square(pk_lin)
pk_gm_approx = b1*pk_nl + prefac_nonlin(k,b2,g2,g3,omch2,h,ns)*np.square(pk_lin)



## compare power spectra in plot
zbin=0

plt.xlabel("k")
plt.ylabel("P(k)")
plt.grid()
plt.xlim(1.e-4,10.)
plt.ylim(80.,1.e5)
plt.xscale('log')
plt.yscale('log')

plt.plot(k,pk_gm[zbin],ls="-", label="PT")
plt.plot(k,pk_gm_approx_pt[zbin],ls="-", label="approx")
plt.plot(k,pk_gm_approx[zbin],ls="--", label="approx, halofit")

plt.legend(loc='lower left')
plt.savefig(fitpath+"/comparison.png")
#plt.show()

plt.clf()
plt.xlabel("k",fontsize=20)
plt.ylabel("P(k)/(b_1*P_{matter,pt}(k))",fontsize=20)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(1.e-3,5.)
#plt.ylim(0.9,1.8)
plt.ylim(0.0,1.1)
plt.xscale('log')

plt.plot(k,pk_gm[zbin]/(b1*pk_mm[zbin]),ls="-", label="PT, BOSS model")
plt.plot(k,pk_gm_approx_pt[zbin]/(b1*pk_mm[zbin]),ls="-", label="approx")
plt.plot(k,pk_gm_approx[zbin]/(b1*pk_mm[zbin]),ls="--", label="approx, halofit")
plt.plot(k,pk_nl[zbin]/pk_mm[zbin],ls="--", label="linear bias, halofit")

#plt.legend(loc='upper left',fontsize=20)
plt.savefig(fitpath+"/comparison_rel.png")
plt.show()
