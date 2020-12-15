import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(Om0=0.265, H0=71.0)

mr = 26.5
z = np.linspace(10.0**(-12.0),3.0,101)

fluxlim = mr - 5.*np.log10(cosmo.luminosity_distance(z).value) - 25. -5.*np.log10(cosmo.H0.value/100.)

np.savetxt("/unix/atlas4/akorn/LSST/cosmosis/cosmosis/modules/euclid_ias/demos/without_sigma8_rescale_Aug_18/joint_new_halo_model/abs_mag_z.txt", np.c_[z, -26.0*np.ones(len(z)), fluxlim])
