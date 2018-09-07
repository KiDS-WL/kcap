MODULE constants

  IMPLICIT NONE

  !Mathematical constants
  REAL, PARAMETER :: pi=3.14159265359 !pi
  REAL, PARAMETER :: twopi=2.*pi !2pi or tau
  REAL, PARAMETER :: em=0.5772156649 !Euler–Mascheroni
  REAL, PARAMETER :: zero=0. !zero
  REAL, PARAMETER :: one=1. !one

  !Physical constants
  REAL, PARAMETER :: kb=1.38065e-23 !Boltzmann constant in SI  
  REAL, PARAMETER :: mp=1.6726219e-27 !Proton mass in kg
  REAL, PARAMETER :: bigG=6.67408e-11 !Gravitational constant [kg^-1 m^3 s^-2]
  REAL, PARAMETER :: eV=1.60218e-19 !Electronvolt in Joules
  REAL, PARAMETER :: cm=0.01 !Centimetre in metres
  REAL, PARAMETER :: rad2deg=180./pi !Radians-to-degrees conversion
  REAL, PARAMETER :: SBconst=5.670367e-8 !Steffan-Boltzmann constant [kg s^-3 K^-4]
  REAL, PARAMETER :: c_light=2.99792458e8 !Speed of light [m/s]
  
  !Cosmology constants
  REAL, PARAMETER :: Hdist=2997.9 !Hubble parameter distance (c/H0) [Mpc/h]
  REAL, PARAMETER :: Htime=9.7776 !Hubble time (1/H0) [Gyrs/h]
  REAL, PARAMETER :: H0=3.243e-18 !Hubble parameter [h/s]
  REAL, PARAMETER :: critical_density=2.7755e11 !Critical density at z=0 in (M_sun/h)/(Mpc/h)^3 (3*H0^2 / 8piG)
  REAL, PARAMETER :: dc0=(3./20.)*(12.*pi)**(2./3.) !Einstein-de Sitter linear collapse density ~1.686
  REAL, PARAMETER :: Dv0=18.*pi**2 !Einsten-de Sitter virialised collapse threshold ~178
  REAL, PARAMETER :: Msun=1.989e30 ! kg/Msun
  REAL, PARAMETER :: Mpc=3.086e22 ! m/Mpc
  REAL, PARAMETER :: yfac=8.125561e-16 !sigma_T/m_e*c^2 in SI

  !Weirdly specific things
  !REAL, PARAMETER :: fh=0.76 !Hydrogen mass fraction
  !REAL, PARAMETER :: mup=4./(5.*fh+3.) !Nuclear mass per particle (~0.588 if fh=0.76)
  !REAL, PARAMETER :: mue=2./(1.+fh) !Nuclear mass per electron (~1.136 if fh=0.76)
  !REAL, PARAMETER :: epfac=(5.*fh+3.)/(2.*(fh+1.)) !Electrons per total number of particles (P_th = P_e * epfac; ~1.932 if fh=0.76)
  
END MODULE constants
