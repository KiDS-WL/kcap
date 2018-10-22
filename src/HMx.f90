MODULE HMx

  ! Module usage statements
  USE constants
  USE array_operations
  USE solve_equations
  USE special_functions
  USE string_operations
  USE calculus_table
  USE cosmology_functions
  
  IMPLICIT NONE

  PRIVATE

  ! Main routines
  PUBLIC :: halomod ! Type
  PUBLIC :: assign_halomod
  PUBLIC :: init_halomod
  PUBLIC :: print_halomod
  PUBLIC :: calculate_HMx
  PUBLIC :: calculate_HMx_a
  PUBLIC :: calculate_HMcode_a
  PUBLIC :: set_halo_type
  PUBLIC :: halo_type
  PUBLIC :: M_nu
  PUBLIC :: virial_radius
  PUBLIC :: convert_mass_definitions
  PUBLIC :: win_type
  PUBLIC :: UPP ! TODO: Retire
  PUBLIC :: p_1void ! TODO: Retire
  PUBLIC :: halo_HI_fraction ! TODO: Retire

  ! Diagnostics
  PUBLIC :: halo_definitions
  PUBLIC :: halo_diagnostics
  PUBLIC :: halo_properties
  PUBLIC :: write_halo_profiles
  PUBLIC :: write_mass_fractions

  ! Winint functions
  PUBLIC :: winint_diagnostics
  PUBLIC :: winint_speed_tests

  ! Mass functions and bias
  ! TODO: Probably most of these should be private
  PUBLIC :: mass_function
  PUBLIC :: b_ps
  PUBLIC :: b_st
  PUBLIC :: b_Tinker
  PUBLIC :: g_ps
  PUBLIC :: g_st
  PUBLIC :: g_Tinker

  ! HMx functions
  PUBLIC :: HMx_alpha
  PUBLIC :: HMx_eps
  PUBLIC :: HMx_Gamma
  PUBLIC :: HMx_M0
  PUBLIC :: HMx_Astar
  PUBLIC :: HMx_Twhim

  ! Fields
  PUBLIC :: n_fields
  PUBLIC :: field_dmonly
  PUBLIC :: field_matter
  PUBLIC :: field_cdm
  PUBLIC :: field_gas
  PUBLIC :: field_star
  PUBLIC :: field_bound_gas
  PUBLIC :: field_free_gas
  PUBLIC :: field_electron_pressure
  PUBLIC :: field_void
  PUBLIC :: field_compensated_void
  PUBLIC :: field_central_galaxies
  PUBLIC :: field_satellite_galaxies
  PUBLIC :: field_galaxies
  PUBLIC :: field_HI
  PUBLIC :: field_cold_gas

  ! Halo-model stuff that needs to be recalculated for each new z
  TYPE halomod
     INTEGER :: ip2h, ibias, imf, iconc, iDolag, iAs, ip2h_corr, ikb
     INTEGER :: idc, iDv, ieta, ikstar, i2hdamp, i1hdamp, itrans
     LOGICAL :: voids
     REAL :: z, a, dc, Dv
     REAL :: alpha, eps, Gamma, M0, Astar, Twhim, cstar, sstar, mstar, Theat, fcold, alphap, Gammap, cstarp ! HMx baryon parameters
     REAL :: A_alpha, B_alpha, C_alpha, D_alpha, E_alpha
     REAL :: A_eps, B_eps, C_eps, D_eps
     REAL :: A_Gamma, B_Gamma, C_Gamma, D_Gamma, E_gamma
     REAL :: A_M0, B_M0, C_M0, D_M0, E_M0
     REAL :: A_Astar, B_Astar, C_Astar, D_Astar
     REAL :: A_Twhim, B_Twhim, C_Twhim, D_Twhim
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:), log_m(:)
     REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
     REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
     REAL, ALLOCATABLE :: k(:), wk(:,:,:)
     INTEGER :: nk
     REAL :: sigv, sigv100, c3, knl, rnl, mnl, neff, sig8z
     REAL :: gmin, gmax, gbmin, gbmax
     REAL :: n_c, n_s, n_g, rho_HI, dlnc
     REAL :: Dv0, Dv1, dc0, dc1, eta0, eta1, f0, f1, ks, As, alp0, alp1 ! HMcode parameters
     REAL :: mgal, HImin, HImax, rcore, hmass
     INTEGER :: n
     INTEGER :: halo_DMONLY, halo_CDM, halo_boundgas, halo_coldgas, halo_freegas, halo_star, halo_HI
     INTEGER :: halo_void, halo_compensated_void, electron_pressure
     INTEGER :: frac_boundgas, frac_star, frac_HI, frac_coldgas
     LOGICAL :: one_parameter_baryons
     LOGICAL :: has_HI, has_galaxies, has_mass_conversions, safe_negative, has_dewiggle
     LOGICAL :: fixed_HMx, response
     REAL :: acc_HMx, large_nu
     CHARACTER(len=256) :: name
     REAL, ALLOCATABLE :: log_k_pdamp(:), log_pdamp(:)
     INTEGER :: n_pdamp
  END TYPE halomod

  ! Window integration
  REAL, PARAMETER :: acc_win=1e-3           ! Window-function integration accuracy parameter
  INTEGER, PARAMETER :: imeth_win=3         ! Window-function integration method
  INTEGER, PARAMETER :: winint_order=3      ! Window-function integration order
  REAL, PARAMETER :: winint_test_seconds=2. ! Approximately how many seconds should each timing test take
    
  ! Halomodel
  LOGICAL, PARAMETER :: slow_hmod=.FALSE.            ! Choose to do the slower hmod initialisation (unnecessary calculations)  
  REAL, PARAMETER :: mmin_response=1e7               ! For HMcode if doing a response
  REAL, PARAMETER :: mmax_response=1e17              ! For HMcode if doing a response
  LOGICAL, PARAMETER :: verbose_convert_mass=.FALSE. ! Verbosity when running the mass-definition-conversion routines
  REAL, PARAMETER :: zinf_Dolag=100.                 ! An approximate infinite z; I changed this from 10 -> 100 to make it more infinite, I used zinf=10 for Mead (2017)!

  ! HMcode
  REAL, PARAMETER :: mmin_hmcode=1e7          ! Minimum mass to consider for one-halo integration
  REAL, PARAMETER :: mmax_hmcode=1e17         ! Maximum mass to consider for one-halo integration
  REAL, PARAMETER :: fdamp_min=1e-3           ! Minimum value for f_damp parameter
  REAL, PARAMETER :: fdamp_max=0.99           ! Maximum value for f_damp parameter
  REAL, PARAMETER :: alpha_transition_min=0.5 ! Minimum value for alpha transition parameter
  REAL, PARAMETER :: alpha_transition_max=2.0 ! Maximum value for alpha transition parameter

  ! HMx
  REAL, PARAMETER :: HMx_alpha_min=1e-2 ! Minimum alpha parameter; needs to be set at not zero
  REAL, PARAMETER :: HMx_Gamma_min=1.05 ! Minimum polytropic index
  REAL, PARAMETER :: HMx_Gamma_max=2.00 ! Maximum polytropic index
  REAL, PARAMETER :: HMx_Astar_min=1e-4 ! Minimum halo star fraction; needs to be set at not zero

  ! Halo types
  LOGICAL, PARAMETER :: verbose_galaxies=.FALSE. ! Verbosity when doing the galaxies initialisation
  LOGICAL, PARAMETER :: verbose_HI=.TRUE.       ! Verbosity when doing the HI initialisation

  ! Field types
  INTEGER, PARAMETER :: n_fields=15
  INTEGER, PARAMETER :: field_dmonly=-1
  INTEGER, PARAMETER :: field_matter=0
  INTEGER, PARAMETER :: field_cdm=1
  INTEGER, PARAMETER :: field_gas=2
  INTEGER, PARAMETER :: field_star=3
  INTEGER, PARAMETER :: field_bound_gas=4
  INTEGER, PARAMETER :: field_free_gas=5
  INTEGER, PARAMETER :: field_electron_pressure=6
  INTEGER, PARAMETER :: field_void=7
  INTEGER, PARAMETER :: field_compensated_void=8
  INTEGER, PARAMETER :: field_central_galaxies=9
  INTEGER, PARAMETER :: field_satellite_galaxies=10
  INTEGER, PARAMETER :: field_galaxies=11
  INTEGER, PARAMETER :: field_HI=12
  INTEGER, PARAMETER :: field_cold_gas=13

CONTAINS

   SUBROUTINE assign_halomod(ihm,hmod,verbose)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ihm
    TYPE(halomod), INTENT(OUT) :: hmod
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    ! Names of pre-defined halo models
    INTEGER, PARAMETER :: nhalomod=27 ! Total number of pre-defined halo-model types (TODO: this is stupid)
    CHARACTER(len=256):: names(1:nhalomod)    
    names(1)='Accurate HMcode (Mead et al. 2016)'
    names(2)='Basic halo-model (Two-halo term is linear)'
    names(3)='Standard halo-model (Seljak 2000)'
    names(4)='Standard halo-model but with Mead et al. (2015) transition'
    names(5)='Standard halo-model but with Delta_v=200 and delta_c=1.686 and Bullock c(M)'
    names(6)='Half-accurate HMcode (Mead et al. 2015, 2016)'
    names(7)='Accurate HMcode (Mead et al. 2015)'
    names(8)='Including scatter in halo properties at fixed mass'
    names(9)='Parameters for CCL tests (high accuracy)'
    names(10)='Comparison of mass conversions with Wayne Hu code'
    names(11)='Standard halo-model (Seljak 2000) but with UPP for electron pressure'
    names(12)='Spherical collapse used for Mead (2017) results'
    names(13)='Experimental log-tanh transition'
    names(14)='Experimental scale-dependent halo bias'
    names(15)='Accurate HMcode (Mead et al. 2018)'
    names(16)='Halo-void model'
    names(17)='HMx - AGN 7.6'
    names(18)='HMx - AGN 7.8'
    names(19)='HMx - AGN 8.0'
    names(20)='Standard halo-model (Seljak 2000) in response'
    names(21)='Cored profile model'
    names(22)='Delta function-NFW star profile model response'
    names(23)='Tinker mass function and bias'
    names(24)='Universal pressure profile'
    names(25)='Villaescusa-Navarro HI halo model'
    names(26)='Delta-function mass function'
    names(27)='Press & Schecter mass function'

    IF(verbose) WRITE(*,*) 'ASSIGN_HALOMOD: Assigning halo model'
    
    ! Default options

    ! Number of points in integration (128 is okay, 1024 is better)
    hmod%n=128

    ! Accuracy for continuous integrals (1e-3 is okay, 1e-4 is better)
    hmod%acc_HMx=1e-3

    ! A large value for nu (6 is okay, corrections are suppressed by exp(-large_nu^2)
    hmod%large_nu=6.

    ! Two-halo term
    ! 1 - Linear theory
    ! 2 - Standard from Seljak (2000)
    ! 3 - Linear theory with damped wiggles
    ! 4 - No two-halo term
    hmod%ip2h=2

    ! Method to correct the two-halo integral
    ! NB. This cannot be a parameter here because the value needs to be changed if doing cumulative distributions of power with mass
    ! 1 - Do nothing
    ! 2 - Add value of missing integral assuming that W(k)=1
    ! 3 - Put the missing part of the integrand as a delta function at lower mass limit
    hmod%ip2h_corr=3

    ! Scale dependent halo bias
    ! 1 - None
    ! 2 - From Fedeli (2014b)
    ! 3 - My own experimental model
    hmod%ikb=1

    ! Order of halo bias to go to
    ! 1 - Linear order (standard)
    ! 2 - Second order
    hmod%ibias=1

    ! One-halo term large-scale damping
    ! 1 - No damping
    ! 2 - Mead et al. (2015)
    ! 3 - k^4 at large scales
    hmod%i1hdamp=1

    ! Mass and halo bias function pair
    ! 1 - Press & Schecter (1974)
    ! 2 - Sheth & Tormen (1999)
    ! 3 - Tinker et al. (2010)
    ! 4 - Delta function in mass
    hmod%imf=2

    ! Concentration-mass relation
    ! 1 - Full Bullock et al. (2001; astro-ph/9909159)
    ! 2 - Simple Bullock et al. (2001; astro-ph/9909159)
    ! 3 - Duffy et al. (2008; astro-ph/0804.2486): mean
    ! 4 - Duffy et al. (2008; astro-ph/0804.2486): virial
    ! 5 - Duffy et al. (2008; astro-ph/0804.2486): relaxed-200
    hmod%iconc=4

    ! Linear collapse threshold delta_c
    ! 1 - Fixed 1.686
    ! 2 - Nakamura & Suto (1997) fitting function
    ! 3 - Mead et al. (2015)
    ! 4 - Mead (2017) fitting function
    ! 5 - Spherical-collapse calculation
    hmod%idc=2

    ! Virial density Delta_v
    ! 1 - Fixed 200
    ! 2 - Bryan & Norman (1998; arXiv:astro-ph/9710107) fitting function
    ! 3 - Mead et al. (2015)
    ! 4 - Mead (2017) fitting function
    ! 5 - Spherical-collapse calculation
    ! 6 - Fixed to unity to give Lagrangian radius
    hmod%iDv=2

    ! eta for halo window function
    ! 1 - No
    ! 2 - Mead et al. (2015)
    hmod%ieta=1

    ! k* for one-halo term large-scale damping
    ! 1 - No
    ! 2 - Mead et al. (2015)
    hmod%ikstar=1

    ! Concentration-mass rescaling
    ! 1 - No
    ! 2 - Mead et al. (2015, 2016)
    hmod%iAs=1

    ! fdamp for two-halo term damping
    ! 1 - No
    ! 2 - Mead et al. (2015)
    ! 3 - Mead et al. (2016)
    hmod%i2hdamp=1

    ! alpha for two- to one-halo transition region
    ! 1 - No
    ! 2 - Smoothed transition with alphas
    ! 4 - New HMx transition
    ! 5 - Tanh transition
    hmod%itrans=1

    ! Use the Dolag c(M) correction for dark energy?
    ! 1 - No
    ! 2 - Yes, exactly as in Dolag et al. (2004)
    ! 3 - Yes, as in Dolag et al. (2004) but with a ^1.5 power
    hmod%iDolag=2

    ! Halo gas fraction
    ! 1 - Fedeli (2014a) bound gas model
    ! 2 - Schneider (2015) bound gas
    ! 3 - Universal baryon fraction
    hmod%frac_boundgas=2

    ! Halo cold gas fraction
    ! 1 - Constant fraction of halo gas
    hmod%frac_coldgas=1
    
    ! Halo star fraction
    ! 1 - Fedeli (2014)
    ! 2 - Constant stellar fraction
    ! 3 - Fedeli (2014) but saturates at high halo mass
    ! 4 - No stars
    hmod%frac_star=3

    ! Halo HI fraction
    ! 1 - Simple model
    ! 2 - Villaescusa-Navarro et al. (2018; 1804.09180)
    hmod%frac_HI=1

    ! DMONLY halo profile
    ! 1 - Analyical NFW
    ! 2 - Non-analytical NFW (good for testing W(k) functions)
    ! 3 - Tophat
    ! 4 - Delta function
    ! 5 - Cored NFW
    ! 6 - Isothermal
    ! 7 - Shell
    hmod%halo_DMONLY=1

    ! CDM halo profile
    ! 1 - NFW
    hmod%halo_CDM=1

    ! Bound gas halo profile
    ! 1 - Simplified Komatsu & Seljak (2001) gas model
    ! 2 - Isothermal beta model
    ! 3 - Full Komatsu & Seljak (2001) gas model
    hmod%halo_boundgas=1

    ! Cold gas halo profile
    ! 1 - Delta function
    hmod%halo_coldgas=1

    ! Free gas halo profile
    ! 1 - Isothermal model (out to 2rv)
    ! 2 - Ejected gas model from Schneider (2015)
    ! 3 - Isothermal shell that connects electron pressure and density to boundgas at rv
    ! 4 - Komatsu-Seljak continuation
    ! 5 - Power-law continuation
    ! 6 - Cubic profile
    ! 7 - Smoothly distributed
    ! 8 - Delta function
    hmod%halo_freegas=7

    ! Stars halo profile
    ! 1 - Fedeli (2014) stellar distribution
    ! 2 - Schneider (2015) stellar distribution
    ! 3 - Delta function
    ! 4 - Delta function at low mass and NFW at high mass
    hmod%halo_star=1

    ! HI halo profile
    ! 1 - NFW
    ! 2 - Delta function
    ! 3 - Polynomial with internal exponential hole
    ! 4 - NFW with internal exponential hole
    ! 5 - Modified NFW
    hmod%halo_HI=1

    ! Electron pressure
    ! 1 - UPP (Arnaud et al. 2010)
    ! 2 - Bound and free gas profiles
    hmod%electron_pressure=2

    ! Do voids?
    hmod%voids=.FALSE.

    ! Set the void model
    ! 1 - Top-hat void
    hmod%halo_void=1

    ! Set the compensated void model
    ! 1 - Top-hat void
    hmod%halo_compensated_void=1

    ! Safeguard against negative terms in cross correlations
    hmod%safe_negative=.FALSE.

    ! HMcode parameters
    hmod%Dv0=418.
    hmod%Dv1=-0.352
    hmod%dc0=1.59
    hmod%dc1=0.0314
    hmod%eta0=0.603
    hmod%eta1=0.300
    hmod%f0=0.0095
    hmod%f1=1.37
    hmod%ks=0.584
    hmod%As=3.13
    hmod%alp0=3.24
    hmod%alp1=1.85
    hmod%one_parameter_baryons=.FALSE.

    ! HMx parameters
    hmod%fixed_HMx=.TRUE.
    hmod%alpha=0.33333 ! Non-virial temperature thing
    hmod%eps=1.        ! Concentration modification
    hmod%Gamma=1.17    ! Polytropic gas index
    hmod%M0=1e14       ! Halo mass that has lost half gas
    hmod%Astar=0.03    ! Maximum star-formation efficiency
    hmod%Twhim=1e6     ! WHIM temperature [K]
    hmod%cstar=10.     ! Stellar concentration r_* = rv/c
    hmod%sstar=1.2     ! sigma_* for f_* distribution
    hmod%Mstar=5e12    ! M* for most efficient halo mass for star formation
    hmod%fcold=0.0     ! Fraction of cold gas, in addition to bound gas
    hmod%alphap=0.0    ! Power-law index of alpha with halo mass
    hmod%Gammap=0.0    ! Power-law index of Gamma with halo mass
    hmod%cstarp=0.0    ! Power-law index of c* with halo mass

    ! $\alpha$ z and Theat variation
    hmod%A_alpha=-0.005
    hmod%B_alpha=0.022
    hmod%C_alpha=0.865
    hmod%D_alpha=-13.565
    hmod%E_alpha=52.516

    ! $\log_{10} \epsilon_c$ z and Theat variation
    hmod%A_eps=-0.289
    hmod%B_eps=2.147
    hmod%C_eps=0.129
    hmod%D_eps=-0.867

    ! $\Gamma$ z and Theat variation
    hmod%A_Gamma=0.026
    hmod%B_Gamma=-0.064
    hmod%C_Gamma=1.150
    hmod%D_Gamma=-17.011
    hmod%E_Gamma=66.289

    ! $\log_{10}M_0$ z and Theat variation
    hmod%A_M0=-0.007
    hmod%B_M0=0.018
    hmod%C_M0=6.788
    hmod%D_M0=-103.453
    hmod%E_M0=406.705

    ! $A_*$ z and Theat variation
    hmod%A_Astar=0.004
    hmod%B_Astar=-0.034
    hmod%C_Astar=-0.017
    hmod%D_Astar=0.165

    ! $\log_{10}T_{WHIM}$ z and Theat variation
    hmod%A_Twhim=-0.024
    hmod%B_Twhim=0.077
    hmod%C_Twhim=0.454
    hmod%D_Twhim=1.717

    ! Do we treat the halomodel as a response model (multiply by HMcode) or not
    hmod%response=.FALSE.

    ! Halo mass if the mass function is a delta function
    hmod%hmass=1e13

    ! Scatter in ln(c)
    hmod%dlnc=0.

    ! HOD parameters
    hmod%mgal=1e13

    ! HI parameters
    hmod%HImin=1e9
    hmod%HImax=1e12

    ! NFW core radius
    hmod%rcore=0.1

    IF(ihm==-1) THEN
       WRITE(*,*) 'ASSIGN_HALOMOD: Choose your halo model'
       DO i=1,nhalomod
          WRITE(*,*) i, TRIM(names(i))
       END DO
       READ(*,*) ihm
       WRITE(*,*)
    END IF
       
    IF(ihm==1 .OR. ihm==7 .OR. ihm==15) THEN
       !  1 - Accurate HMcode  (Mead et al. 2016)
       !  7 - Accurate HMcode  (Mead et al. 2015)
       ! 15 - Accurate HMcode  (Mead et al. 2018)
       hmod%ip2h=1
       hmod%i1hdamp=2
       hmod%iconc=1
       hmod%idc=3
       hmod%iDv=3
       hmod%ieta=2
       hmod%ikstar=2
       hmod%iAs=2
       hmod%i2hdamp=3
       hmod%itrans=2
       hmod%iDolag=3
       IF(ihm==7) THEN
          ! Mead et al. (2015)
          hmod%i2hdamp=2
          hmod%itrans=2
          hmod%f0=0.188
          hmod%f1=4.29
          hmod%alp0=2.93
          hmod%alp1=1.77
          hmod%iDolag=2
       ELSE IF(ihm==15) THEN
          ! Mead et al. (2018)
          hmod%i1hdamp=3 ! k^4 at large scales for one-halo term
          hmod%ip2h=3    ! Linear theory with damped wiggles
          hmod%i2hdamp=2 ! Change back to Mead (2015) model for two-halo damping
          hmod%Dv0=446.8
          hmod%Dv1=-0.3237
          hmod%dc0=1.6145
          hmod%dc1=0.007774
          hmod%eta0=0.5403
          hmod%eta1=0.2345
          hmod%f0=0.09563
          hmod%f1=4.1524
          hmod%ks=0.6902
          hmod%As=3.020
          hmod%alp0=3.072
          hmod%alp1=1.848
       END IF       
    ELSE IF(ihm==2) THEN
       ! Basic halo model with linear two halo term (Delta_v = 200, delta_c = 1.686))
       hmod%ip2h=1
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
    ELSE IF(ihm==3) THEN
       ! Standard halo-model calculation (Seljak 2000)
       ! This is the default, so do nothing here
    ELSE IF(ihm==4) THEN
       ! Standard halo-model calculation but with Mead et al. (2015) smoothed two- to one-halo transition and one-halo damping       
       hmod%itrans=4
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%safe_negative=.TRUE.
    ELSE IF(ihm==5) THEN
       ! Standard halo-model calculation but with Delta_v = 200 and delta_c = 1.686 fixed and Bullock c(M)
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
    ELSE IF(ihm==6) THEN
       ! Half-accurate halo-model calculation, inspired by (Mead et al. 2015, 2016)
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%i2hdamp=3
       hmod%safe_negative=.TRUE.
       hmod%idc=3
       hmod%iDv=3
       hmod%ieta=2
       hmod%iAs=2
       hmod%itrans=2
       hmod%iconc=1
       hmod%iDolag=3
    ELSE IF(ihm==8) THEN
       ! Include scatter in halo properties
       hmod%dlnc=0.25
    ELSE IF(ihm==9) THEN
       ! For CCL comparison and benchmark generation
       hmod%n=2048 ! Increase accuracy for the CCL benchmarks
       hmod%acc_HMx=1e-5 ! Increase accuracy for the CCL benchmarks
       hmod%ip2h=2
       hmod%ip2h_corr=3
       hmod%ibias=1
       hmod%i1hdamp=1
       hmod%imf=2
       hmod%iconc=4 ! Virial Duffy relation
       hmod%idc=2   ! Virial dc
       hmod%iDv=2   ! Virial Dv
       hmod%ieta=1
       hmod%ikstar=1
       hmod%iAs=1
       hmod%i2hdamp=1
       hmod%itrans=1
       hmod%iDolag=1
    ELSE IF(ihm==10) THEN
       ! For mass conversions comparison with Wayne Hu's code
       hmod%iconc=2
       hmod%idc=1
       hmod%iDv=1
    ELSE IF(ihm==11) THEN
       ! UPP
       hmod%electron_pressure=1
    ELSE IF(ihm==12) THEN
       ! Spherical-collapse model to produce Mead (2017) results
       hmod%iconc=1
       hmod%idc=5
       hmod%iDv=5
       hmod%iDolag=2 ! This seems not to be important for these results
    ELSE IF(ihm==13) THEN
       ! Experimental log-tanh transition
       hmod%itrans=5
    ELSE IF(ihm==14) THEN
       ! Experimental scale-dependent halo bias
       hmod%ikb=3
       hmod%ikstar=2
       hmod%i1hdamp=3
    ELSE IF(ihm==16) THEN
       ! Halo-void model
       hmod%voids=.TRUE.
    ELSE IF(ihm==17 .OR. ihm==18 .OR. ihm==19) THEN
       ! 17 - HMx AGN 7.6
       ! 18 - HMx AGN 7.8
       ! 19 - HMx AGN 8.0
       hmod%fixed_HMx=.FALSE.
       hmod%itrans=4
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%safe_negative=.TRUE.
       hmod%response=.TRUE.
       IF(ihm==17) THEN
          ! AGN 7.6
          hmod%Theat=10**7.6
       ELSE IF(ihm==18) THEN
          ! AGN 7.8
          hmod%Theat=10**7.8
       ELSE IF(ihm==19) THEN
          ! AGN 8.0
          hmod%Theat=10**8.0
       END IF
    ELSE IF(ihm==20) THEN
       ! Standard halo model but as response with HMcode
       hmod%response=.TRUE.
    ELSE IF(ihm==21) THEN
       ! Cored NFW halo profile model
       hmod%halo_DMONLY=5 ! Cored profile
    ELSE IF(ihm==22) THEN
       ! Different stellar profile
       hmod%halo_star=2 ! Schneider & Teyssier (2015)
       hmod%response=.TRUE.
    ELSE IF(ihm==23) THEN
       ! Tinker mass function and bias
       hmod%imf=3
    ELSE IF(ihm==24) THEN
       ! UPP for electron pressure
       hmod%electron_pressure=1
    ELSE IF(ihm==25) THEN
       ! Villaescusa-Navarro HI halo model
       hmod%imf=3     ! Tinker mass function
       hmod%frac_HI=3 ! HI mass fraction with z evolution (Villaescusa-Navarro et al. 1804.09180)
       !hmod%halo_HI=3 ! Exponentially cored polynomial (Villaescusa-Navarro et al. 1804.09180)
       !hmod%halo_HI=2 ! Delta function
       hmod%halo_HI=5 ! Modified NFW (Padmanabhan & Refregier 2017; 1607.01021)
    ELSE IF(ihm==26) THEN
       ! Delta function mass function
       hmod%ip2h=2        ! Standard two-halo term
       hmod%iDv=6         ! Lagrangian radius haloes
       hmod%imf=4         ! Delta function mass function
       hmod%halo_DMONLY=3 ! Top-hat halo profile
    ELSE IF(ihm==27) THEN
       hmod%imf=1 ! Press & Schecter (1974) mass function
    ELSE
       STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
    END IF
    hmod%name=names(ihm)

    IF(verbose) THEN
       WRITE(*,*) 'ASSIGN_HALOMOD: ', TRIM(names(ihm))
       WRITE(*,*) 'ASSIGN_HALOMOD: Done'
       WRITE(*,*)
    END IF
    
  END SUBROUTINE assign_halomod
  
  SUBROUTINE init_halomod(mmin,mmax,a,hmod,cosm,verbose)

    ! Halo-model initialisation routine
    ! The computes other tables necessary for the one-halo integral
    IMPLICIT NONE
    REAL, INTENT(IN) :: mmin, mmax
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i
    REAL :: Dv, dc, m, nu, R, sig, A0, z  

    ! Set the redshift (this routine needs to be called anew for each z)
    hmod%a=a
    z=redshift_a(a)
    hmod%z=z

    IF(ALLOCATED(hmod%log_m)) CALL deallocate_HMOD(hmod)
    CALL allocate_HMOD(hmod)

    ! Set flags to false
    hmod%has_galaxies=.FALSE.
    hmod%has_HI=.FALSE.
    hmod%has_mass_conversions=.FALSE.
    hmod%has_dewiggle=.FALSE.

    ! Find value of sigma_V
    hmod%sigv=sigmaV(0.,a,cosm)
    IF(hmod%i2hdamp==3) hmod%sigv100=sigmaV(100.,a,cosm)
    hmod%sig8z=sigma(8.,a,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Filling look-up tables'
       WRITE(*,*) 'INIT_HALOMOD: Number of entries:', hmod%n
       WRITE(*,*) 'INIT_HALOMOD: Tables being filled at redshift:', REAL(z)
       WRITE(*,*) 'INIT_HALOMOD: Tables being filled at scale-factor:', REAL(a)
       WRITE(*,*) 'INIT_HALOMOD: sigma_V [Mpc/h]:', REAL(hmod%sigv)
       IF(hmod%i2hdamp==3) WRITE(*,*) 'INIT_HALOMOD: sigmaV_100 [Mpc/h]:', REAL(hmod%sigv100)
       WRITE(*,*) 'INIT_HALOMOD: sigma_8(z):', REAL(hmod%sig8z)
    END IF

    ! Loop over halo masses and fill arrays
    DO i=1,hmod%n

       m=exp(progression(log(mmin),log(mmax),i,hmod%n))
       R=radius_m(m,cosm)
       sig=sigma(R,a,cosm)
       nu=nu_R(R,hmod,cosm)

       hmod%m(i)=m
       hmod%rr(i)=R
       hmod%sig(i)=sig
       hmod%nu(i)=nu

    END DO

    ! log-mass array
    hmod%log_m=log(hmod%m)

    IF(verbose) WRITE(*,*) 'INIT_HALOMOD: M, R, nu, sigma tables filled'

    ! Get delta_c
    dc=delta_c(hmod,cosm)
    hmod%dc=dc

    ! Fill virial radius table using real radius table
    Dv=Delta_v(hmod,cosm)
    hmod%Dv=Dv
    hmod%rv=hmod%rr/(Dv**(1./3.))

    ! Write some useful information to the screen
    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: virial radius tables filled'
       WRITE(*,*) 'INIT_HALOMOD: Delta_v:', REAL(Dv)
       WRITE(*,*) 'INIT_HALOMOD: delta_c:', REAL(dc)
       WRITE(*,*) 'INIT_HALOMOD: Minimum nu:', REAL(hmod%nu(1))
       WRITE(*,*) 'INIT_HALOMOD: Maximum nu:', REAL(hmod%nu(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Minimum R_v [Mpc/h]:', REAL(hmod%rv(1))
       WRITE(*,*) 'INIT_HALOMOD: Maximum R_v [Mpc/h]:', REAL(hmod%rv(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Minimum log10(M) [Msun/h]:', REAL(log10(hmod%m(1)))
       WRITE(*,*) 'INIT_HALOMOD: Maximum log10(M) [Msun/h]:', REAL(log10(hmod%m(hmod%n)))
    END IF
    
    IF(hmod%imf==4) THEN

       ! Do nothing
       
    ELSE

       ! Calculate missing mass things if necessary

       IF(hmod%ip2h_corr==2 .OR. hmod%ip2h_corr==3) THEN

          IF(slow_hmod) hmod%gmin=1.-integrate_hmod(hmod%nu(1),hmod%large_nu,g_nu,hmod,hmod%acc_HMx,3)
          IF(slow_hmod) hmod%gmax=integrate_hmod(hmod%nu(hmod%n),hmod%large_nu,g_nu,hmod,hmod%acc_HMx,3)
          hmod%gbmin=1.-integrate_hmod(hmod%nu(1),hmod%large_nu,gb_nu,hmod,hmod%acc_HMx,3)
          IF(slow_hmod) hmod%gbmax=integrate_hmod(hmod%nu(hmod%n),hmod%large_nu,gb_nu,hmod,hmod%acc_HMx,3)

          IF(verbose) THEN          
             IF(slow_hmod) WRITE(*,*) 'INIT_HALOMOD: Missing g(nu) at low end:', REAL(hmod%gmin)
             IF(slow_hmod) WRITE(*,*) 'INIT_HALOMOD: Missing g(nu) at high end:', REAL(hmod%gmax)
             WRITE(*,*) 'INIT_HALOMOD: Missing g(nu)b(nu) at low end:', REAL(hmod%gbmin)
             IF(slow_hmod) WRITE(*,*) 'INIT_HALOMOD: Missing g(nu)b(nu) at high end:', REAL(hmod%gbmax)          
          END IF

          IF(slow_hmod .AND. hmod%gmin<0.) STOP 'INIT_HALOMOD: Error, missing g(nu) at low end is less than zero'       
          IF(slow_hmod .AND. hmod%gmax<-1e-4) STOP 'INIT_HALOMOD: Error, missing g(nu) at high end is less than zero'
          IF(hmod%gbmin<0.) STOP 'INIT_HALOMOD: Error, missing g(nu)b(nu) at low end is less than zero'
          IF(slow_hmod .AND. hmod%gbmax<-1e-4) STOP 'INIT_HALOMOD: Error, missing g(nu)b(nu) at high end is less than zero'       

       END IF

    END IF
       
    ! Calculate the total stellar mass fraction
    IF(slow_hmod .AND. verbose) WRITE(*,*) 'INIT_HALOMOD: Omega_stars:', Omega_stars(hmod,cosm)

    ! Find non-linear radius and scale
    ! This is defined as nu(M_star)=1 *not* sigma(M_star)=1, so depends on delta_c
    hmod%rnl=r_nl(hmod)
    hmod%mnl=mass_r(hmod%rnl,cosm)
    hmod%knl=1./hmod%rnl

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Non-linear mass [log10(M*) [Msun/h]]:', REAL(log10(hmod%mnl))
       WRITE(*,*) 'INIT_HALOMOD: Non-linear halo virial radius [Mpc/h]:', REAL(virial_radius(hmod%mnl,hmod,cosm))
       WRITE(*,*) 'INIT_HALOMOD: Non-linear Lagrangian radius [Mpc/h]:', REAL(hmod%rnl)
       WRITE(*,*) 'INIT_HALOMOD: Non-linear wavenumber [h/Mpc]:', REAL(hmod%knl)
    END IF

    hmod%neff=effective_index(hmod,cosm)

    IF(verbose) WRITE(*,*) 'INIT_HALOMOD: Collapse n_eff:', REAL(hmod%neff)

    CALL fill_halo_concentration(hmod,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Halo concentration tables filled'
       WRITE(*,*) 'INIT_HALOMOD: Minimum concentration:', REAL(hmod%c(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Maximum concentration:', REAL(hmod%c(1))
    END IF

    IF(slow_hmod) THEN
       A0=one_halo_amplitude(hmod,cosm)
       IF(verbose) THEN
          WRITE(*,*) 'INIT_HALOMOD: One-halo amplitude [Mpc/h]^3:', REAL(A0)
          WRITE(*,*) 'INIT_HALOMOD: One-halo amplitude [log10(M) [Msun/h]]:', REAL(log10(A0*comoving_matter_density(cosm)))
       END IF
    END IF

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE init_halomod

  SUBROUTINE print_halomod(hmod,cosm,verbose)

    !This subroutine writes out the physical halo-model parameters at some redshift 
    !(e.g., Delta_v) rather than the model parameters
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose

    IF(verbose) THEN

       WRITE(*,*) '============================================'
       WRITE(*,*) 'HALOMODEL: ', TRIM(hmod%name)
       WRITE(*,*) '============================================'
       WRITE(*,*) 'HALOMODEL: Accuracy parameters'
       WRITE(*,*) '============================================'
       WRITE(*,*) 'Number of points in look-up tables:', hmod%n
       WRITE(*,*) 'Halo-model accuracy parameter:', hmod%acc_HMx
       WRITE(*,*) 'Large value of nu:', hmod%large_nu
       WRITE(*,*) '============================================'

       ! Form of the two-halo term
       IF(hmod%ip2h==1) WRITE(*,*) 'HALOMODEL: Linear two-halo term'
       IF(hmod%ip2h==2) WRITE(*,*) 'HALOMODEL: Standard two-halo term (Seljak 2000)'
       IF(hmod%ip2h==3) WRITE(*,*) 'HALOMODEL: Linear two-halo term with damped wiggles'
       IF(hmod%ip2h==4) WRITE(*,*) 'HALOMODEL: No two-halo term'

       ! Order to go to in halo bias
       IF(hmod%ip2h==2) THEN
          IF(hmod%ibias==1) WRITE(*,*) 'HALOMODEL: Linear halo bias'
          IF(hmod%ibias==2) WRITE(*,*) 'HALOMODEL: Second-order halo bias'
       END IF

       ! Correction for missing low-mass haloes
       IF(hmod%ip2h==2) THEN
          IF(hmod%ip2h_corr==1) WRITE(*,*) 'HALOMODEL: No two-halo correction applied for missing low-mass haloes'
          IF(hmod%ip2h_corr==2) WRITE(*,*) 'HALOMODEL: Two-halo term corrected by adding missing g(nu)b(nu)' 
          IF(hmod%ip2h_corr==3) WRITE(*,*) 'HALOMODEL: Two-halo term corrected via delta function at low mass end'      
       END IF

       ! 'Scale-dependent halo bias' - not really a halo bias, this is more of a fudge
       IF(hmod%ikb==1) WRITE(*,*) 'HALOMODEL: Scale-independent halo bias'
       IF(hmod%ikb==2) WRITE(*,*) 'HALOMODEL: Scale-dependent halo bias fudge from Fedeli (2014b)'
       IF(hmod%ikb==3) WRITE(*,*) 'HALOMODEL: Scale-dependent halo bias fudge'

       ! Halo mass function
       IF(hmod%imf==1) WRITE(*,*) 'HALOMODEL: Press & Schecter (1974) mass function'
       IF(hmod%imf==2) WRITE(*,*) 'HALOMODEL: Sheth & Tormen (1999) mass function'
       IF(hmod%imf==3) WRITE(*,*) 'HALOMODEL: Tinker et al. (2010) mass function'
       IF(hmod%imf==4) WRITE(*,*) 'HALOMODEL: Delta function mass function'

       ! Concentration-mass relation
       IF(hmod%iconc==1) WRITE(*,*) 'HALOMODEL: Full Bullock et al. (2001) concentration-mass relation'
       IF(hmod%iconc==2) WRITE(*,*) 'HALOMODEL: Simple Bullock et al. (2001) concentration-mass relation'
       IF(hmod%iconc==3) WRITE(*,*) 'HALOMODEL: Full-sample mean density Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==4) WRITE(*,*) 'HALOMODEL: Full-sample virial denity Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==5) WRITE(*,*) 'HALOMODEL: Relaxed-sample x200 critical density Duffy et al. (2008) concentration-mass relation'

       ! Concentration-mass relation correction
       IF(hmod%iDolag==1) WRITE(*,*) 'HALOMODEL: No concentration-mass correction for dark energy'
       IF(hmod%iDolag==2) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction'
       IF(hmod%iDolag==3) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction with 1.5 exponent'

       ! Bound gas fraction
       IF(hmod%frac_boundgas==1) WRITE(*,*) 'HALOMODEL: Halo bound gas fraction: Fedeli (2014a)'
       IF(hmod%frac_boundgas==2) WRITE(*,*) 'HALOMODEL: Halo bound gas fraction: Schneider (2015)'
       IF(hmod%frac_boundgas==3) WRITE(*,*) 'HALOMODEL: Halo bound gas fraction: Universal baryon fraction'

       ! Cold gas fraction
       IF(hmod%frac_coldgas==1) WRITE(*,*) 'HALOMODEL: Cold gas fraction: Constant fraction of halo gas'

       ! Star fraction
       IF(hmod%frac_star==1) WRITE(*,*) 'HALOMODEL: Halo star fraction: Fedeli (2014)'
       IF(hmod%frac_star==2) WRITE(*,*) 'HALOMODEL: Halo star fraction: Constant'
       IF(hmod%frac_star==3) WRITE(*,*) 'HALOMODEL: Halo star fraction: Fedeli (2014) but saturated at high halo mass'
       IF(hmod%frac_star==4) WRITE(*,*) 'HALOMODEL: Halo star fraction: No stars in haloes'

       ! HI fraction
       IF(hmod%frac_HI==1) WRITE(*,*) 'HALOMODEL: Halo HI fraction: Simple model'
       IF(hmod%frac_HI==2) WRITE(*,*) 'HALOMODEL: Halo HI fraction: Villaescusa-Navarro et al. (2018; 1804.09180)'

       ! DMONLY halo model
       IF(hmod%halo_DMONLY==1) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: NFW'
       IF(hmod%halo_DMONLY==2) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Non-analytical NFW (good for testing W(k) functions)'
       IF(hmod%halo_DMONLY==3) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Tophat'
       IF(hmod%halo_DMONLY==4) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Delta function'
       IF(hmod%halo_DMONLY==5) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Cored NFW'
       IF(hmod%halo_DMONLY==6) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Isothermal'
       IF(hmod%halo_DMONLY==7) WRITE(*,*) 'HALOMODEL: DMONLY halo profile: Shell'
       
       ! CDM halo profile
       IF(hmod%halo_CDM==1) WRITE(*,*) 'HALOMODEL: CDM halo profile: NFW'

       ! Bound gas halo profile
       IF(hmod%halo_boundgas==1) WRITE(*,*) 'HALOMODEL: Bound gas profile: Simplified Komatsu & Seljak (2001)'
       IF(hmod%halo_boundgas==2) WRITE(*,*) 'HALOMODEL: Bound gas profile: Isothermal beta profile'
       IF(hmod%halo_boundgas==3) WRITE(*,*) 'HALOMODEL: Bound gas profile: Full Komatsu & Seljak (2001)'

       ! Cold gas profile
       IF(hmod%halo_coldgas==1) WRITE(*,*) 'HALOMODEL: Cold gas profile: Delta function'

       ! Free gas halo profile
       IF(hmod%halo_freegas==1) WRITE(*,*) 'HALOMODEL: Free gas profile: Isothermal model (out to 2rv)'
       IF(hmod%halo_freegas==2) WRITE(*,*) 'HALOMODEL: Free gas profile: Ejected gas model from Schneider (2015)'
       IF(hmod%halo_freegas==3) WRITE(*,*) 'HALOMODEL: Free gas profile: Isothermal shell that connects electron pressure and density to boundgas at rv'
       IF(hmod%halo_freegas==4) WRITE(*,*) 'HALOMODEL: Free gas profile: Komatsu-Seljak continuation'
       IF(hmod%halo_freegas==5) WRITE(*,*) 'HALOMODEL: Free gas profile: Power-law continuation'
       IF(hmod%halo_freegas==6) WRITE(*,*) 'HALOMODEL: Free gas profile: Cubic profile'
       IF(hmod%halo_freegas==7) WRITE(*,*) 'HALOMODEL: Free gas profile: Smoothly distributed'
       IF(hmod%halo_freegas==8) WRITE(*,*) 'HALOMODEL: Free gas profile: Delta function'

       ! Star halo profile
       IF(hmod%halo_star==1) WRITE(*,*) 'HALOMODEL: Star profile: Fedeli (2014) stellar distribution'
       IF(hmod%halo_star==2) WRITE(*,*) 'HALOMODEL: Star profile: Schneider (2015) stellar distribution'
       IF(hmod%halo_star==3) WRITE(*,*) 'HALOMODEL: Star profile: Delta function'
       IF(hmod%halo_star==4) WRITE(*,*) 'HALOMODEL: Star profile: Delta function at low mass and NFW at high mass'

       ! HI halo profile
       IF(hmod%halo_HI==1) WRITE(*,*) 'HALOMODEL: HI profile: NFW profile'
       IF(hmod%halo_HI==2) WRITE(*,*) 'HALOMODEL: HI profile: Delta function'
       IF(hmod%halo_HI==3) WRITE(*,*) 'HALOMODEL: HI profile: Polynomial with exponential hole (Villaescusa-Navarro et al. 2018)'
       IF(hmod%halo_HI==4) WRITE(*,*) 'HALOMODEL: HI profile: Modified NFW with exponential hole (Villaescusa-Navarro et al. 2018)'
       IF(hmod%halo_HI==5) WRITE(*,*) 'HALOMODEL: HI profile: Modified NFW (Padmanabhan & Refregier 2017)'

       ! Electron pressure profile
       IF(hmod%electron_pressure==1) WRITE(*,*) 'HALOMODEL: Electron pressure: Using UPP'
       IF(hmod%electron_pressure==2) WRITE(*,*) 'HALOMODEL: Electron pressure: Using gas profiles'

       ! Are voids being done?
       IF(hmod%voids) THEN
          
          WRITE(*,*) 'HALOMODEL: Considering voids'
          
          ! Void model
          IF(hmod%halo_void==1) WRITE(*,*) 'HALOMODEL: Void model: Top-hat'

          ! Compensated void model
          IF(hmod%halo_compensated_void==1) WRITE(*,*) 'HALOMODEL: Compensated void model: Top-hat'
          
       END IF

       ! delta_c
       IF(hmod%idc==1) WRITE(*,*) 'HALOMODEL: Fixed delta_c = 1.686'
       IF(hmod%idc==2) WRITE(*,*) 'HALOMODEL: delta_c from Nakamura & Suto (1997) fitting function'
       IF(hmod%idc==3) WRITE(*,*) 'HALOMODEL: delta_c from Mead et al. (2015, 2016) power spectrum fit'
       IF(hmod%idc==4) WRITE(*,*) 'HALOMODEL: delta_c from Mead (2017) fitting function'
       IF(hmod%idc==5) WRITE(*,*) 'HALOMODEL: delta_c from spherical-collapse calculation'

       ! Delta_v
       IF(hmod%iDv==1) WRITE(*,*) 'HALOMODEL: Fixed Delta_v = 200'
       IF(hmod%iDv==2) WRITE(*,*) 'HALOMODEL: Delta_v from Bryan & Norman (1998) fitting function'
       IF(hmod%iDv==3) WRITE(*,*) 'HALOMODEL: Delta_v from Mead et al. (2015, 2016) power spectrum fit'
       IF(hmod%iDv==4) WRITE(*,*) 'HALOMODEL: Delta_v from Mead (2017) fitting function'
       IF(hmod%iDv==5) WRITE(*,*) 'HALOMODEL: Delta_v from spherical-collapse calculation'
       IF(hmod%iDv==6) WRITE(*,*) 'HALOMODEL: Delta_v to give haloes Lagrangian radius'

       ! eta for halo window function
       IF(hmod%ieta==1) WRITE(*,*) 'HALOMODEL: eta = 0 fixed'
       IF(hmod%ieta==2) WRITE(*,*) 'HALOMODEL: eta from Mead et al. (2015, 2016) power spectrum fit'

       ! Small-scale two-halo term damping coefficient
       IF(hmod%i2hdamp==1) WRITE(*,*) 'HALOMODEL: No two-halo term damping at small scales'
       IF(hmod%i2hdamp==2) WRITE(*,*) 'HALOMODEL: Two-halo term damping from Mead et al. (2015)'
       IF(hmod%i2hdamp==3) WRITE(*,*) 'HALOMODEL: Two-halo term damping from Mead et al. (2016)'

       ! Large-scale one-halo term damping function
       IF(hmod%i1hdamp==1) WRITE(*,*) 'HALOMODEL: No damping in one-halo term at large scales'
       IF(hmod%i1hdamp==2) WRITE(*,*) 'HALOMODEL: One-halo term large-scale damping via an exponential'
       IF(hmod%i1hdamp==3) WRITE(*,*) 'HALOMODEL: One-halo term large-scale damping like Delta^2 ~ k^7'

       ! Large-scale one-halo term damping coefficient
       IF(hmod%i1hdamp .NE. 1) THEN
          IF(hmod%ikstar==1) WRITE(*,*) 'HALOMODEL: No damping in one-halo term at large scales'
          IF(hmod%ikstar==2) WRITE(*,*) 'HALOMODEL: One-halo term damping function from Mead et al. (2015, 2016)'
       END IF

       ! Concentration-mass scaling
       IF(hmod%iAs==1) WRITE(*,*) 'HALOMODEL: No rescaling of concentration-mass relation'
       IF(hmod%iAs==2) WRITE(*,*) 'HALOMODEL: Concentration-mass relation rescaled mass independetly (Mead et al. 2015, 2016)'

       ! Two- to one-halo transition region
       IF(hmod%itrans==1) WRITE(*,*) 'HALOMODEL: Standard sum of two- and one-halo terms'
       IF(hmod%itrans==2) WRITE(*,*) 'HALOMODEL: Smoothed transition using alpha'
       IF(hmod%itrans==4) WRITE(*,*) 'HALOMODEL: Experimental smoothed transition for HMx'
       IF(hmod%itrans==5) WRITE(*,*) 'HALOMODEL: Tanh transition with k_nl'

       ! Response
       IF(hmod%response) WRITE(*,*) 'HALOMODEL: Treated as response with HMcode multiplication'

       ! Numerical parameters
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: Standard parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'redshift:', hmod%z
       WRITE(*,fmt='(A30,F10.5)') 'scale factor:', hmod%a
       WRITE(*,fmt='(A30,F10.5)') 'Dv:', Delta_v(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'dc:', delta_c(hmod,cosm)      
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMcode parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'eta:', eta(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'k*:', kstar(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'A:', As(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'fdamp:', fdamp(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'alpha:', alpha_transition(hmod,cosm)
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMx parameters'
       WRITE(*,*) '======================================='
       IF(hmod%fixed_HMx) THEN
          WRITE(*,fmt='(A30,F10.5)') 'alpha:', hmod%alpha
          WRITE(*,fmt='(A30,F10.5)') 'alpha index:', hmod%alphap
          WRITE(*,fmt='(A30,F10.5)') 'epsilon:', hmod%eps
          WRITE(*,fmt='(A30,F10.5)') 'Gammma:', hmod%Gamma
          WRITE(*,fmt='(A30,F10.5)') 'Gammma index:', hmod%Gammap
          WRITE(*,fmt='(A30,F10.5)') 'alpha:', hmod%alpha
          WRITE(*,fmt='(A30,F10.5)') 'log10(M0) [Msun/h]:', log10(hmod%M0)
          WRITE(*,fmt='(A30,F10.5)') 'A*:', hmod%Astar
          WRITE(*,fmt='(A30,F10.5)') 'log10(T_WHIM) [K]:', log10(hmod%Twhim)
          WRITE(*,fmt='(A30,F10.5)') 'c*:', hmod%cstar
          WRITE(*,fmt='(A30,F10.5)') 'c* index:', hmod%cstarp
       ELSE
          WRITE(*,fmt='(A30,F10.5)') 'log10(T_heat) [K]:', log10(hmod%Theat)
          WRITE(*,fmt='(A30,F10.5)') 'alpha:', HMx_alpha(1e14,hmod)
          WRITE(*,fmt='(A30,F10.5)') 'epsilon:', HMx_eps(hmod)
          WRITE(*,fmt='(A30,F10.5)') 'Gamma:', HMx_Gamma(1e14,hmod)
          WRITE(*,fmt='(A30,F10.5)') 'log10(M0) [Msun/h]:', log10(HMx_M0(hmod))
          WRITE(*,fmt='(A30,F10.5)') 'A*:', HMx_Astar(hmod)
          WRITE(*,fmt='(A30,F10.5)') 'log10(T_WHIM) [K]:', log10(HMx_Twhim(hmod))
          WRITE(*,fmt='(A30,F10.5)') 'c*:', HMx_cstar(1e14,hmod)
       END IF
       WRITE(*,fmt='(A30,F10.5)') 'sigma*:', hmod%sstar
       WRITE(*,fmt='(A30,F10.5)') 'log10(M*) [Msun/h]:', log10(hmod%Mstar)
       WRITE(*,fmt='(A30,F10.5)') 'f_cold:', hmod%fcold
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HOD parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_gal) [Msun/h]:', log10(hmod%mgal)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_HI_min) [Msun/h]:', log10(hmod%HImin)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_HI_max) [Msun/h]:', log10(hmod%HImax)
       WRITE(*,*) '======================================='
       IF(hmod%halo_DMONLY==5) WRITE(*,fmt='(A30,F10.5)') 'r_core [Mpc/h]:', hmod%rcore
       IF(hmod%dlnc .NE. 0.) WRITE(*,fmt='(A30,F10.5)') 'delta ln(c):', hmod%dlnc
       IF(hmod%imf==4) WRITE(*,*) 'Halo mass log10(M) [Msun/h]:', log10(hmod%hmass)
       WRITE(*,*)
       
    END IF

  END SUBROUTINE print_halomod
  
  SUBROUTINE dewiggle_init(hmod,cosm)

    ! Initialise the dewiggled power spectrum
    IMPLICIT NONE
    !REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: kv(4), pv(4), sigv, a
    REAL, ALLOCATABLE :: Pk(:), Pkraw(:), k(:)
    INTEGER :: i, nk

    IF(.NOT. ALLOCATED(cosm%log_k_plin)) STOP 'DEWIGGLE_INIT: Error, P(k) needs to be tabulated for this to work'

    nk=cosm%n_plin
    hmod%n_pdamp=nk
    sigv=hmod%sigv
    a=hmod%a

    ! Allocate a raw array so as to plot the processed and unprocessed results
    ALLOCATE(Pkraw(nk), Pk(nk), k(nk))

    ! Allocate the internal arrays from the cosmology arrays
    k=exp(cosm%log_k_plin)
    Pk=exp(cosm%log_plin)

    ! Fixed k values - CAMB!
    kv(1)=0.008
    kv(2)=0.01
    kv(3)=0.8
    kv(4)=1.0

    ! Fix p values
    DO i=1,4
       pv(i)=exp(find(log(kv(i)),log(k),log(Pk),nk,3,3,2))      
    END DO

    ! Create new 'raw' spectrum which has no wiggles
    DO i=1,nk
       IF(k(i)<=kv(1) .OR. k(i)>=kv(4)) THEN
          Pkraw(i)=Pk(i)
       ELSE
          Pkraw(i)=exp(Lagrange_polynomial(log(k(i)),3,log(kv),log(pv)))
       END IF
    END DO

    ! Isolate just the wiggles
    Pk=Pk-Pkraw

    ! Damp the wiggles
    Pk=Pk*exp(-(sigv*k)**2)

    ! Add the damped wiggles back in
    Pk=Pk+Pkraw

    ! Create the damped power array
    IF(ALLOCATED(hmod%log_k_pdamp)) DEALLOCATE(hmod%log_k_pdamp)
    IF(ALLOCATED(hmod%log_pdamp))   DEALLOCATE(hmod%log_pdamp)
    ALLOCATE(hmod%log_k_pdamp(nk),hmod%log_pdamp(nk))

    ! Fill the k array
    hmod%log_k_pdamp=log(k)
    
    ! Grow damped power to the correct redshift and fill array
    hmod%log_pdamp=log(Pk*grow(a,cosm)**2)

    ! Set the flag
    hmod%has_dewiggle=.TRUE.

  END SUBROUTINE dewiggle_init

  REAL FUNCTION p_dewiggle(k,hmod,cosm)

    ! Call the dewiggled power spectrum
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%has_dewiggle .EQV. .FALSE.) CALL dewiggle_init(hmod,cosm)
    p_dewiggle=exp(find(log(k),hmod%log_k_pdamp,hmod%log_pdamp,hmod%n_pdamp,3,3,2))

  END FUNCTION p_dewiggle

  FUNCTION halo_type(i)

    ! Name function for halo types
    ! TODO: This must be able to be combined with set_halo_type
    ! TODO: Can this be in the header?
    IMPLICIT NONE
    CHARACTER(len=256) :: halo_type
    INTEGER :: i
    
    halo_type=''
    IF(i==field_dmonly)             halo_type='DMONLY'
    IF(i==field_matter)             halo_type='Matter'
    IF(i==field_cdm)                halo_type='CDM'
    IF(i==field_gas)                halo_type='Gas'
    IF(i==field_star)               halo_type='Star'
    IF(i==field_bound_gas)          halo_type='Bound gas'
    IF(i==field_free_gas)           halo_type='Free gas'
    IF(i==field_electron_pressure)  halo_type='Electron pressure'
    IF(i==field_void)               halo_type='Void'
    IF(i==field_compensated_void)   halo_type='Compensated void'
    IF(i==field_central_galaxies)   halo_type='Central galaxies'
    IF(i==field_satellite_galaxies) halo_type='Satellite galaxies'
    IF(i==field_galaxies)           halo_type='Galaxies'
    IF(i==field_HI)                 halo_type='HI'
    IF(i==field_cold_gas)           halo_type='Cold gas'
    IF(halo_type=='') STOP 'HALO_TYPE: Error, i not specified correctly'
    
  END FUNCTION halo_type

  SUBROUTINE set_halo_type(ip)

    ! Set the halo types
    ! TODO: This must be able to be combined with halo_type
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ip
    INTEGER :: i

    INTEGER, PARAMETER :: halo_min_i=-1
    INTEGER, PARAMETER :: halo_max_i=12

    WRITE(*,*) 'SET_HALO_TYPE: Choose halo type'
    WRITE(*,*) '==============================='
    DO i=halo_min_i,halo_max_i
       WRITE(*,fmt='(I3,A3,A26)') i, '- ', TRIM(halo_type(i))
    END DO
    READ(*,*) ip
    WRITE(*,*) '==============================='
    WRITE(*,*)

    IF(ip<halo_min_i .OR. ip>halo_max_i) STOP 'SET_HALO_TYPE: Error, you have chosen a bad halo'

  END SUBROUTINE set_halo_type

  SUBROUTINE calculate_HMx(itype,nt,mmin,mmax,k,nk,a,na,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response)

    ! Public facing function, calculates the halo model power for the desired 'k' and 'a' range
    ! TODO: Change :,:,i to i,:,: for speed
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype(nt)
    INTEGER, INTENT(IN) :: nt
    REAL, INTENT(IN) :: mmin, mmax
    REAL, INTENT(IN) :: k(nk)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: a(na)
    INTEGER, INTENT(IN) :: na
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_li(:,:)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_2h(:,:,:,:)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_1h(:,:,:,:)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow_hm(:,:,:,:)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    LOGICAL, INTENT(IN) :: response
    REAL :: z
    INTEGER :: i, j
    LOGICAL :: verbose2

    ! To avoid splurge of stuff printed to screen
    verbose2=verbose

    ! Deallocate arrays if they are already allocated
    IF(ALLOCATED(pow_li)) DEALLOCATE(pow_li)
    IF(ALLOCATED(pow_2h)) DEALLOCATE(pow_2h)
    IF(ALLOCATED(pow_1h)) DEALLOCATE(pow_1h)
    IF(ALLOCATED(pow_hm)) DEALLOCATE(pow_hm)

    ! Allocate power arrays
    ALLOCATE(pow_li(nk,na),pow_2h(nt,nt,nk,na),pow_1h(nt,nt,nk,na),pow_hm(nt,nt,nk,na))

    ! Do the halo-model calculation
    DO i=na,1,-1
       z=redshift_a(a(i))
       CALL init_halomod(mmin,mmax,a(i),hmod,cosm,verbose2)
       CALL print_halomod(hmod,cosm,verbose2)
       CALL calculate_HMx_a(itype,nt,k,nk,pow_li(:,i),pow_2h(:,:,:,i),pow_1h(:,:,:,i),pow_hm(:,:,:,i),hmod,cosm,verbose2,response) ! Slow
       IF(i==na .and. verbose) THEN
          WRITE(*,*) 'CALCULATE_HMx: Doing calculation'
          DO j=1,nt
             WRITE(*,*) 'CALCULATE_HMx: Haloes:', itype(j), TRIM(halo_type(itype(j)))
          END DO
          WRITE(*,*) '======================================='
          WRITE(*,*) '                            a         z'
          WRITE(*,*) '======================================='
       END IF
       IF(verbose) WRITE(*,fmt='(A15,I5,2F10.3)') 'CALCULATE_HMx:', i, a(i), z
       verbose2=.FALSE.
    END DO
    IF(verbose) THEN
       WRITE(*,*) '======================================='
       WRITE(*,*) 'CALCULATE_HMx: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_HMx

  SUBROUTINE calculate_HMcode_a(k,a,Pk,nk,cosm)

    ! Get the HMcode prediction at this z for this cosmology
    IMPLICIT NONE
    REAL, INTENT(IN) :: k(nk)
    REAL, INTENT(IN) :: a
    REAL, INTENT(OUT) :: Pk(nk)
    INTEGER, INTENT(IN) :: nk
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: pow_lin(nk), pow_2h(nk), pow_1h(nk)
    TYPE(halomod) :: hmod
    
    INTEGER :: ihm=1 ! Set HMcode, could/should be parameter
    INTEGER, PARAMETER :: dmonly(1)=field_dmonly ! DMONLY

    ! Do an HMcode run
    CALL assign_halomod(ihm,hmod,verbose=.FALSE.)
    CALL init_halomod(mmin_hmcode,mmax_hmcode,a,hmod,cosm,verbose=.FALSE.)
    CALL calculate_HMx_a(dmonly,1,k,nk,pow_lin,pow_2h,pow_1h,Pk,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

  END SUBROUTINE calculate_HMcode_a

  SUBROUTINE calculate_HMx_a(itype,nt,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype(nt)
    INTEGER, INTENT(IN) :: nt
    REAL, INTENT(IN) :: k(nk)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(OUT) :: pow_li(nk)
    REAL, INTENT(OUT) :: pow_2h(nt,nt,nk)
    REAL, INTENT(OUT) :: pow_1h(nt,nt,nk)
    REAL, INTENT(OUT) :: pow_hm(nt,nt,nk)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    LOGICAL, INTENT(IN) :: response
    REAL :: plin
    REAL :: powg_2h(nk), powg_1h(nk), powg_hm(nk)
    REAL :: hmcode_2h(nk), hmcode_1h(nk), hmcode_hm(nk)
    INTEGER :: i, ihmcode
    TYPE(halomod) :: hmcode
    
    INTEGER, PARAMETER :: dmonly(1)=field_dmonly ! DMONLY
    
    ! Write to screen
    IF(verbose) THEN
       DO i=1,nt
          WRITE(*,*) 'CALCULATE_HMX_A: Halo type:', itype(i), TRIM(halo_type(itype(i)))
       END DO
       WRITE(*,*) 'CALCULATE_HMX_A: k min [h/Mpc]:', REAL(k(1))
       WRITE(*,*) 'CALCULATE_HMX_A: k max [h/Mpc]:', REAL(k(nk))
       WRITE(*,*) 'CALCULATE_HMX_A: number of k:', nk
       WRITE(*,*) 'CALCULATE_HMX_A: a:', REAL(hmod%a)
       WRITE(*,*) 'CALCULATE_HMX_A: z:', REAL(hmod%z)
       WRITE(*,*) 'CALCULATE_HMX_A: Calculating halo-model power spectrum'
       WRITE(*,*)
    END IF

    ! Do an HMcode calculation for multiplying the response
    ! This should probably be moved to init_halomod (recursive...)
    ! Annoying having to fix mmin and mmax separately here
    IF(hmod%response) THEN
       ihmcode=1
       CALL assign_halomod(ihmcode,hmcode,verbose=.FALSE.)
       CALL init_halomod(mmin_response,mmax_response,hmod%a,hmcode,cosm,verbose=.FALSE.)
    END IF
    
    ! Loop over k values
    !TODO: add OMP support properly. What is private and what is shared? CHECK THIS!
!!$OMP PARALLEL DO DEFAULT(SHARED)!, private(k,plin,pow_2h,pow_1h,pow,pow_lin)
!!$OMP PARALLEL DO DEFAULT(PRIVATE)
!!$OMP PARALLEL DO FIRSTPRIVATE(nk,cosm,compute_p_lin,k,a,pow_lin,plin,itype1,itype2,z,pow_2h,pow_1h,pow,hmod)
!!$OMP PARALLEL DO
    DO i=1,nk

       ! Get the linear power
       plin=p_lin(k(i),hmod%a,cosm)
       pow_li(i)=plin

       ! Do the halo model calculation
       CALL calculate_HMx_ka(itype,nt,k(i),plin,pow_2h(:,:,i),pow_1h(:,:,i),pow_hm(:,:,i),hmod,cosm) ! (slow array accessing)
          
       IF(response .OR. hmod%response) THEN

          ! If doing a response then calculate a DMONLY prediction too
          CALL calculate_HMx_ka(dmonly,1,k(i),plin,powg_2h(i),powg_1h(i),powg_hm(i),hmod,cosm)
          pow_li(i)=1.                               ! This is just linear-over-linear, which is one
          pow_2h(:,:,i)=pow_2h(:,:,i)/powg_2h(i) ! Two-halo response (slow array accessing)
          pow_1h(:,:,i)=pow_1h(:,:,i)/powg_1h(i) ! One-halo response (slow array accessing)
          pow_hm(:,:,i)=pow_hm(:,:,i)/powg_hm(i) ! Full model response (slow array accessing)

          IF((.NOT. response) .AND. hmod%response) THEN

             ! If multiplying the response by an 'accurate' HMcode prediction
             CALL calculate_HMx_ka(dmonly,1,k(i),plin,hmcode_2h(i),hmcode_1h(i),hmcode_hm(i),hmcode,cosm)
             pow_li(i)=plin                               ! Linear power is just linear power again
             pow_2h(:,:,i)=pow_2h(:,:,i)*hmcode_2h(i) ! Multiply two-halo response through by HMcode two-halo term (slow array accessing)
             pow_1h(:,:,i)=pow_1h(:,:,i)*hmcode_1h(i) ! Multiply one-halo response through by HMcode one-halo term (slow array accessing)
             pow_hm(:,:,i)=pow_hm(:,:,i)*hmcode_hm(i) ! Multiply response through by HMcode (slow array accessing)
             
          END IF
          
       END IF

    END DO
!!$OMP END PARALLEL DO

  END SUBROUTINE calculate_HMx_a

  SUBROUTINE calculate_HMx_ka(itype,nt,k,plin,pow_2h,pow_1h,pow_hm,hmod,cosm)

    ! Gets the one- and two-halo terms and combines them
    ! TODO: Re-support scatter
    ! TODO: include scatter in two-halo term
    ! TODO: Can I avoid calling Window twice for smooth gas etc. ?
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype(nt)
    INTEGER, INTENT(IN) :: nt
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: plin
    REAL, INTENT(OUT) :: pow_2h(nt,nt)
    REAL, INTENT(OUT) :: pow_1h(nt,nt)
    REAL, INTENT(OUT) :: pow_hm(nt,nt)    
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: wk(hmod%n,nt), wk2(hmod%n,2), wk_squared(hmod%n)
    INTEGER :: i, j, j1, j2, ih(2)

    ! Calls expressions for one- and two-halo terms and then combines to form the full power spectrum
    IF(k==0.) THEN

       ! This should really never be called for k=0
       pow_2h=0.
       pow_1h=0.

    ELSE
       
          ! Get the window functions
          CALL init_windows(k,itype,nt,wk,hmod%n,hmod,cosm)

          ! Loop over fields and get the one-halo term for each pair
          DO j1=1,nt
             DO j2=j1,nt
                IF(hmod%dlnc==0.) THEN
                   wk_squared=wk(:,j1)*wk(:,j2)
                ELSE
                   ih(1)=itype(j1)
                   ih(2)=itype(j2)
                   wk_squared=wk_squared_scatter(hmod%n,ih,k,hmod,cosm)
                END IF
                pow_1h(j1,j2)=p_1h(wk_squared,k,hmod,cosm)            
             END DO
          END DO

       ! If linear theory is used for two-halo term we need to recalculate the window functions for the two-halo term with k=0 fixed
       ! TODO: Include what to do with ip2h=1,3 properly
       IF(hmod%ip2h==1 .OR. hmod%ip2h==3 .OR. hmod%halo_freegas==7) THEN

          CALL reinit_windows(itype,nt,wk,hmod%n,hmod,cosm)
          
       END IF

       ! Get the two-halo term
       DO i=1,nt
          DO j=i,nt
             ih(1)=itype(i)
             ih(2)=itype(j)
             wk2(:,1)=wk(:,i)
             wk2(:,2)=wk(:,j)
             pow_2h(i,j)=p_2h(ih,wk2,hmod%n,k,plin,hmod,cosm)
          END DO
       END DO

    END IF

    ! Loop over fields and get the total halo-model power
    DO i=1,nt
       DO j=i,nt
          pow_hm(i,j)=p_hm(k,pow_2h(i,j),pow_1h(i,j),hmod,cosm)
       END DO
    END DO

    ! Construct symmetric parts using ij=ji symmetry of spectra
    DO j=1,nt
       DO i=j+1,nt
          pow_1h(i,j)=pow_1h(j,i)
          pow_2h(i,j)=pow_2h(j,i)
          pow_hm(i,j)=pow_hm(j,i)
       END DO
    END DO

  END SUBROUTINE calculate_HMx_ka

  SUBROUTINE init_windows(k,fields,nf,wk,nm,hmod,cosm)

    ! Fill the window functions for all the different fields
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    INTEGER, INTENT(IN) :: fields(nf)
    REAL, INTENT(OUT) :: wk(nm,nf)
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nm
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    REAL :: m, rv, c, rs, nu, et
    INTEGER :: i_all, i_cdm, i_gas, i_sta
    LOGICAL :: quick_matter

    ! This should be set to false initially
    quick_matter=.FALSE.

    ! Get the array positions corresponding to all, cdm, gas, stars if they exist
    i_all=array_position(field_matter,fields,nf)
    i_cdm=array_position(field_cdm,fields,nf)
    i_gas=array_position(field_gas,fields,nf)
    i_sta=array_position(field_star,fields,nf)

    ! If all, cdm, gas and stars exist then activate the quick-matter mode
    IF((i_all .NE. 0) .AND. (i_cdm .NE. 0) .AND. (i_gas .NE. 0) .AND. (i_sta .NE. 0)) THEN
       quick_matter=.TRUE.
    END IF

    ! Get eta
    et=eta(hmod,cosm)

    ! Calculate the halo window functions for each field
    DO j=1,nf

       IF(quick_matter .AND. j==i_all) CYCLE

       ! Loop over masses to fill window-function array
       DO i=1,nm
          
          m=hmod%m(i)
          rv=hmod%rv(i)
          c=hmod%c(i)
          rs=rv/c
          nu=hmod%nu(i)
          
          wk(i,j)=win_type(.FALSE.,fields(j),k*nu**et,m,rv,rs,hmod,cosm)
          
       END DO

    END DO

    ! If quick-matter mode is active then create the total matter window by summing contributions
    IF(quick_matter) THEN
       wk(:,i_all)=wk(:,i_cdm)+wk(:,i_gas)+wk(:,i_sta)
    END IF

  END SUBROUTINE init_windows

  SUBROUTINE reinit_windows(fields,nf,wk,nm,hmod,cosm)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: fields(nf)
    REAL, INTENT(INOUT) :: wk(nm,nf)
    INTEGER, INTENT(IN) :: nf
    INTEGER, INTENT(IN) :: nm
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m, rv, c, rs, nu, et, fc, pc
    REAL :: rho0, T0
    INTEGER :: i
    INTEGER :: i_all, i_gas, i_pre
    
    ! Get the array positions corresponding to all, cdm, gas, stars if they exist
    i_all=array_position(0,fields,nf)
    i_gas=array_position(2,fields,nf)
    i_pre=array_position(6,fields,nf)

    IF(i_all==0 .AND. i_gas==0 .AND. i_pre==0) THEN

       ! Do nothing

    ELSE

       ! Get eta
       et=eta(hmod,cosm)
       
       ! Loop over mass and apply corrections

       DO i=1,hmod%n

          ! Halo variables
          m=hmod%m(i)
          rv=hmod%rv(i)
          c=hmod%c(i)
          rs=rv/c
          nu=hmod%nu(i)
          
          IF((i_all .NE. 0) .OR. (i_gas .NE. 0)) THEN

             ! Correction factor for the gas density profiles
             fc=halo_freegas_fraction(m,hmod,cosm)*m/comoving_matter_density(cosm)

             ! Apply the correction factor
             IF(i_all .NE. 0) wk(i,i_all)=wk(i,i_all)+fc
             IF(i_gas .NE. 0) wk(i,i_gas)=wk(i,i_gas)+fc
             
          END IF

          IF(i_pre .NE. 0) THEN

             ! TODO: The units here are a mess
          
             ! Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
             rho0=m*halo_freegas_fraction(m,hmod,cosm)!/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density) ! rho0 in [(Msun/h)/(Mpc/h)^3]
             rho0=rho0*msun/Mpc/Mpc/Mpc ! Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
             rho0=rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

             ! This is the total thermal pressure of the WHIM
             T0=HMx_Twhim(hmod) ! [K]

             ! Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
             pc=(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
             pc=pc/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
             pc=pc*cosm%mue/cosm%mup ! Convert from total thermal pressure to electron pressure

             ! Apply the correction factor
             wk(i,i_pre)=wk(i,i_pre)+pc

          END IF
          
       END DO

    END IF
    
  END SUBROUTINE reinit_windows

  FUNCTION wk_squared_scatter(n,itype,k,hmod,cosm)

    IMPLICIT NONE
    REAL :: wk_squared_scatter(n)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: itype(2)
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: m, rv, c

    DO i=1,n
       m=hmod%m(i)
       rv=hmod%rv(i)
       c=hmod%c(i)
       wk_squared_scatter(i)=integrate_scatter(c,hmod%dlnc,itype,k,m,rv,hmod,cosm,hmod%acc_HMx,3)
    END DO

  END FUNCTION wk_squared_scatter

  REAL FUNCTION p_2h(ih,wk,n,k,plin,hmod,cosm)

    ! Produces the 'two-halo' power
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: wk(n,2)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: plin
    TYPE(halomod), INTENT(INOUT) :: hmod    
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: sigv, frac, rhom
    REAL :: nu, m, m0, wki(2), b0, wk0(2), nu0
    REAL :: integrand1(n,2), integrand2(n,2)
    REAL :: sum1(2), sum2(2)
    INTEGER :: i, j

    rhom=comoving_matter_density(cosm)

    IF(hmod%ip2h==1) THEN

       ! Simply linear theory
       ! TODO: This should be multiplied by a W(k->0) thing to work for non-matter fields
       p_2h=plin

    ELSE IF(hmod%ip2h==3) THEN

       ! Damped BAO linear theory
       ! TODO: This should be multiplied by a W(k->0) thing to work for non-matter fields
       p_2h=p_dewiggle(k,hmod,cosm)

    ELSE IF(hmod%ip2h==4) THEN

       ! No two-halo term
       p_2h=0.
       
    ELSE

       IF(hmod%imf==4) THEN

          ! In this case the mass function is a delta function...
          
          m0=hmod%hmass
          nu0=nu_M(m0,hmod,cosm)
          b0=b_nu(nu0,hmod)
          wk0(1)=find(log(m0),hmod%log_m,wk(:,1),n,3,3,2)
          wk0(2)=find(log(m0),hmod%log_m,wk(:,2),n,3,3,2)
          sum1(1)=rhom*b0*wk0(1)/m0
          sum1(2)=rhom*b0*wk0(2)/m0

       ELSE

          ! ...otherwise you need to do an integral
          
          DO i=1,n

             ! Some variables to make equations cleaner below
             m=hmod%m(i)
             nu=hmod%nu(i)
             wki=wk(i,:) ! (slow array accessing)

             DO j=1,2

                ! Linear bias term, standard two-halo term integral
                integrand1(i,j)=rhom*g_nu(nu,hmod)*b_nu(nu,hmod)*wki(j)/m

                IF(hmod%ibias==2) THEN
                   ! Second-order bias term
                   integrand2(i,j)=rhom*g_nu(nu,hmod)*b2_nu(nu,hmod)*wki(j)/m
                END IF

             END DO

          END DO

          ! Evaluate these integrals from the tabulated values
          DO j=1,2
             sum1(j)=integrate_table(hmod%nu,integrand1(:,j),n,1,n,3)
          END DO

          IF(hmod%ip2h_corr==1) THEN
             ! Do nothing in this case
             ! There will be large errors if any signal is from low-mass haloes
             ! e.g., for the matter power spectrum
          ELSE IF(hmod%ip2h_corr==2) THEN
             ! Add on the value of integral b(nu)*g(nu) assuming W(k)=1
             ! Advised by Yoo et al. (????) and Cacciato et al. (2012)
             STOP 'P_2H: This will not work for fields that do not have mass fractions defined'
             DO j=1,2
                sum1(j)=sum1(j)+hmod%gbmin*halo_fraction(ih(j),m,hmod,cosm)
             END DO
          ELSE IF(hmod%ip2h_corr==3) THEN
             ! Put the missing part of the integrand as a delta function at the low-mass limit of the integral
             ! I think this is the best thing to do
             m0=hmod%m(1)
             wki=wk(1,:) ! (slow array accessing)
             DO j=1,2             
                sum1(j)=sum1(j)+rhom*hmod%gbmin*wki(j)/m0
             END DO
          ELSE
             STOP 'P_2h: Error, ip2h_corr not specified correctly'
          END IF

       END IF

       p_2h=plin*sum1(1)*sum1(2)

       ! Second-order bias correction
       ! This needs to have the property that \int f(nu)b2(nu) du = 0
       ! This means it is hard to check that the normalisation is correct
       ! e.g., how much do low mass haloes matter
       IF(hmod%ibias==2) THEN
                  
          ! Varying mmin *does* make a difference to the values of the integrals
          DO j=1,2
             sum2(j)=integrate_table(hmod%nu,integrand2(:,j),n,1,n,3)
          END DO
          p_2h=p_2h+(plin**2)*sum2(1)*sum2(2)*rhom**2
          
       END IF

    END IF

    ! Apply the damping to the two-halo term
    IF(hmod%i2hdamp .NE. 1) THEN
       ! Two-halo damping parameters
       sigv=hmod%sigv
       frac=fdamp(hmod,cosm)
       IF(frac==0.) THEN
          p_2h=p_2h
       ELSE
          p_2h=p_2h*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2)
       END IF
    END IF
       
    IF(hmod%ikb==2) THEN
       ! Fedeli (2014b) 'scale-dependent halo bias'
       p_2h=p_2h*(1.+k)**0.54
    ELSE IF(hmod%ikb==3) THEN
       ! My own experimental model
       p_2h=p_2h*(1.+(k/hmod%knl))**0.5
    END IF

    ! For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0
    IF(p_2h<0.) THEN
       WRITE(*,*) 'P_2H: Halo type 1:', ih(1)
       WRITE(*,*) 'P_2H: Halo type 1:', ih(2)
       WRITE(*,*) 'P_2H: k [h/Mpc]:', k
       WRITE(*,*) 'P_2H: z:', hmod%z
       WRITE(*,*) 'P_2H: Delta^2_{2H}:', p_2h
       WRITE(*,*) 'P_2H: Caution! P_2h < 0, this was previously fixed by setting P_2h = 0 explicitly'
       !p_2h=0.
       STOP
    END IF

  END FUNCTION p_2h

  REAL FUNCTION p_1h(wk2,k,hmod,cosm)

    ! Calculates the one-halo term
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL, INTENT(IN) :: wk2(hmod%n) ! This needs to be defined after halomod type, annoyingly
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m, g, fac, ks, wk20, m0, rhom
    REAL, ALLOCATABLE :: integrand(:)
    INTEGER :: i

    INTEGER, PARAMETER :: iorder=1 ! Use basic trapezium rule because the integrand is messy due to possible rapid oscillations in W(k)

    ! Matter density
    rhom=comoving_matter_density(cosm)

    IF(hmod%imf==4) THEN

       ! In this case the mass function is a delta function...
       
       m0=hmod%hmass
       wk20=find(log(m0),hmod%log_m,wk2,hmod%n,3,3,2)
       p_1h=rhom*wk20/m0

    ELSE

       ! ...otherwise you need to do an integral

       ALLOCATE(integrand(hmod%n))
       integrand=0.

       ! Calculates the value of the integrand at all nu values!
       DO i=1,hmod%n
          g=g_nu(hmod%nu(i),hmod)
          m=hmod%m(i)
          integrand(i)=g*wk2(i)/m
       END DO

       ! Carries out the integration   
       p_1h=rhom*integrate_table(hmod%nu,integrand,hmod%n,iorder,hmod%n,1)

       DEALLOCATE(integrand)

    END IF

    ! Convert from P(k) -> Delta^2(k)
    p_1h=p_1h*(4.*pi)*(k/twopi)**3

    ! Damping of the 1-halo term at very large scales
    ks=kstar(hmod,cosm)       

    IF(ks>0.) THEN

       IF(hmod%i1hdamp==1) THEN
          ! Do nothing in this case
       ELSE IF(hmod%i1hdamp==2) THEN
          IF((k/ks)**2>7.) THEN
             ! Prevents problems if k/ks is very large
             fac=0.
          ELSE
             fac=exp(-((k/ks)**2))
          END IF
          p_1h=p_1h*(1.-fac)
       ELSE IF(hmod%i1hdamp==3) THEN
          ! Note that the power here should be 4 because it multiplies Delta^2(k) ~ k^3 at low k (NOT 7)
          ! Want f(k<<ks) ~ k^4; f(k>>ks) = 1
          fac=1./(1.+(ks/k)**4)
          p_1h=p_1h*fac
       ELSE
          STOP 'P_1H: Error, i1hdamp not specified correctly'          
       END IF

    END IF

  END FUNCTION p_1h

  REAL FUNCTION p_hm(k,pow_2h,pow_1h,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: pow_2h
    REAL, INTENT(IN) :: pow_1h
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: alpha

    ! alpha is set to one sometimes, which is just the standard halo-model sum of terms
    ! No need to have an IF statement around this
    IF(hmod%itrans==2 .OR. hmod%itrans==4) THEN

       ! If either term is less than zero then we need to be careful
       IF(pow_2h<0. .OR. pow_1h<0.) THEN

          ! Either the safe option and just sum the components
          IF(hmod%safe_negative) THEN
             p_hm=pow_2h+pow_1h
          ELSE
             WRITE(*,*) 'P_HM: Two-halo term:', pow_2h
             WRITE(*,*) 'P_1h: One-halo term:', pow_1h
             STOP 'CALCULATE_HALOMOD_K: Error, either pow_2h or pow_1h is less than zero, which is a problem for the smoothed transition'
          END IF
          
       ELSE

          ! Do the standard smoothed transition
          alpha=alpha_transition(hmod,cosm)
          p_hm=(pow_2h**alpha+pow_1h**alpha)**(1./alpha)
          
       END IF
       
    ELSE IF(hmod%itrans==5) THEN

       ! Sigmoid transition
       p_hm=pow_2h+sigmoid_log(1.*k/hmod%knl,1.)*(pow_1h-pow_2h)
       
    ELSE

       ! Standard sum
       p_hm=pow_2h+pow_1h
       
    END IF

    ! If we are adding in compensated voids
    IF(hmod%voids) THEN
       p_hm=p_hm+p_1void(k,hmod)
    END IF

  END FUNCTION p_hm

  REAL FUNCTION p_1void(k,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc, wk, V, rvoid, rcomp, nu
    REAL :: integrand(hmod%n)
    INTEGER :: i, n

    REAL, PARAMETER :: dv=-3.
    REAL, PARAMETER :: fvoid=1.1
    LOGICAL, PARAMETER :: compensate=.TRUE.
    LOGICAL, PARAMETER :: simple=.FALSE.

    IF(simple) THEN
       n=1
    ELSE
       n=hmod%n
    END IF

    DO i=1,n

       !Get the void radius and compensation radius
       IF(simple) THEn
          rvoid=10.
       ELSE         
          rvoid=hmod%rr(i)
          nu=hmod%nu(i)        
       END IF
       rcomp=fvoid*rvoid

       !Calculate the compensation over-density
       dc=-dv*rvoid**3/(rcomp**3-rvoid**3)

       !Calculate the void Fourier transform
       IF(compensate) THEN
          wk=(4.*pi/3.)*((dv-dc)*wk_tophat(k*rvoid)*rvoid**3+dc*wk_tophat(k*rcomp)*rcomp**3)
       ELSE
          wk=(4.*pi/3.)*dv*wk_tophat(k*rvoid)*rvoid**3
       END IF

       !Calculate the void volume
       IF(compensate) THEN
          V=rcomp**3
       ELSE
          V=rvoid**3
       END IF

       IF(simple .EQV. .FALSE.) THEN
          integrand(i)=g_nu(nu,hmod)*wk**2/V
       END IF

    END DO

    !Calculate the void one-halo term
    IF(simple) THEN
       p_1void=wk**2/V
    ELSE
       p_1void=integrate_table(hmod%nu,integrand,n,1,n,1)
    END IF

    p_1void=p_1void*(4.*pi)*(k/twopi)**3

  END FUNCTION p_1void
    
  SUBROUTINE halo_diagnostics(hmod,cosm,dir)

    ! Writes out to file a whole set of halo diagnostics
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    CHARACTER(len=*), INTENT(IN) :: dir
    REAL :: mass    
    CHARACTER(len=64) :: ext
    CHARACTER(len=512) :: base
    CHARACTER(len=1024) :: outfile
    INTEGER :: m, mi1, mi2

    ! Integer 10^m to produce haloes between
    REAL, PARAMETER :: m1=1e10
    REAL, PARAMETER :: m2=1e16

    WRITE(*,*) 'HALO_DIAGNOSTICS: Outputting diagnostics'

    outfile=TRIM(dir)//'/mass_fractions.dat'
    WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
    CALL write_mass_fractions(hmod,cosm,outfile)

    IF(hmod%z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(hmod%z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(hmod%z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(hmod%z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_DIAGNOSTICS: Error, need to make this better with z'
    END IF

    mi1=NINT(log10(m1))
    mi2=NINT(log10(m2))
    
    DO m=mi1,mi2

       mass=10.**m

       base=TRIM(dir)//'/halo_profile_m'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_profiles(mass,hmod,cosm,outfile)

       base=TRIM(dir)//'/halo_window_m'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_transforms(mass,hmod,cosm,outfile)

    END DO

    WRITE(*,*) 'HALO_DIAGNOSTICS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_diagnostics

  SUBROUTINE halo_definitions(hmod,cosm,dir)

    ! Writes out to files the different halo definitions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: dir
    CHARACTER(len=256) :: fradius, fmass, fconc
    CHARACTER(len=64) :: ext
    INTEGER :: i

    WRITE(*,*) 'HALO_DEFINITIONS: Outputting definitions'

    IF(hmod%z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(hmod%z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(hmod%z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(hmod%z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_DEFINITIONS: Error, need to make this better with z'
    END IF

    fradius=TRIM(dir)//'/radius'//TRIM(ext)
    fmass=TRIM(dir)//'/mass'//TRIM(ext)
    fconc=TRIM(dir)//'/concentration'//TRIM(ext)

    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fradius)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fmass)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fconc)
    
    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod,cosm)

    OPEN(7,file=fradius)
    OPEN(8,file=fmass)
    OPEN(9,file=fconc)
    DO i=1,hmod%n
       WRITE(7,*) hmod%rv(i), hmod%r200(i), hmod%r500(i), hmod%r200c(i), hmod%r500c(i)
       WRITE(8,*) hmod%m(i),  hmod%m200(i), hmod%m500(i), hmod%m200c(i), hmod%m500c(i)
       WRITE(9,*) hmod%c(i),  hmod%c200(i), hmod%c500(i), hmod%c200c(i), hmod%c500c(i)
    END DO
    CLOSE(7)
    CLOSE(8)
    CLOSE(9)

    WRITE(*,*) 'HALO_DEFINITIONS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_definitions

  SUBROUTINE halo_properties(hmod,dir)

    !Writes out to files the different halo definitions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    CHARACTER(len=*), INTENT(IN) :: dir
    CHARACTER(len=256) :: output
    CHARACTER(len=64) :: ext
    INTEGER :: i

    WRITE(*,*) 'HALO_PROPERTIES: Outputting definitions'

    IF(hmod%z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(hmod%z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(hmod%z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(hmod%z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_PROPERTIES: Error, need to make this better with z'
    END IF

    output=TRIM(dir)//'/properies'//TRIM(ext)
    WRITE(*,*) 'HALO_PROPERTIES: ', TRIM(output)

    OPEN(7,file=output)
    DO i=1,hmod%n
       WRITE(7,*) hmod%m(i), hmod%rr(i), hmod%rv(i), hmod%nu(i), hmod%c(i),hmod%sig(i), hmod%sigf(i), hmod%zc(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'HALO_PROPERTIES: Done'
    WRITE(*,*)

  END SUBROUTINE halo_properties

  SUBROUTINE write_mass_fractions(hmod,cosm,outfile)

    ! Writes out the halo mass fractions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile   
    REAL :: m
    INTEGER :: i, j

    REAL, PARAMETER :: mmin=1e10
    REAL, PARAMETER :: mmax=1e16
    INTEGER, PARAMETER :: n=101

    OPEN(7,file=outfile)
    DO i=1,n
       m=exp(progression(log(mmin),log(mmax),i,n))
       WRITE(7,*) m, (halo_fraction(j,m,hmod,cosm), j=1,5)
    END DO
    CLOSE(7)

  END SUBROUTINE write_mass_fractions

  SUBROUTINE write_halo_profiles(m,hmod,cosm,outfile)

    ! Writes out the halo density profiles
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile    
    REAL :: r, rv, rs, c
    INTEGER :: i, j, nf
    INTEGER, ALLOCATABLE :: fields(:)

    REAL, PARAMETER :: rmin=1e-3     ! Mininum r/rv
    REAL, PARAMETER :: rmax=1.1e0    ! Maximum r/rv
    INTEGER, PARAMETER :: n=512      ! Number of points
    LOGICAL, PARAMETER :: real_space=.TRUE. ! Real profiles

    ! Calculate halo attributes
    rv=exp(find(log(m),hmod%log_m,log(hmod%rv),hmod%n,3,3,2))
    c=find(log(m),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c

    ! Field types
    nf=5
    ALLOCATE(fields(nf))
    fields(1)=field_matter
    fields(2)=field_cdm
    fields(3)=field_gas
    fields(4)=field_star
    fields(5)=field_electron_pressure

    ! Write file
    OPEN(7,file=outfile)
    DO i=1,n
       r=exp(progression(log(rmin),log(rmax),i,n))
       r=r*rv
       WRITE(7,*) r/rv, (win_type(real_space,fields(j),r,m,rv,rs,hmod,cosm)*rv**3, j=1,nf) ! rv**3 here is from r^2 dr in integral
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_profiles

  SUBROUTINE write_halo_transforms(m,hmod,cosm,outfile)

    !Writes out to file the Fourier transform of the halo density profiles
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile  
    REAL :: x, rv, c, rs, k, rhobar
    INTEGER :: i, j, nf
    INTEGER, ALLOCATABLE :: fields(:)

    REAL, PARAMETER :: xmin=1e-1      ! Minimum r/rv
    REAL, PARAMETER :: xmax=1e2       ! Maximum r/rv
    INTEGER, PARAMETER :: n=512       ! Number of points
    LOGICAL, PARAMETER :: rsp=.FALSE. ! Fourier profiles

    !Calculate halo attributes
    rv=exp(find(log(m),hmod%log_m,log(hmod%rv),hmod%n,3,3,2))
    c=find(log(m),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c

    ! Field types
    nf=5
    ALLOCATE(fields(nf))
    fields(1)=field_matter
    fields(2)=field_cdm
    fields(3)=field_gas
    fields(4)=field_star
    fields(5)=field_electron_pressure

    !Need mean density
    rhobar=comoving_matter_density(cosm)

    !Write file
    OPEN(7,file=outfile)
    DO i=1,n
       x=exp(progression(log(xmin),log(xmax),i,n))
       k=x/rv
       WRITE(7,*) x, (win_type(rsp,fields(j),k,m,rv,rs,hmod,cosm)*rhobar/m, j=1,nf)
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_transforms

  REAL FUNCTION delta_c(hmod,cosm)

    !Linear collapse density
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    a=hmod%a

    IF(hmod%idc==1) THEN
       !Fixed value
       delta_c=1.686
    ELSE IF(hmod%idc==2) THEN
       !From Nakamura & Suto (1997) LCDM fitting function
       delta_c=dc_NakamuraSuto(a,cosm)
    ELSE IF(hmod%idc==3) THEN
       !From Mead et al. (2015, 2016)
       delta_c=hmod%dc0+hmod%dc1*log(sigma(8.,a,cosm))
       delta_c=delta_c*(dc_NakamuraSuto(a,cosm)/dc0)
    ELSE IF(hmod%idc==4) THEN
       !From Mead (2017) fitting function
       delta_c=dc_Mead(a,cosm)
    ELSE IF(hmod%idc==5) THEN
       !From spheircal-collapse calculation
       delta_c=dc_spherical(a,cosm)
    ELSE
       WRITE(*,*) 'DELTA_C: idc:', hmod%idc
       STOP 'DELTA_C: Error, idc defined incorrectly'
    END IF

  END FUNCTION delta_c

  REAL FUNCTION Delta_v(hmod,cosm)

    ! Virialised overdensity
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    a=hmod%a

    IF(hmod%iDv==1) THEN
       ! Fixed value
       Delta_v=200.
    ELSE IF(hmod%iDv==2) THEN
       ! From Bryan & Norman (1998; arXiv:astro-ph/9710107) fitting functions
       Delta_v=Dv_BryanNorman(a,cosm)    
    ELSE IF(hmod%iDv==3) THEN
       ! From Mead et al. (2015, 2016)
       Delta_v=hmod%Dv0*Omega_m(a,cosm)**hmod%Dv1
    ELSE IF(hmod%iDv==4) THEN
       ! From Mead (2017) fitting function
       Delta_v=Dv_Mead(a,cosm)
    ELSE IF(hmod%iDv==5) THEN
       ! From spheircal-collapse calculation
       Delta_v=Dv_spherical(a,cosm)
    ELSE IF(hmod%iDv==6) THEN
       ! Lagrangian radius
       Delta_v=1.
    ELSE
       STOP 'DELTA_V: Error, iDv defined incorrectly'
    END IF

  END FUNCTION Delta_v

  REAL FUNCTION eta(hmod,cosm)

    ! Calculates the eta that comes into the bastardised one-halo term
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: eta0

    IF(hmod%ieta==1) THEN
       eta=0.
    ELSE IF(hmod%ieta==2) THEN
       ! From Mead et al. (2015; arXiv 1505.07833, 2016)
       IF(hmod%one_parameter_baryons) THEN
          eta0=0.98-hmod%As*0.12
       ELSE
          eta0=hmod%eta0
       END IF
       eta=eta0-hmod%eta1*(sigma(8.,hmod%a,cosm))
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION eta

  FUNCTION kstar(hmod,cosm)

    ! Calculates the one-halo damping wave number
    IMPLICIT NONE
    REAL :: kstar
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%ikstar==1) THEN
       ! Set to zero for the standard Poisson one-halo term
       kstar=0.
    ELSE IF(hmod%ikstar==2) THEN
       ! One-halo cut-off wavenumber from Mead et al. (2015, 2016)
       kstar=hmod%ks/hmod%sigv
    ELSE
       STOP 'KSTAR: Error, ihm defined incorrectly'
    END IF

  END FUNCTION kstar

  REAL FUNCTION As(hmod,cosm)

    ! Halo concentration pre-factor
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A
    
    IF(hmod%iAs==1) THEN
       ! Set to 4 for the standard Bullock value
       As=4.
    ELSE IF(hmod%iAs==2) THEN
       ! This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833, 2016)
       As=hmod%As
    ELSE
       STOP 'AS: Error, iconc defined incorrectly'
    END IF

    ! Now this is divided by 4 so as to be relative to the Bullock base result
    As=As/4.

  END FUNCTION As

  REAL FUNCTION fdamp(hmod,cosm)

    ! Calculates the linear-theory damping factor
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%i2hdamp==1) THEN
       ! Set to 0 for the standard linear theory two halo term
       fdamp=0.
    ELSE IF(hmod%i2hdamp==2) THEN
       ! Mead et al. (2015)
       fdamp=hmod%f0*hmod%sig8z**hmod%f1
    ELSE IF(hmod%i2hdamp==3) THEN
       ! Mead et al. (2016)
       fdamp=hmod%f0*hmod%sigv100**hmod%f1
    ELSE
       STOP 'FDAMP: Error, i2hdamp defined incorrectly'
    END IF

    ! Catches extreme values of fdamp that occur for ridiculous cosmologies
    IF(fdamp<fdamp_min) fdamp=0.
    IF(fdamp>fdamp_max) fdamp=fdamp_max

  END FUNCTION fdamp

  FUNCTION alpha_transition(hmod,cosm)

    ! Calculates the alpha to smooth the two- to one-halo transition
    IMPLICIT NONE
    REAL :: alpha_transition
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    ! To prevent compile-time warnings
    crap=cosm%A
    
    IF(hmod%itrans==2) THEN
       ! From Mead et al. (2015, 2016)   
       alpha_transition=hmod%alp0*hmod%alp1**hmod%neff       
    ELSE IF(hmod%itrans==4) THEN
       ! Specially for HMx, exponentiated Mead et al. (2016) result
       alpha_transition=(hmod%alp0*hmod%alp1**hmod%neff)**2.5
    ELSE
       alpha_transition=1.
    END IF

    ! Catches values of alpha that are crazy
    IF(alpha_transition<alpha_transition_min) alpha_transition=alpha_transition_min
    IF(alpha_transition>alpha_transition_max) alpha_transition=alpha_transition_max 

  END FUNCTION alpha_transition

  REAL FUNCTION HMx_alpha(m,hmod)

    ! TODO: Should Mp be M* ?
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, E
    REAL, PARAMETER :: Mp=1e14

    IF(hmod%fixed_HMx) THEN

       HMx_alpha=hmod%alpha*(m/Mp)**hmod%alphap

    ELSE

       A=hmod%A_alpha
       B=hmod%B_alpha
       C=hmod%C_alpha
       D=hmod%D_alpha
       E=hmod%E_alpha

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_alpha=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_alpha=(A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
       HMx_alpha=HMx_alpha/(1.+z) ! NEW: (1+z) accounts for previous wrong definition of temperature

    END IF

    IF(HMx_alpha<HMx_alpha_min) HMx_alpha=HMx_alpha_min
       
  END FUNCTION HMx_alpha
  
  REAL FUNCTION HMx_eps(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%fixed_HMx) THEN

       HMx_eps=hmod%eps

    ELSE

       A=hmod%A_eps
       B=hmod%B_eps
       C=hmod%C_eps
       D=hmod%D_eps

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_eps=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_eps=10**HMx_eps

    END IF
       
  END FUNCTION HMx_eps

  REAL FUNCTION HMx_Gamma(m,hmod)

    ! TODO: Should Mp be M* ?
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, E
    REAL, PARAMETER :: Mp=1e14

    IF(hmod%fixed_HMx) THEN

       HMx_Gamma=hmod%Gamma*(m/Mp)**hmod%Gammap

    ELSE

       A=hmod%A_Gamma
       B=hmod%B_Gamma
       C=hmod%C_Gamma
       D=hmod%D_Gamma
       E=hmod%E_Gamma

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_Gamma=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_Gamma=(A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E

    END IF

    !IF(HMx_Gamma<=1.) STOP 'HMx_GAMMA: Error, Gamma <= 1'
    IF(HMx_Gamma<HMx_Gamma_min) HMx_Gamma=HMx_Gamma_min
    IF(HMx_Gamma>HMx_Gamma_max) HMx_Gamma=HMx_Gamma_max
       
  END FUNCTION HMx_Gamma

  REAL FUNCTION HMx_M0(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D, E

    IF(hmod%fixed_HMx) THEN

       HMx_M0=hmod%M0

    ELSE

       A=hmod%A_M0
       B=hmod%B_M0
       C=hmod%C_M0
       D=hmod%D_M0
       E=hmod%E_M0

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_M0=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_M0=(A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
       HMx_M0=10**HMx_M0

    END IF
       
  END FUNCTION HMx_M0

  REAL FUNCTION HMx_Astar(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%fixed_HMx) THEN

       HMx_Astar=hmod%Astar

    ELSE

       A=hmod%A_Astar
       B=hmod%B_Astar
       C=hmod%C_Astar
       D=hmod%D_Astar

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_Astar=A*(1.+z)*T+B*(1.+z)+C*T+D

    END IF

    IF(HMx_Astar<HMx_Astar_min) HMx_Astar=HMx_Astar_min
       
  END FUNCTION HMx_Astar

  REAL FUNCTION HMx_Twhim(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D
    REAL :: zs(4), Ts(4)
    LOGICAL, PARAMETER :: simple_model=.TRUE.

    IF(hmod%fixed_HMx) THEN

       IF(simple_model) THEN
       
          HMx_Twhim=hmod%Twhim

       ELSE

          ! From my fitting
          
          zs(1)=0.0
          zs(2)=0.5
          zs(3)=1.0
          zs(4)=2.0
          
          Ts(1)=5.894
          Ts(2)=5.899
          Ts(3)=5.791
          Ts(4)=5.474
          
          z=hmod%z
          HMx_Twhim=10**find(z,zs,Ts,4,3,3,2)

       END IF

    ELSE

       A=hmod%A_Twhim
       B=hmod%B_Twhim
       C=hmod%C_Twhim
       D=hmod%D_Twhim

       z=hmod%z
       T=log10(hmod%Theat)
       !HMx_Twhim=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_Twhim=(A*(1.+z)**2+B*(1.+z)+C)*T+D
       HMx_Twhim=10**HMx_Twhim

    END IF
       
  END FUNCTION HMx_Twhim

  REAL FUNCTION HMx_cstar(m,hmod)

    ! TODO: Should Mp be M* ?
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL, PARAMETER :: Mp=1e14

    HMx_cstar=hmod%cstar*(m/Mp)**hmod%cstarp
       
  END FUNCTION HMx_cstar

  FUNCTION r_nl(hmod)

    ! Calculates R_nl where nu(R_nl)=1.
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: r_nl  

    IF(hmod%nu(1)>1.) THEN
       ! This catches some very strange values
       r_nl=hmod%rr(1)
    ELSE
       r_nl=exp(find(log(1.),log(hmod%nu),log(hmod%rr),hmod%n,3,3,2))
    END IF

  END FUNCTION r_nl

  SUBROUTINE allocate_HMOD(hmod)

    ! Allocates memory for the look-up tables
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    INTEGER :: n

    ! Number of entries in look-up table
    n=hmod%n

    ALLOCATE(hmod%log_m(n))
    ALLOCATE(hmod%zc(n),hmod%m(n),hmod%c(n),hmod%rv(n))
    ALLOCATE(hmod%nu(n),hmod%rr(n),hmod%sigf(n),hmod%sig(n))
    ALLOCATE(hmod%m500(n),hmod%r500(n),hmod%c500(n))
    ALLOCATE(hmod%m500c(n),hmod%r500c(n),hmod%c500c(n))
    ALLOCATE(hmod%m200(n),hmod%r200(n),hmod%c200(n))
    ALLOCATE(hmod%m200c(n),hmod%r200c(n),hmod%c200c(n))

    hmod%log_m=0.
    hmod%zc=0.
    hmod%m=0.
    hmod%c=0.
    hmod%rv=0.
    hmod%nu=0.
    hmod%rr=0.
    hmod%sigf=0.
    hmod%sig=0.

    hmod%m500=0.
    hmod%r500=0.
    hmod%c500=0.

    hmod%m500c=0.
    hmod%r500c=0.
    hmod%c500c=0.

    hmod%m200=0.
    hmod%r200=0.
    hmod%c200=0.

    hmod%m200c=0.
    hmod%r200c=0.
    hmod%c200c=0.

  END SUBROUTINE allocate_HMOD

  SUBROUTINE deallocate_HMOD(hmod)

    !Deallocates the look-up tables
    IMPLICIT NONE
    TYPE(halomod) :: hmod

    !Deallocates look-up tables
    DEALLOCATE(hmod%log_m)
    DEALLOCATE(hmod%zc,hmod%m,hmod%c,hmod%rv,hmod%nu,hmod%rr,hmod%sigf,hmod%sig)
    DEALLOCATE(hmod%m500,hmod%r500,hmod%c500,hmod%m500c,hmod%r500c,hmod%c500c)
    DEALLOCATE(hmod%m200,hmod%r200,hmod%c200,hmod%m200c,hmod%r200c,hmod%c200c)

    !Deallocate experimental window tables
    !DEALLOCATE(hmod%log_win,hmod%log_k)

  END SUBROUTINE deallocate_HMOD

  REAL FUNCTION Omega_stars(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_stars=rhobar(hmod%nu(1),hmod%large_nu,rhobar_star_integrand,hmod,cosm)
    Omega_stars=Omega_stars/comoving_critical_density(hmod%a,cosm)

  END FUNCTION Omega_stars

  SUBROUTINE init_galaxies(hmod,cosm)

    !Calcuate the number densities of galaxies
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu_min

    nu_min=nu_M(hmod%mgal,hmod,cosm)
    hmod%n_c=rhobar(nu_min,hmod%large_nu,rhobar_central_integrand,hmod,cosm)
    hmod%n_s=rhobar(nu_min,hmod%large_nu,rhobar_satellite_integrand,hmod,cosm)
    hmod%n_g=hmod%n_c+hmod%n_s
    IF(verbose_galaxies) THEN
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of central galaxies [(Mpc/h)^-3]:', REAL(hmod%n_c)
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of satellite galaxies [(Mpc/h)^-3]:', REAL(hmod%n_s)
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of all galaxies [(Mpc/h)^-3]:', REAL(hmod%n_g)
       WRITE(*,*)
    END IF

    hmod%has_galaxies=.TRUE.
  
  END SUBROUTINE init_galaxies

  SUBROUTINE init_HI(hmod,cosm)

    ! Gets the background HI density by integrating the HI mass function
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu_min, nu_max

    nu_min=hmod%nu(1)
    nu_max=hmod%nu(hmod%n)
    hmod%rho_HI=rhobar(nu_min,hmod%large_nu,rhobar_HI_integrand,hmod,cosm)
    IF(verbose_HI) THEN
       WRITE(*,*) 'INIT_HI: z:', hmod%z
       WRITE(*,*) 'INIT_HI: HI density [log10(rho/(Msun/h)/(Mpc/h)^3)]:', REAL(log10(hmod%rho_HI))
       WRITE(*,*) 'INIT_HI: Omega_HI(z):', hmod%rho_HI/comoving_critical_density(hmod%a,cosm)
       WRITE(*,*) 'INIT_HI: Omega_HI(z) relative to z=0 critical density:', hmod%rho_HI/comoving_critical_density(1.,cosm)
       WRITE(*,*) 'INIT_HI: rho_HI / rho_matter:', hmod%rho_HI/comoving_matter_density(cosm)
       WRITE(*,*)
    END IF

    hmod%has_HI=.TRUE.

  END SUBROUTINE init_HI
  
  REAL FUNCTION nu_R(R,hmod,cosm)

    ! Calculates nu(R) where R is the comoving Lagrangian halo radius
    IMPLICIT NONE
    REAL, INTENT(IN) :: R
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    a=hmod%a
    nu_R=delta_c(hmod,cosm)/sigma(R,a,cosm)

  END FUNCTION nu_R

  REAL FUNCTION nu_M(M,hmod,cosm)

    ! Calculates nu(M) where M is the halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: M
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: R

    R=radius_m(M,cosm)
    nu_M=nu_R(R,hmod,cosm)
    
  END FUNCTION nu_M

  FUNCTION M_nu(nu,hmod)

    !Calculates M(nu) where M is the halo mass and nu is the peak height
    IMPLICIT NONE
    REAL :: M_nu
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    M_nu=exp(find(nu,hmod%nu,hmod%log_m,hmod%n,3,3,2))

  END FUNCTION M_nu

  REAL FUNCTION rhobar_central_integrand(nu,hmod,cosm)

    !Integrand for the number density of central galaxies
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Could remove
    REAL :: M, crap

    crap=cosm%A

    M=M_nu(nu,hmod)    
    rhobar_central_integrand=N_centrals(M,hmod)*g_nu(nu,hmod)/M
    
  END FUNCTION rhobar_central_integrand

  REAL FUNCTION rhobar_satellite_integrand(nu,hmod,cosm)

    !Integrand for the number density of satellite galaxies
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Could remove
    REAL :: M, crap

    crap=cosm%A

    M=M_nu(nu,hmod)    
    rhobar_satellite_integrand=N_satellites(M,hmod)*g_nu(nu,hmod)/M
    
  END FUNCTION rhobar_satellite_integrand

  REAL FUNCTION rhobar_star_integrand(nu,hmod,cosm)

    !Integrand for the matter density of stars
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: M

    M=M_nu(nu,hmod)    
    rhobar_star_integrand=halo_star_fraction(M,hmod,cosm)*g_nu(nu,hmod)
    
  END FUNCTION rhobar_star_integrand

  REAL FUNCTION rhobar_HI_integrand(nu,hmod,cosm)

    !Integrand for the HI mass density
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: M

    M=M_nu(nu,hmod)    
    rhobar_HI_integrand=halo_HI_fraction(M,hmod,cosm)*g_nu(nu,hmod)
    
  END FUNCTION rhobar_HI_integrand

  REAL FUNCTION rhobar(nu_min,nu_max,integrand,hmod,cosm)

    !Calculate the mean density of a tracer
    !Integrand here is a function of mass, i.e. I(M); R = rho * Int I(M)dM
    !TODO: This function is very slow and should be accelerated somehow
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu_min, nu_max
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    INTERFACE
       REAL FUNCTION integrand(M,hmod,cosm)
         IMPORT :: halomod, cosmology
         REAL, INTENT(IN) :: M
         TYPE(halomod), INTENT(INOUT) :: hmod
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION integrand
    END INTERFACE
    
    rhobar=comoving_matter_density(cosm)*integrate_hmod_cosm_exp(log(nu_min),log(nu_max),integrand,hmod,cosm,hmod%acc_HMx,3)
    
  END FUNCTION rhobar
    
  FUNCTION one_halo_amplitude(hmod,cosm)

    !Calculates the amplitude of the shot-noise plateau of the one-halo term
    IMPLICIT NONE
    REAL :: one_halo_amplitude
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: integrand(hmod%n), g, m
    INTEGER :: i

    !Calculates the value of the integrand at all nu values!
    DO i=1,hmod%n
       g=g_nu(hmod%nu(i),hmod)
       m=hmod%m(i)
       integrand(i)=g*m
    END DO

    one_halo_amplitude=integrate_table(hmod%nu,integrand,hmod%n,1,hmod%n,1)/comoving_matter_density(cosm)

  END FUNCTION one_halo_amplitude

  REAL FUNCTION mass_function(m,hmod,cosm)

    ! Calculates the halo mass function, what I call n(M)
    IMPLICIT NONE
    REAL, INTENT(IN) :: m ! Halo mass
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: nu, dnu_dm

    nu=nu_M(m,hmod,cosm)
    dnu_dm=derivative_table(m,hmod%m,hmod%nu,hmod%n,3,3)
    mass_function=comoving_matter_density(cosm)*g_nu(nu,hmod)*dnu_dm/m
    
  END FUNCTION mass_function

  SUBROUTINE convert_mass_definitions(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rhom, rhoc, Dv, a

    ! Scale factor
    a=hmod%a

    ! Get the densities
    rhom=comoving_matter_density(cosm)
    rhoc=comoving_critical_density(a,cosm)
    Dv=Delta_v(hmod,cosm)

    ! Calculate Delta = 200, 500 and Delta_c = 200, 500 quantities
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,1.,hmod%r500,hmod%c500,hmod%m500,500.,1.,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,1.,hmod%r200,hmod%c200,hmod%m200,200.,1.,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,rhom,hmod%r500c,hmod%c500c,hmod%m500c,500.,rhoc,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,rhom,hmod%r200c,hmod%c200c,hmod%m200c,200.,rhoc,hmod%n)

    hmod%has_mass_conversions=.TRUE.
    
  END SUBROUTINE convert_mass_definitions

  SUBROUTINE convert_mass_definition(r1,c1,m1,D1,rho1,r2,c2,m2,D2,rho2,n)

    !Converts mass definition from Delta_1 rho_1 overdense to Delta_2 rho_2 overdense
    IMPLICIT NONE
    REAL, INTENT(IN) :: r1(n) ! Array of initial virial radii [Mpc/h]
    REAL, INTENT(IN) :: c1(n) ! Array of initial halo concentration
    REAL, INTENT(IN) :: m1(n) ! Array of initial halo mass [Msun/h]
    REAL, INTENT(IN) :: D1 ! Initial halo overdensity definition (e.g., 200, 500)
    REAL, INTENT(IN) :: rho1 ! Initial halo overdensity defintion (critical or mass)
    REAL, INTENT(OUT) :: r2(n) ! Output array of virial radii with new definition [Mpc/h]
    REAL, INTENT(OUT) :: c2(n) ! Output array of halo concentration with new definition
    REAL, INTENT(OUT) :: m2(n) ! Output array of halo mass with new definition [Msun/h]
    REAL, INTENT(IN) :: D2 ! Final halo overdensity definition (e.g., 200, 500)
    REAL, INTENT(IN) :: rho2 ! Final halo overdensity defintion (critical or mass)
    INTEGER, INTENT(IN) :: n  ! Number of entries in tables
    REAL :: f(n)
    REAL :: rmin, rmax, rs
    REAL, ALLOCATABLE :: r(:)
    INTEGER :: i, j

    IF(verbose_convert_mass) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Converting mass definitions:'
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Initial overdensity:', D1
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Final overdensity:', D2
    END IF

    !Make an array of general 'r' values for the solution later
    rmin=r1(1)/10. !Should be sufficient for reasonable cosmologies
    rmax=r1(n)*10. !Should be sufficient for reasonable cosmologies
    CALL fill_array(log(rmin),log(rmax),r,n) !Necessary to log space
    r=exp(r)
    
    IF(verbose_convert_mass) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: rmin:', rmin
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: rmax:', rmax
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: nr:', n
    END IF
    
    !Now use the find algorithm to invert L(r_i)=R(r_j) so that r_j=R^{-1}[L(r_i)]
    IF(verbose_convert_mass) THEN
       WRITE(*,*) '========================================================================================================'
       WRITE(*,*) '         M_old         rv_old          c_old    M_new/M_old    r_new/r_old    c_new/c_old  M`_new/M`_old' 
       WRITE(*,*) '========================================================================================================'
    END IF
    DO i=1,n

       !Calculate the halo scale radius
       rs=r1(i)/c1(i)

       !Fill up the f(r) table which needs to be solved for R | f(R)=0 
       DO j=1,n
          f(j)=r1(i)**3*rho1*D1/NFW_factor(r1(i)/rs)-r(j)**3*rho2*D2/NFW_factor(r(j)/rs) !NFW
          !f(j)=r1(i)**2*rho1*D1-r(j)**2*rho2*D2 !Isothemral sphere
       END DO
       
       !First find the radius R | f(R)=0; I am fairly certain that I can use log on 'r' here
       r2(i)=exp(find_solve(0.,log(r),f,n))

       !Now do the concentration and mass conversions
       c2(i)=r2(i)/rs
       m2(i)=m1(i)*(rho2*D2*r2(i)**3)/(rho1*D1*r1(i)**3)
       !m2(i)=m1(i)*NFW_factor(c2(i))/NFW_factor(c1(i))

       IF(verbose_convert_mass) THEN
          WRITE(*,fmt='(ES15.7,6F15.7)') m1(i), r1(i), c1(i), m2(i)/m1(i), r2(i)/r1(i), c2(i)/c1(i), NFW_factor(c2(i))/NFW_factor(c1(i))
       END IF

    END DO

    IF(verbose_convert_mass) WRITE(*,*) '========================================================================================================'

    IF(verbose_convert_mass) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Done'
       WRITE(*,*)
       STOP
    END IF

  END SUBROUTINE convert_mass_definition

  REAL FUNCTION NFW_factor(x)

    ! The NFW 'mass' factor that crops up all the time
    ! This is X(c) in M(r) = M X(r/rs) / X(c)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    NFW_factor=log(1.+x)-x/(1.+x)

  END FUNCTION NFW_factor

  FUNCTION radius_m(m,cosm)

    ! The comoving radius corresponding to mass M in a homogeneous universe
    IMPLICIT NONE
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(INOUT) :: cosm

    radius_m=(3.*m/(4.*pi*comoving_matter_density(cosm)))**(1./3.)

  END FUNCTION radius_m

  FUNCTION virial_radius(m,hmod,cosm)

    ! The comoving halo virial radius
    IMPLICIT NONE
    REAL :: virial_radius
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    virial_radius=(3.*m/(4.*pi*comoving_matter_density(cosm)*Delta_v(hmod,cosm)))**(1./3.)

  END FUNCTION virial_radius

  REAL FUNCTION effective_index(hmod,cosm)

    ! Power spectrum slope a the non-linear scale
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod

    ! Numerical differentiation to find effective index at collapse
    effective_index=-3.-derivative_table(log(hmod%rnl),log(hmod%rr),log(hmod%sig**2),hmod%n,3,3)

    ! For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
    ! In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    IF(effective_index<cosm%n-4.) effective_index=cosm%n-4.
    IF(effective_index>cosm%n)      effective_index=cosm%n

  END FUNCTION effective_index

  SUBROUTINE fill_halo_concentration(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: mnl, m, zc, z
    INTEGER :: i

    !iconc = 1: Full Bullock et al. (2001)
    !iconc = 2: Simple Bullock et al. (2001)
    !iconc = 3: Duffy et al. (2008): mean
    !iconc = 4: Duffy et al. (2008): virial

    ! Get the redshift
    z=hmod%z

    !Any initialisation for the c(M) relation goes here
    IF(hmod%iconc==1) THEN
       !Fill the collapse z look-up table
       CALL zcoll_Bullock(z,hmod,cosm)
    ELSE IF(hmod%iconc==2) THEN
       mnl=hmod%mnl
    END IF

    !Fill concentration-mass for all halo masses
    DO i=1,hmod%n

       !Halo mass
       m=hmod%m(i)

       !Choose concentration-mass relation
       IF(hmod%iconc==1) THEN
          zc=hmod%zc(i)
          hmod%c(i)=conc_Bullock(z,zc)
       ELSE IF(hmod%iconc==2) THEN         
          hmod%c(i)=conc_Bullock_simple(m,mnl)
       ELSE IF(hmod%iconc==3) THEN
          hmod%c(i)=conc_Duffy_mean(m,z)
       ELSE IF(hmod%iconc==4) THEN
          hmod%c(i)=conc_Duffy_virial(m,z)
       ELSE IF(hmod%iconc==5) THEN
          hmod%c(i)=conc_Duffy_relaxed200(m,z)
       ELSE
          STOP 'FILL_HALO_CONCENTRATION: Error, iconc specified incorrectly'
       END IF

       ! Rescale halo concentrations via the 'A' HMcode parameter
       hmod%c(i)=hmod%c(i)*As(hmod,cosm)   

    END DO

    ! Dolag2004 prescription for adding DE dependence
    IF(hmod%iDolag==2 .OR. hmod%iDolag==3) CALL Dolag_correction(hmod,cosm)

    ! Rescale the concentration-mass relation for gas the epsilon parameter
    ! This only rescales the concentrations of haloes that *contain* substantial amounts of gas
    DO i=1,hmod%n
       m=hmod%m(i)
       hmod%c(i)=hmod%c(i)*gas_correction(m,hmod,cosm)
    END DO
    
  END SUBROUTINE fill_halo_concentration

  REAL FUNCTION gas_correction(m,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    gas_correction=(1.+(HMx_eps(hmod)-1.)*halo_boundgas_fraction(m,hmod,cosm)/(cosm%Om_b/cosm%Om_m))
    
  END FUNCTION gas_correction

  SUBROUTINE Dolag_correction(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: g_LCDM, g_wCDM, ainf, f
    TYPE(cosmology) :: cosm_LCDM

    ! The 'infinite' scale factor
    ainf=scale_factor_z(zinf_Dolag)

    ! Save the growth function in the current cosmology
    g_wCDM=grow(ainf,cosm)

    ! Make a flat LCDM cosmology and calculate growth
    cosm_LCDM=cosm
    cosm_LCDM%iw=1
    cosm_LCDM%w=-1.
    cosm_LCDM%wa=0.
    cosm_LCDM%Om_w=0.
    cosm_LCDM%Om_v=1.-cosm%Om_m ! Added this so that 'making a LCDM cosmology' works for curved models.
    cosm_LCDM%verbose=.FALSE.
    CALL init_cosmology(cosm_LCDM) ! This is **essential**
    
    g_LCDM=growth_Linder(ainf,cosm_LCDM)
    f=g_wCDM/g_LCDM
    
    !WRITE(*,*) 'DOLAG_CORRECTION:', f
    !STOP

    ! Changed this to a power of 1.5, which produces more accurate results for extreme DE
    IF(hmod%iDolag==2) THEN
       hmod%c=hmod%c*f
    ELSE IF(hmod%iDolag==3) THEN
       hmod%c=hmod%c*f**1.5
    END IF
    
  END SUBROUTINE Dolag_correction

  REAL FUNCTION conc_Bullock(z,zc)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    REAL, INTENT(IN) :: zc

    REAL, PARAMETER :: A=4. !Pre-factor for Bullock relation
    
    conc_Bullock=A*(1.+zc)/(1.+z)

  END FUNCTION conc_Bullock

  SUBROUTINE zcoll_Bullock(z,hmod,cosm)

    !This fills up the halo collapse redshift table as per Bullock relations   
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: dc
    REAL :: af, zf, RHS, a, growz
    INTEGER :: i

    REAL, PARAMETER :: f=0.01**(1./3.) ! This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)

    a=scale_factor_z(z)

    ! Fills up a table for sigma(fM) for Bullock c(m) relation    
    IF(hmod%iconc==1) THEN
       DO i=1,hmod%n
          hmod%sigf(i)=sigma(hmod%rr(i)*f,a,cosm)
       END DO
    END IF

    ! Do numerical inversion
    DO i=1,hmod%n

       ! I don't think this is really consistent with dc varying as a function of z
       ! but the change will *probably* be *very* small
       dc=delta_c(hmod,cosm)

       RHS=dc*grow(a,cosm)/hmod%sigf(i)
       
       ! growz=find(a,af_tab,grow_tab,cosm%ng,3,3,2)
       ! growz=exp(find(log(a),cosm%a_growth,cosm%growth,cosm%n_growth,3,3,2))
       growz=grow(a,cosm)

       IF(RHS>growz) THEN
          zf=z
       ELSE
          af=exp(find(log(RHS),cosm%log_growth,cosm%log_a_growth,cosm%n_growth,3,3,2))
          zf=redshift_a(af)
       END IF

       hmod%zc(i)=zf

    END DO

  END SUBROUTINE zcoll_Bullock

  FUNCTION conc_Bullock_simple(m,mstar)

    !The simple concentration-mass relation from Bullock et al. (2001; astro-ph/9908159v3 equation 18)
    IMPLICIT NONE
    REAL :: conc_Bullock_simple
    REAL, INTENT(IN) :: m, mstar

    conc_Bullock_simple=9.*(m/mstar)**(-0.13)
    
  END FUNCTION conc_Bullock_simple

  REAL FUNCTION conc_Duffy_mean(m,z)

    !Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z
    
    REAL, PARAMETER :: m_piv=2e12 !Pivot mass in Msun/h
    REAL, PARAMETER :: A=10.14
    REAL, PARAMETER :: B=-0.081
    REAL, PARAMETER :: C=-1.01

    !Equation (4) in 0804.2486, parameters from 10th row of Table 1
    conc_Duffy_mean=A*(m/m_piv)**B*(1.+z)**C
    
  END FUNCTION conc_Duffy_mean

  REAL FUNCTION conc_Duffy_virial(m,z)

    !Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z
    
    REAL, PARAMETER :: m_piv=2e12 !Pivot mass in Msun/h
    REAL, PARAMETER :: A=7.85
    REAL, PARAMETER :: B=-0.081
    REAL, PARAMETER :: C=-0.71

    !Equation (4) in 0804.2486, parameters from 6th row of Table 1
    conc_Duffy_virial=A*(m/m_piv)**B*(1.+z)**C
    
  END FUNCTION conc_Duffy_virial

  REAL FUNCTION conc_Duffy_relaxed200(m,z)

    !Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z
    
    REAL, PARAMETER :: m_piv=2e12 !Pivot mass in Msun/h
    REAL, PARAMETER :: A=6.71
    REAL, PARAMETER :: B=-0.091
    REAL, PARAMETER :: C=-0.44

    !Equation (4) in 0804.2486, parameters from 4th row of Table 1
    conc_Duffy_relaxed200=A*(m/m_piv)**B*(1.+z)**C
    
  END FUNCTION conc_Duffy_relaxed200

  REAL FUNCTION mass_r(r,cosm)

    !Calcuates the mass contains in a sphere of comoving radius 'r' in a homogeneous universe
    IMPLICIT NONE
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !Relation between mean cosmological mass and radius
    mass_r=(4.*pi/3.)*comoving_matter_density(cosm)*(r**3)

  END FUNCTION mass_r

  REAL FUNCTION win_type(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Selects the halo profile type
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(itype==field_dmonly) THEN
       ! Overdensity if all the matter were CDM
       win_type=win_DMONLY(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_matter) THEN
       ! Matter overdensity (sum of CDM, gas, stars)
       win_type=win_total(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_cdm) THEN
       ! CDM overdensity
       win_type=win_CDM(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_gas) THEN
       ! All gas, both bound and free overdensity
       win_type=win_gas(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_star) THEN
       ! Stellar overdensity
       win_type=win_star(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_bound_gas) THEN
       ! Bound gas overdensity
       win_type=win_boundgas(real_space,itype,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_free_gas) THEN
       ! Free gas overdensity
       win_type=win_freegas(real_space,itype,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_electron_pressure) THEN
       ! Electron pressure
       win_type=win_electron_pressure(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_void) THEN
       ! Void
       win_type=win_void(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_compensated_void) THEN
       ! Compensated void
       win_type=win_compensated_void(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_central_galaxies) THEN
       ! Central galaxies
       win_type=win_centrals(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_satellite_galaxies) THEN
       ! Satellite galaxies
       win_type=win_satellites(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_galaxies) THEN
       ! All galaxies
       win_type=win_galaxies(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_HI) THEN
       ! Neutral hydrogen - HI
       win_type=win_HI(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(itype==field_cold_gas) THEN
       ! Cold gas
       win_type=win_coldgas(real_space,itype,k,m,rv,rs,hmod,cosm)
    ELSE
       WRITE(*,*) 'WIN_TYPE: itype:', itype
       STOP 'WIN_TYPE: Error, itype not specified correclty' 
    END IF

  END FUNCTION win_type

  REAL FUNCTION win_total(real_space,k,m,rv,rs,hmod,cosm)

    ! The halo profile of all the matter
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    win_total=win_CDM(real_space,k,m,rv,rs,hmod,cosm)+win_gas(real_space,k,m,rv,rs,hmod,cosm)+win_star(real_space,k,m,rv,rs,hmod,cosm)

  END FUNCTION win_total

  REAL FUNCTION win_DMONLY(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for all matter under the assumption that it is all CDM
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, c, rss, p1, p2

    ! Set additional halo parameters to zero
    p1=0.
    p2=0.

    rmin=0.
    rmax=rv

    IF(hmod%halo_DMONLY==1) THEN
       ! Analytical NFW
       irho=5 
    ELSE IF(hmod%halo_DMONLY==2) THEN
       ! Non-analyical NFW
       irho=4 
    ELSE IF(hmod%halo_DMONLY==3) THEN
       ! Tophat
       irho=2 
    ELSE IF(hmod%halo_DMONLY==4) THEN
       ! Delta function
       irho=0 
    ELSE IF(hmod%halo_DMONLY==5) THEN
       ! Cored NFW
       irho=24
       p1=hmod%rcore
    ELSE IF(hmod%halo_DMONLY==6) THEN
       ! Isothermal
       irho=1
    ELSE IF(hmod%halo_DMONLY==7) THEN
       ! Shell
       irho=28
       !STOP 'WIN_DMONLY: Error, shell profile is not implemented yet'
    ELSE
       STOP 'WIN_DMONLY: Error, halo_DMONLY specified incorrectly'
    END IF   

    ! Force it to use the gravity-only concentration relation (unapply gas correction)
    ! TODO: This is super ugly and should be improved somehow; also is unncessary calculations
    ! TODO: Somehow should store modified and unmodified halo concentrations
    c=rv/rs
    c=c/gas_correction(m,hmod,cosm)
    rss=rv/c
    !rss=rs

    IF(real_space) THEN
       r=k
       win_DMONLY=rho(r,rmin,rmax,rv,rss,p1,p2,irho)
       win_DMONLY=win_DMONLY/normalisation(rmin,rmax,rv,rss,p1,p2,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_DMONLY=m*win_norm(k,rmin,rmax,rv,rss,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_DMONLY

  REAL FUNCTION win_CDM(real_space,k,m,rv,rs,hmod,cosm)

    ! The halo profile for CDM
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2

    ! Default maximum and minimum radii
    rmin=0.
    rmax=rv

    ! Default additional halo parameters
    p1=0.
    p2=0.

    IF(hmod%halo_CDM==1) THEN      
       irho=5 ! Analytical NFW
    ELSE
       STOP 'WIN_CDM: Error, halo_CDM specified incorrectly'
    END IF   

    IF(real_space) THEN
       r=k
       win_CDM=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_CDM=win_CDM/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_CDM=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

    win_CDM=halo_CDM_fraction(m,hmod,cosm)*win_CDM

  END FUNCTION win_CDM

  REAL FUNCTION win_gas(real_space,k,m,rv,rs,hmod,cosm)

    !Halo profile for gas
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    INTEGER, PARAMETER :: itype_gas=2

    win_gas=win_boundgas(real_space,itype_gas,k,m,rv,rs,hmod,cosm)+win_freegas(real_space,itype_gas,k,m,rv,rs,hmod,cosm)+win_coldgas(real_space,itype_gas,k,m,rv,rs,hmod,cosm)

  END FUNCTION win_gas

  REAL FUNCTION win_boundgas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for the bound gas component
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rho0, T0, r, a
    REAL :: rmin, rmax, p1, p2
    INTEGER :: irho_density, irho_electron_pressure

    ! Initially set the halo parameters to zero
    p1=0.
    p2=0.

    ! Set maximum and minimum integration radius
    rmin=0.
    rmax=rv

    IF(hmod%halo_boundgas==1 .OR. hmod%halo_boundgas==3) THEN
       
       ! Komatsu & Seljak (2001) profile
       IF(hmod%halo_boundgas==1) THEN
          ! Simplified KS model
          irho_density=11
          irho_electron_pressure=13
       ELSE IF(hmod%halo_boundgas==3) THEN
          ! Full KS model
          irho_density=21
          irho_electron_pressure=23
       END IF    
       p1=HMx_Gamma(m,hmod)
       
    ELSE IF(hmod%halo_boundgas==2) THEN
       
       ! Set cored isothermal profile with beta=2/3 
       irho_density=6 
       irho_electron_pressure=irho_density ! okay to use density for electron pressure because temperature is constant
       
    ELSE       
       STOP 'WIN_BOUNDGAS: Error, halo_boundgas not specified correctly'
    END IF

    IF(itype==0 .OR. itype==2 .OR. itype==4) THEN

       ! Density profile of bound gas
       IF(real_space) THEN
          r=k
          win_boundgas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_density)
          win_boundgas=win_boundgas/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
       ELSE
          ! Properly normalise and convert to overdensity
          win_boundgas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_density)/comoving_matter_density(cosm)
       END IF

       win_boundgas=halo_boundgas_fraction(m,hmod,cosm)*win_boundgas

    ELSE IF(itype==6) THEN

       ! Electron pressure profile of bound gas
       IF(real_space) THEN
          r=k
          win_boundgas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
       ELSE
          ! The electron pressure window is T(r) x rho_e(r), we want unnormalised, so multiply through by normalisation
          ! TODO: Can I make the code more efficient here by having an unnorm window function?
          win_boundgas=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)*normalisation(rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
          !win_boundgas=winint(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
       END IF

       ! Calculate the value of the density profile prefactor and change units from cosmological to SI
       rho0=m*halo_boundgas_fraction(m,hmod,cosm)/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
       rho0=rho0*msun/mpc/mpc/mpc ! Overflow with REAL*4 if you use mpc**3
       rho0=rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

       ! Calculate the value of the temperature prefactor [K]
       a=hmod%a
       T0=HMx_alpha(m,hmod)*virial_temperature(m,rv,hmod%a,cosm)

       ! Convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
       win_boundgas=win_boundgas*(rho0/(mp*cosm%mue))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
       win_boundgas=win_boundgas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
       win_boundgas=win_boundgas*cosm%mue/cosm%mup ! Convert from total thermal pressure to electron pressure

    ELSE

       STOP 'WIN_BOUNDGAS: Error, itype not specified correctly'

    END IF

  END FUNCTION win_boundgas

  REAL FUNCTION win_freegas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for the free gas component
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: re, rmin, rmax, r, A, rho0, rhov, T0, p1, p2, beta, c, thing, m0
    INTEGER :: irho_density, irho_electron_pressure

    ! Enable to force the electron pressure to be matched at the virial radius
    ! This is enabled by default for some halo gas/pressure models
    LOGICAL :: match_electron_pressure=.FALSE.

       ! Set the halo 'parameter' variables to zero initially
       p1=0.
       p2=0.

       ! Set default min and max radii
       rmin=0.
       rmax=rv

       IF(halo_freegas_fraction(m,hmod,cosm)==0.) THEN

          ! Sometimes the freegas fraction will be zero, in which case this avoids problems
          win_freegas=0.

       ELSE

          IF(hmod%halo_freegas==1) THEN

             ! Simple isothermal model, motivated by constant velocity and rate expulsion
             irho_density=1
             irho_electron_pressure=irho_density ! Okay because T is constant
             rmin=0.
             rmax=2.*rv

          ELSE IF(hmod%halo_freegas==2) THEN

             ! Ejected gas model from Schneider (2015)
             irho_density=10
             irho_electron_pressure=irho_density ! Okay because T is constant
             rmin=rv
             re=rv
             p1=re
             rmax=15.*re ! Needs to be such that integral converges (15rf seems okay)

          ELSE IF(hmod%halo_freegas==3) THEN

             ! Now do isothermal shell connected to the KS profile continuously
             irho_density=16
             irho_electron_pressure=irho_density ! Okay because T is constant

             ! Isothermal model with continuous link to KS
             rhov=win_boundgas(.TRUE.,2,rv,m,rv,rs,hmod,cosm) ! This is the value of the density at the halo boundary for the bound gas           
             A=rhov/rho(rv,0.,rv,rv,rs,p1,p2,irho_density) ! This is A, as in A/r^2

             rmin=rv
             rmax=rv+halo_freegas_fraction(m,hmod,cosm)/(4.*pi*A) ! This ensures density continuity and mass conservation

             c=10. ! How many times larger than the virial radius can the gas cloud go?          
             IF(rmax>c*rv) rmax=c*rv ! This needs to be set otherwise get huge decrement in gas power at large scales
             match_electron_pressure=.TRUE. ! Match the electron pressure at the boundary

          ELSE IF(hmod%halo_freegas==4) THEN

             ! Ejected gas is a continuation of the KS profile
             irho_density=11 ! KS
             irho_electron_pressure=13 ! KS
             rmin=rv
             rmax=2.*rv
             p1=HMx_Gamma(m,hmod)

          ELSE IF(hmod%halo_freegas==5) THEN

             m0=1e14

             IF(m<m0) THEN

                irho_density=0
                irho_electron_pressure=irho_density
                rmin=0.
                rmax=rv

             ELSE

                ! Set the density profile to be the power-law profile
                irho_density=17
                irho_electron_pressure=irho_density ! Not okay

                ! Calculate the KS index at the virial radius
                c=rv/rs
                beta=(c-(1.+c)*log(1.+c))/((1.+c)*log(1.+c))
                beta=beta/(HMx_Gamma(m,hmod)-1.) ! This is the power-law index at the virial radius for the KS gas profile
                p1=beta
                !WRITE(*,*) 'Beta:', beta, log10(m)
                IF(beta<=-3.) beta=-2.9 ! If beta<-3 then there is only a finite amount of gas allowed in the free component

                ! Calculate the density at the boundary of the KS profile
                rhov=win_boundgas(.TRUE.,2,rv,m,rv,rs,hmod,cosm)
                !WRITE(*,*) 'rho_v:', rhov

                ! Calculate A as in rho(r)=A*r**beta
                A=rhov/rho(rv,0.,rv,rv,rs,p1,p2,irho_density)
                !WRITE(*,*) 'A:', A

                ! Set the minimum radius for the power-law to be the virial radius
                rmin=rv
                !WRITE(*,*) 'rmin:', rmin

                ! Set the maximum radius so that it joins to KS profile seamlessly
                thing=(beta+3.)*halo_freegas_fraction(m,hmod,cosm)/(4.*pi*A)+(rhov*rv**3)/A
                !WRITE(*,*) 'thing:', thing
                IF(thing>0.) THEN
                   ! This then fixes the condition of contiunity in amplitude and gradient
                   rmax=thing**(1./(beta+3.))
                ELSE
                   ! If there are no sohmodions then fix to 10rv and accept discontinuity
                   ! There may be no sohmodion if there is a lot of free gas and if beta<-3
                   rmax=10.*rv
                END IF
                !WRITE(*,*) 'rmax 2:', rmax

             END IF

          ELSE IF(hmod%halo_freegas==6) THEN

             ! Cubic profile
             rmin=rv
             rmax=3.*rv
             irho_density=18
             irho_electron_pressure=irho_density

          ELSE IF(hmod%halo_freegas==7) THEN

             ! Smooth profile (rho=0)
             rmin=0.
             rmax=rv
             irho_density=19
             irho_electron_pressure=irho_density

          ELSE IF(hmod%halo_freegas==8) THEN

             ! Delta function
             rmin=0.
             rmax=rv
             irho_density=0
             irho_electron_pressure=irho_density

          ELSE
             STOP 'WIN_FREEGAS: Error, halo_freegas specified incorrectly'
          END IF

          ! Density profile
          IF(itype==0 .OR. itype==2 .OR. itype==5) THEN

             ! Density profile of free gas
             IF(real_space) THEN
                r=k
                win_freegas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_density)
                win_freegas=win_freegas/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
             ELSE
                ! Properly normalise and convert to overdensity
                win_freegas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_density)/comoving_matter_density(cosm)
             END IF

             win_freegas=halo_freegas_fraction(m,hmod,cosm)*win_freegas

          ! Electron pressure profile
          ELSE IF(itype==6) THEN

             ! If we are applying a pressure-matching condition
             IF(match_electron_pressure) THEN

                STOP 'WIN_FREEGAS: Check the units and stuff here *very* carefully'

                r=k
                IF(r>rmin .AND. r<rmax) THEN
                   ! Only works for isothermal profile
                   win_freegas=win_boundgas(.TRUE.,6,rv,m,rv,rs,hmod,cosm)*(r/rv)**(-2)
                ELSE
                   win_freegas=0.
                END IF

             ELSE

                ! Electron pressure profile of free gas
                IF(real_space) THEN
                   r=k
                   win_freegas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
                ELSE  
                   win_freegas=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)*normalisation(rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
                END IF

                ! Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
                rho0=m*halo_freegas_fraction(m,hmod,cosm)/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density) ! rho0 in [(Msun/h)/(Mpc/h)^3]
                rho0=rho0*msun/Mpc/Mpc/Mpc ! Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
                rho0=rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

                ! This is the total thermal pressure of the WHIM
                T0=HMx_Twhim(hmod) ! [K]

                ! Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
                win_freegas=win_freegas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
                win_freegas=win_freegas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
                win_freegas=win_freegas*cosm%mue/cosm%mup ! Convert from total thermal pressure to electron pressure

             END IF

          ELSE

             STOP 'WIN_FREEGAS: Error, itype not specified correctly'

          END IF

       END IF

!!$    END IF

  END FUNCTION win_freegas

  REAL FUNCTION win_coldgas(real_space,itype,k,m,rv,rs,hmod,cosm)

    ! Halo profile for the cold gas component
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r, rmin, rmax, p1, p2
    INTEGER :: irho

    ! Initially set the halo parameters to zero
    p1=0.
    p2=0.

    ! Set maximum and minimum integration radius
    rmin=0.
    rmax=rv

    IF(hmod%halo_coldgas==1) THEN
        ! Delta function
       irho=0       
    ELSE    
       STOP 'WIN_COLDGAS: Error, halo_coldgas not specified correctly'
    END IF

    IF(itype==0 .OR. itype==2 .OR. itype==13) THEN

       ! Density profile of cold gas
       IF(real_space) THEN
          r=k
          win_coldgas=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
          win_coldgas=win_coldgas/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       ELSE
          ! Properly normalise and convert to overdensity
          win_coldgas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
       END IF

       win_coldgas=halo_coldgas_fraction(m,hmod,cosm)*win_coldgas

    ELSE IF(itype==6) THEN

       ! No electron-pressure contribution from the cold gas
       win_coldgas=0.

    ELSE

       STOP 'WIN_COLDGAS: Error, itype not specified correctly'

    END IF

  END FUNCTION win_coldgas

  REAL FUNCTION win_star(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for stars
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: rstar, r, rmin, rmax, p1, p2

    ! Initially set p1, p2
    p1=0.
    p2=0.

    ! Set default maximuim and minimum radii
    rmin=0.
    rmax=rv

    IF(hmod%halo_star==1) THEN
       ! Fedeli (2014)
       irho=7
       rstar=rv/HMx_cstar(m,hmod)
       p1=rstar
       rmax=rv ! Set so that not too much bigger than rstar, otherwise bumps integration goes tits
    ELSE IF(hmod%halo_star==2) THEN
       ! Schneider (2015), following Mohammed (2014)
       irho=9
       !rstar=0.01*rv
       rstar=rv/HMx_cstar(m,hmod)
       p1=rstar
       rmax=10.*rstar ! Set so that not too much bigger than rstar, otherwise bumps integration goes crazy
    ELSE IF(hmod%halo_star==3) THEN
       ! Delta function
       irho=0
    ELSE IF(hmod%halo_star==4) THEN
       ! Transition mass between NFW and delta function
       ! NOTE: mstar here is the same as in the stellar halo-mass fraction. It should probably not be this
       IF(m<hmod%mstar) THEN
          irho=0 ! Delta function
       ELSE
          irho=5 ! NFW
       END IF
    ELSE
       STOP 'WIN_STAR: Error, halo_star specified incorrectly'
    END IF

    IF(real_space) THEN
       r=k
       win_star=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_star=win_star/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE
       ! Properly normalise and convert to overdensity
       win_star=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

    win_star=halo_star_fraction(m,hmod,cosm)*win_star

  END FUNCTION win_star

  REAL FUNCTION win_electron_pressure(real_space,k,m,rv,rs,hmod,cosm)

    !Halo electron pressure profile function for the sum of bound + unbound electron gas
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    INTEGER, PARAMETER :: itype_pressure=6

    IF(hmod%electron_pressure==1) THEN
       ! This overrides everything and just uses the UPP
       win_electron_pressure=UPP(real_space,k,m,rv,rs,hmod,cosm)
    ELSE IF(hmod%electron_pressure==2) THEN
       ! Otherwise use...
       win_electron_pressure=win_boundgas(real_space,itype_pressure,k,m,rv,rs,hmod,cosm)+win_freegas(real_space,itype_pressure,k,m,rv,rs,hmod,cosm)
    ELSE
       STOP 'WIN_ELECTRON_PRESSURE: Error, electron_pressure specified incorrectly'
    END IF

  END FUNCTION win_electron_pressure

  REAL FUNCTION win_void(real_space,k,m,rv,rs,hmod,cosm)

    !Void profile
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2

    ! Initially set p1, p2
    p1=0.
    p2=0.

    ! Set default maximuim and minimum radii
    rmin=0.
    rmax=rv

    IF(hmod%halo_void==1) THEN
       !Top-hat void
       irho=2
       rmin=0.
       rmax=10.*rv
    ELSE
       STOP 'WIN_VOID: Error, halo_void specified incorrectly'
    END IF

    IF(real_space) THEN
       r=k
       win_void=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_void=win_void/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE       
       win_void=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_void

  REAL FUNCTION win_compensated_void(real_space,k,m,rv,rs,hmod,cosm)

    !Profile for compensated voids
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2

    ! Initially set p1, p2
    p1=0.
    p2=0.

    ! Set default maximuim and minimum radii
    rmin=0.
    rmax=rv

    IF(hmod%halo_compensated_void==1) THEN
       !Top-hat
       irho=2
       rmin=0.
       rmax=10.*rv
    ELSE
       STOP 'WIN_COMPENSATED_VOID: Error, halo_compensated_void specified incorrectly'
    END IF  

    IF(real_space) THEN
       r=k
       win_compensated_void=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_compensated_void=win_compensated_void/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE       
       win_compensated_void=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_compensated_void

  REAL FUNCTION win_galaxies(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for all galaxies
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    win_galaxies=win_centrals(real_space,k,m,rv,rs,hmod,cosm)+win_satellites(real_space,k,m,rv,rs,hmod,cosm)
    
  END FUNCTION win_galaxies

  REAL FUNCTION win_centrals(real_space,k,m,rv,rs,hmod,cosm)

    !Halo profile for central galaxies
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2
   
    IF(hmod%has_galaxies .EQV. .FALSE.) CALL init_galaxies(hmod,cosm)

    ! Default minimum and maximum radii
    rmin=0.
    rmax=rv

    ! Default additional halo parameters
    p1=0.
    p2=0.

    ! Delta functions
    irho=0

    IF(real_space) THEN
       r=k
       win_centrals=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_centrals=win_centrals/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE      
       win_centrals=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/hmod%n_c
    END IF

    win_centrals=N_centrals(m,hmod)*win_centrals

  END FUNCTION win_centrals

  REAL FUNCTION win_satellites(real_space,k,m,rv,rs,hmod,cosm)

    ! Halo profile for satellite galaxies
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2
    
    IF(hmod%has_galaxies .EQV. .FALSE.) CALL init_galaxies(hmod,cosm)

    ! Initially set p1, p2
    p1=0.
    p2=0.

    ! Default maximum and minimum radii
    rmin=0.
    rmax=rv

    ! NFW profile
    irho=5

    IF(real_space) THEN
       r=k
       win_satellites=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_satellites=win_satellites/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE
       win_satellites=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/hmod%n_s
    END IF

    win_satellites=N_satellites(m,hmod)*win_satellites
    
  END FUNCTION win_satellites

  REAL FUNCTION N_centrals(m,hmod)

    ! The number of central galaxies as a function of halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m<hmod%mgal) THEN
       N_centrals=0
    ELSE
       N_centrals=1
    END IF
    
  END FUNCTION N_centrals

  REAL FUNCTION N_satellites(m,hmod)

    ! The number of satellite galxies as a function of halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m<hmod%mgal) THEN
       N_satellites=0
    ELSE
       N_satellites=CEILING(m/hmod%mgal)-1
    END IF
    
  END FUNCTION N_satellites

  REAL FUNCTION N_galaxies(m,hmod)

    ! The number of central galaxies as a function of halo mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    N_galaxies=N_centrals(m,hmod)+N_satellites(m,hmod)
    
  END FUNCTION N_galaxies

  REAL FUNCTION win_HI(real_space,k,m,rv,rs,hmod,cosm)

    ! Returns the real or Fourier space HI halo profile
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space ! Real space or Fourier space
    REAL, INTENT(IN) :: k  ! Comoving wave vector (or radius)
    REAL, INTENT(IN) :: m  ! Halo mass
    REAL, INTENT(IN) :: rv ! Halo virial radius
    REAL, INTENT(IN) :: rs ! Halo scale radius 
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halomodel
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    INTEGER :: irho
    REAL :: r, rmin, rmax, p1, p2
    REAL :: r0, alpha, c_HI, r_HI, z

    IF(hmod%has_HI .EQV. .FALSE.) CALL init_HI(hmod,cosm)

    ! Default minimum and maximum radii
    rmin=0.
    rmax=rv

    ! Default additional halo parameters
    p1=0.
    p2=0.

    IF(hmod%halo_HI==1) THEN
       ! NFW profile
       irho=5
    ELSE IF(hmod%halo_HI==2) THEN
       ! Delta function
       irho=0
    ELSE IF(hmod%halo_HI==3) THEN
       ! Polynomial with exponential cut off Villaescusa-Navarro et al. (1804.09180)
       irho=25
       r0=10**(-2.5) ! 0.003 Mpc/h (really small)
       alpha=3.00 ! Tending to homogeneity  
       p1=r0
       p2=alpha
    ELSE IF(hmod%halo_HI==4) THEN
       ! Modified NFW with exponential cut off Villaescusa-Navarro et al. (1804.09180)
       irho=26
       r0=10**(-2.5) ! 0.003 Mpc/h (really small)
       r_HI=10**(-3.0)    
       p1=r0
       p2=r_HI
    ELSE IF(hmod%halo_HI==5) THEN
       ! Modified NFW from Padmanabhan & Refreiger (1607.01021)
       irho=27
       z=hmod%z
       c_HI=130.
       c_HI=4.*c_HI*((M/1e11)**(-0.109))/(1.+z)
       r_HI=rv/c_HI
       p1=r_HI
    ELSE
       STOP 'win_HI: Error, halo_HI not specified correctly'
    END IF

    IF(real_space) THEN
       r=k
       win_HI=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_HI=win_HI/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE
       !win_HI=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho) ! Wrong, but is what I first sent Richard and Kiyo
       win_HI=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/hmod%rho_HI
    END IF

    win_HI=halo_HI_fraction(m,hmod,cosm)*win_HI

  END FUNCTION win_HI

  REAL FUNCTION virial_temperature(M,rv,a,cosm)

    ! Halo physical virial temperature [K]
    ! Calculates the temperature as if pristine gas falls into the halo
    ! Energy is equally distributed between the particles
    IMPLICIT NONE
    REAL :: M ! virial mass
    REAL, INTENT(IN) :: rv ! virial radius
    REAL, INTENT(IN) :: a ! scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm ! cosmology

    REAL, PARAMETER :: modes=3. ! 1/2 k_BT per mode, 3 modes for 3 dimensions 

    virial_temperature=bigG*((M*msun)*mp*cosm%mup)/(a*rv*mpc) ! NEW: a in denominator converts comoving->physical halo radius
    virial_temperature=virial_temperature/(kb*modes/2.) ! Convert to temperature from energy

  END FUNCTION virial_temperature

  REAL FUNCTION UPP(real_space,k,m,rv,rs,hmod,cosm)

    ! Universal electron pressure profile (Arnaud et al. 2010; arxiv:0910.1234)
    ! Note *very* well that this is for *electron* pressure (see Arnaud 2010 and my notes on this)
    ! Note that is is also the physical pressure, and relates to the comoving pressure via (1+z)^3
    ! The units of the pressure profile are [eV/cm^3]
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    TYPE(halomod), INTENT(INOUT) :: hmod 
    TYPE(cosmology), INTENT(INOUT) :: cosm 
    REAL :: r500c, rmin, rmax, a, z, r, m500c, E
    
    REAL, PARAMETER :: alphap=0.12 ! Exponent correction
    REAL, PARAMETER :: b=0.        ! Hydrostatic mass bias
    INTEGER, PARAMETER :: irho=14  ! Set UPP profile

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod,cosm)

    ! Get r500 for UPP
    r500c=exp(find(log(m),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ! [Mpc/h]

    ! Set the radius range for the profile
    rmin=0.
    rmax=rv
    
    IF(real_space) THEN
       r=k
       UPP=rho(r,rmin,rmax,rv,rs,r500c,zero,irho)
    ELSE
       UPP=winint(k,rmin,rmax,rv,rs,r500c,zero,irho,imeth_win)
    END IF

    ! UPP, P(x), equation 4.1 in Ma et al. (2015)
    m500c=exp(find(log(m),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ![Msun/h]
    m500c=m500c*(1.-b) ![Msun/h]

    ! Dimensionless Hubble parameter
    a=hmod%a ! Scale-factor appears explicitly here
    E=sqrt(Hubble2(a,cosm))

    ! Pre-factors from equation 4.1 in Ma et al. (2015) [eV cm^-3 - no h factors]
    UPP=UPP*((m500c/2.1e14)**(alphap+2./3.))*(E**(8./3.))*1.65*(cosm%h/0.7)**2

    ! The standard UPP is written is in physical units
    ! scale to comoving using (1+z)^3 because pressure ~energy density
    z=hmod%z ! Redshift appears explicitly here
    UPP=UPP/(1.+z)**3

  END FUNCTION UPP

  FUNCTION rho(r,rmin,rmax,rv,rs,p1,p2,irho)

    ! This is an UNNORMALISED halo profile of any sort

    ! Types of profile
    ! ================
    !  0 - Delta function at r=0
    !  1 - Isothermal: r^-2
    !  2 - Top hat: r^0
    !  3 - Moore (1999)
    !  4 - NFW (1997)
    !  5 - Analytic NFW
    !  6 - Beta model with beta=2/3
    !  7 - Star profile
    !  8 - Komatsu & Seljak (2001) according to Schneider (2015)
    !  9 - Stellar profile from Schneider (2015)
    ! 10 - Ejected gas profile (Schneider 2015)
    ! 11 - Simplified Komatsu & Seljak (2001) density
    ! 12 - Simplified Komatsu & Seljak (2001) temperature
    ! 13 - Simplified Komatsu & Seljak (2001) pressure
    ! 14 - Universal pressure profile
    ! 15 - Isothermal beta model, beta=0.86 (Ma et al. 2015)
    ! 16 - Isothermal exterior
    ! 17 - Power-law profile
    ! 18 - Cubic profile: r^-3
    ! 19 - Smooth profile (rho = 0, not really physical)
    ! 20 - Exponential profile
    ! 21 - Full Komatsu & Seljak (2001) density
    ! 22 - Full Komatsu & Seljak (2001) temperature
    ! 23 - Full Komatsu & Seljak (2001) pressure
    ! 24 - Cored NFW profile (Copeland, Taylor & Hall 2018)
    ! 25 - Polynomial with central hole (Villaescusa-Navarro et al. 2018)
    ! 26 - Modified NFW with central hole (Villaescusa-Navarro et al. 2018)
    ! 27 - Modified NFW (Padmanabhan & Refregier 2018)
    ! 28 - Shell

    IMPLICIT NONE
    REAL :: rho
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2 ! Standard profile parameters
    INTEGER, INTENT(IN) :: irho
    REAL :: y, ct, t, c, beta, Gamma, r500c, rt, A, re, rstar, B, rb, r0, alpha, rh, eta0
    REAL :: f1, f2
    REAL :: crap

    ! UPP parameters
    REAL, PARAMETER :: P0=6.41
    REAL, PARAMETER :: c500=1.81
    REAL, PARAMETER :: alpha_UPP=1.33
    REAL, PARAMETER :: beta_UPP=4.13
    REAL, PARAMETER :: gamma_UPP=0.31

    ! To stop compile-time warnings
    crap=p2

    IF(r<rmin .OR. r>rmax) THEN
       ! The profile is considered to be zero outside this region
       rho=0.
    ELSE
       IF(irho==0) THEN
          ! Delta function
          ! Do not assign any value to rho as this gets handled properly elsewhere
          ! STOP 'RHO: You should not be here for a delta-function profile'
          rho=0. 
       ELSE IF(irho==1 .OR. irho==16) THEN
          ! Isothermal
          rho=1./r**2
       ELSE IF(irho==2) THEN
          ! Top hat
          rho=1.
       ELSE IF(irho==3) THEN
          ! Moore (1999)
          y=r/rs
          rho=1./((y**1.5)*(1.+y**1.5))
       ELSE IF(irho==4 .OR. irho==5) THEN
          ! NFW (1997)
          y=r/rs
          rho=1./(y*(1.+y)**2)
       ELSE IF(irho==6) THEN
          ! Isothermal beta model (X-ray gas; SZ profiles; beta=2/3 fixed)
          ! Also known as 'cored isothermal profile'
          y=r/rs
          beta=2./3.
          rho=1./((1.+y**2)**(3.*beta/2.))
       ELSE IF(irho==7) THEN
          ! Stellar profile from Fedeli (2014a)
          rstar=p1
          y=r/rstar
          rho=(1./y)*exp(-y)
       ELSE IF(irho==8) THEN
          ! Komatsu & Seljak (2001) profile with NFW transition radius
          ! VERY slow to calculate the W(k) for some reason
          ! Also creates a weird upturn in P(k) that I do not think can be correct
          STOP 'RHO: This is fucked'
          t=sqrt(5.)
          rt=rv/t
          y=r/rs
          c=rs/rv
          ct=c/t
          Gamma=(1.+3.*ct)*log(1.+ct)/((1.+ct)*log(1.+ct)-ct)
          IF(r<=rt) THEN
             ! Komatsu Seljak in the interior
             rho=(log(1.+y)/y)**Gamma
          ELSE
             ! NFW in the outskirts
             A=((rt/rs)*(1.+rt/rs)**2)*(log(1.+rt/rs)/(rt/rs))**Gamma
             rho=A/(y*(1.+y)**2)
          END IF
       ELSE IF(irho==9) THEN
          ! Stellar profile from Schneider (2015) via Mohammed (2014)
          rstar=p1
          rho=exp(-(r/(2.*rstar))**2)/r**2
          ! Converting to y caused the integration to crash for some reason !?!
          !y=r/rs
          !rho=exp(-(y/2.)**2.)/y**2.
       ELSE IF(irho==10) THEN
          ! Ejected gas profile from Schneider (2015)
          re=p1
          rho=exp(-0.5*(r/re)**2)
       ELSE IF(irho==11 .OR. irho==12 .OR. irho==13 .OR. irho==21 .OR. irho==22 .OR. irho==23) THEN
          ! Komatsu & Seljak (2001) profile
          !Gamma=1.18 ! Recommended by Rabold (2017)          
          IF(irho==11 .OR. irho==12 .OR. irho==13) THEN
             Gamma=p1
             y=r/rs
             rho=log(1.+y)/y
          ELSE IF(irho==21 .OR. irho==22 .OR. irho==23) THEN
             c=rv/rs
             Gamma=p1+0.01*(c-6.5)
             !WRITE(*,*) 'Gamma:', Gamma
             eta0=0.00676*(c-6.5)**2+0.206*(c-6.5)+2.48
             !WRITE(*,*) 'eta0:', eta0
             f1=(3./eta0)*(Gamma-1.)/Gamma
             !WRITE(*,*) 'f1:', f1
             f2=c/NFW_factor(c)
             !WRITE(*,*) 'f2:', f2
             B=f1*f2
             IF(B>1.) B=1.
             !WRITE(*,*) 'B:', B
             y=r/rs
             !WRITE(*,*) 'y:', y
             rho=1.-B*(1.-log(1.+y)/y)
             !WRITE(*,*) 'rho:', rho
          ELSE
            STOP 'RHO: Error, irho specified incorrectly'
          END IF
          IF(irho==11 .OR. irho==21) THEN
             ! KS density profile
             rho=rho**(1./(Gamma-1.))
          ELSE IF(irho==12 .OR. irho==22) THEN
             ! KS temperature profile
             rho=rho
          ELSE IF(irho==13 .OR. irho==23) THEN
             ! KS pressure profile
             rho=rho**(Gamma/(Gamma-1.))
          END IF
       ELSE IF(irho==14) THEN
          ! UPP is in terms of r500c, not rv
          r500c=p1
          ! UPP funny-P(x), equation 4.2 in Ma et al. (2015)
          f1=(c500*r/r500c)**gamma_UPP
          f2=(1.+(c500*r/r500c)**alpha_UPP)**((beta_UPP-gamma_UPP)/alpha_UPP)
          rho=P0/(f1*f2)
       ELSE IF(irho==15) THEN
          ! Isothermal beta model
          !beta=0.86 ! from Ma et al. (2015)
          beta=p1
          rho=(1.+(r/rs)**2)**(-3.*beta/2.)
       ELSE IF(irho==16) THEN
          ! Isothermal exterior
          rho=1./r**2
       ELSE IF(irho==17) THEN
          ! Power-law profile
          beta=p1
          rho=r**beta
       ELSE IF(irho==18) THEN
          ! Cubic profile
          rho=r**(-3)
       ELSE IF(irho==19) THEN
          ! Smooth profile
          rho=0.
       ELSE IF(irho==20) THEN
          ! Exponential profile (HI from Padmanabhan et al. 2017)
          re=p1
          rho=exp(-r/re)
       ELSE IF(irho==24) THEN
          ! Cored NFW (Copeland, Taylor & Hall 2018)
          rb=p1
          f1=(r+rb)/rs
          f2=(1.+r/rs)**2
          rho=1./(f1*f2)
       ELSE IF(irho==25) THEN
          ! polynomial with central exponential hole
          r0=p1
          alpha=p2
          rho=(r**(-alpha))*exp(-r0/r)
       ELSE IF(irho==26) THEN
          ! modified NFW with central exponential hole
          r0=p1
          rh=p2
          rho=(1./((0.75+r/rh)*(1.+r/rh)**2))*exp(-r0/r)
       ELSE IF(irho==27) THEN
          ! modified NFW from Padmanabhan & Refregier (2017; 1607.01021)
          rh=p1
          rho=(1./((0.75+r/rh)*(1.+r/rh)**2))
       ELSE IF(irho==28) THEN
          ! Shell
          rho=0.
       ELSE
          STOP 'RHO: Error, irho not specified correctly'
       END IF

    END IF

  END FUNCTION rho

  REAL FUNCTION rhor2at0(irho)

    ! This is the value of rho(r)*r^2 at r=0
    ! For most profiles this is zero, BUT not if rho(r->0) -> r^-2
    ! Note if rho(r->0) -> r^n with n<-2 then the profile mass would diverge!

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: irho

    IF(irho==0) THEN
       STOP 'RHOR2AT0: You should not be here for a delta-function profile'
    ELSE IF(irho==1 .OR. irho==9) THEN
       !1 - Isothermal
       !9 - Stellar profile from Schneider (2015)
       rhor2at0=1.
    ELSE IF(irho==18) THEN
       STOP 'RHOR2AT0: Error, profile diverges at the origin'
    ELSE
       rhor2at0=0.
    END IF

  END FUNCTION rhor2at0

  FUNCTION normalisation(rmin,rmax,rv,rs,p1,p2,irho)

    ! This calculates the normalisation of a halo
    ! This is the integral of 4pir^2*rho(r)*dr between rmin and rmax

    ! Profile results
    !  0 - Delta function (M = 1)
    !  1 - Isothermal (M = 4pi*rv)
    !  2 - Top hat (M = (4pi/3)*rv^3)
    !  3 - Moore (M = (8pi/3)*rv^3*ln(1+c^1.5)/c^3)
    !  4 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
    !  5 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
    !  6 - Beta model with beta=2/3 (M = 4*pi*rs^3*(rv/rs-atan(rv/rs)))
    !  7 - Fedeli stellar model (M = 4*pi*rstar^2 * [1-exp(-rmax/rstar)*(1.+rmax/rstar)]
    !  8 - No
    !  9 - Stellar profile (Schneider (2015)
    ! 10 - Ejected gas profile (Schneider 2015)
    ! 11 - No
    ! 12 - No
    ! 13 - No
    ! 14 - No
    ! 15 - No
    ! 16 - Isothermal shell (M = 4pi*(rmax-rmin))
    ! 17 - No
    ! 18 - Cubic profile
    ! 19 - Smooth profile
    ! 20 - No
    ! 21 - No
    ! 22 - No
    ! 23 - No
    ! 24 - No
    ! 25 - No
    ! 26 - No
    ! 27 - No
    ! 28 - Shell

    IMPLICIT NONE
    REAL :: normalisation
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    REAL :: cmax, re, rstar, beta, rb, c, b

    IF(irho==0) THEN
       ! Delta function
       normalisation=1.
    ELSE IF(irho==1 .OR. irho==16) THEN
       ! Isothermal
       normalisation=4.*pi*(rmax-rmin)
    ELSE IF(irho==2) THEN
       ! Top hat
       normalisation=4.*pi*(rmax**3-rmin**3)/3.
    ELSE IF(irho==3) THEN
       ! Moore et al. (1999)
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=(2./3.)*4.*pi*(rs**3)*log(1.+cmax**1.5)
       END IF
    ELSE IF(irho==4 .OR. irho==5) THEN
       ! NFW (1997)
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=4.*pi*(rs**3)*NFW_factor(cmax)!(log(1.+cmax)-cmax/(1.+cmax))
       END IF
    ELSE IF(irho==6) THEN
       ! Beta model with beta=2/3
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=4.*pi*(rs**3)*(cmax-atan(cmax))
       END IF
    ELSE IF(irho==7) THEN
       ! Fedeli (2014) stellar model
       IF(rmin .NE. 0) THEN
          ! I could actually derive an analytical expression here if this was ever necessary
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          ! This would be even easier if rmax -> infinity (just 4*pi*rstar^2)
          rstar=p1
          normalisation=4.*pi*(rstar**3)*(1.-exp(-rmax/rstar)*(1.+rmax/rstar))
          !normalisation=4.*pi*rstar**3 ! rmax/rstar -> infinity limit (rmax >> rstar)
       END IF
    ELSE IF(irho==9) THEN
       ! Stellar profile from Schneider & Teyssier (2015)       
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          ! Assumed to go on to r -> infinity
          rstar=p1
          normalisation=4.*(pi**(3./2.))*rstar
       END IF
    ELSE IF(irho==10) THEN
       ! Ejected gas profile from Schneider (2015)
       ! Assumed to go on to r -> infinity
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          ! Assumed to go on to r -> infinity
          re=p1
          normalisation=4.*pi*sqrt(pi/2.)*re**3
       END IF
    ELSE IF(irho==17) THEN
       ! Power-law profile
       beta=p1
       normalisation=(4.*pi/(beta+3.))*(rmax**(beta+3.)-rmin**(beta+3.))
    ELSE IF(irho==18) THEN
       ! Cubic profile
       normalisation=4.*pi*log(rmax/rmin)
    ELSE IF(irho==19) THEN
       ! Smooth profile, needs a normalisation I think
       normalisation=1.
    ELSE IF(irho==24) THEN
       ! Cored NFW profile
       rb=p1
       IF(rb==0.) THEN
          ! This is then the standard NFW case
          c=rv/rs
          normalisation=4.*pi*(rs**3)*NFW_factor(c)
       ELSE
          ! Otherwise there is actually a core
          b=rv/rb
          c=rv/rs
          normalisation=(4.*pi*rs**3)/(b-c)**2
          normalisation=normalisation*(b*(b-2.*c)*NFW_factor(c)+(log(1.+b)-b/(1.+c))*c**2)
       END IF
    ELSE IF(irho==28) THEN
       normalisation=4.*pi*rv**3
    ELSE
       ! Otherwise need to do the integral numerically
       ! k=0 gives normalisation
       normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
    END IF

  END FUNCTION normalisation
  
  REAL FUNCTION win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)

    ! This is an UNNORMALISED halo profile of any sort

    ! Types of profile
    ! ================
    !  0 - Delta function at r=0
    !  1 - Isothermal: r^-2
    !  2 - Top hat: constant
    !  3 - No
    !  4 - No
    !  5 - NFW
    !  6 - No
    !  7 - Star profile
    !  8 - No
    !  9 - Stellar profile from Schneider (2015)
    ! 10 - Ejected gas profile (Schneider 2015)
    ! 11 - No
    ! 12 - No
    ! 13 - No
    ! 14 - No
    ! 15 - No
    ! 16 - Isothermal shell
    ! 17 - No
    ! 18 - No
    ! 19 - Smooth profile
    ! 20 - Exponential profile
    ! 21 - No
    ! 22 - No
    ! 23 - No
    ! 24 - Cored NFW profile (Copeland, Taylor & Hall 2018)
    ! 25 - No
    ! 26 - No
    ! 27 - No
    ! 28 - Shell

    ! Calculates the normalised spherical Fourier Transform of the density profile
    ! Note that this means win_norm(k->0)=1
    ! and that win must be between 0 and 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: rmin
    REAL, INTENT(IN) :: rmax
    REAL, INTENT(IN) :: rv
    REAL, INTENT(IN) :: rs
    REAL, INTENT(IN) :: p1
    REAL, INTENT(IN) :: p2
    INTEGER, INTENT(IN) :: irho
    REAL :: re, f1, f2, rstar, kstar, rb    

    IF(k==0.) THEN

       ! If called for the zero mode (e.g. for the normalisation)
       win_norm=1.

    ELSE

       IF(irho==0) THEN
          ! Delta function profile is not localised in Fourier Space
          win_norm=1.
       ELSE IF(irho==1) THEN
          win_norm=wk_isothermal(k*rmax)
       ELSE IF(irho==2) THEN
          ! Analytic for top hat
          win_norm=wk_tophat(k*rmax)
       ELSE IF(irho==5) THEN
          ! Analytic for NFW
          win_norm=win_NFW(k,rmax,rs)
       ELSE IF(irho==7) THEN
          ! Analytic for Fedeli (2014) stellar profile
          rstar=p1
          kstar=k*rstar
          f1=kstar-exp(-rmax/rstar)*(sin(k*rmax)+kstar*cos(k*rmax))
          f2=kstar*(1.+kstar**2)
          win_norm=f1/f2
          !win_norm=1./(1.+kstar**2) !bigRstar -> infinity limit (rmax >> rstar)
       ELSE IF(irho==9) THEN
          ! Only valid if rmin=0 and rmax=inf
          rstar=p1
          win_norm=(sqrt(pi)/2.)*erf(k*rstar)/(k*rstar)
       ELSE IF(irho==10) THEN
          ! Ejected gas profile
          re=p1
          win_norm=exp(-1.5*(k*re)**2.)
       ELSE IF(irho==16) THEN
          ! Isothermal shells
          win_norm=wk_isothermal_2(k*rmax,k*rmin)
       ELSE IF(irho==19) THEN
          ! Smooth profile
          win_norm=0.
       ELSE IF(irho==20) THEN
          ! Exponential profile
          re=p1
          win_norm=1./(1.+(k*re)**2)**2
       ELSE IF(irho==24) THEN
          ! Cored NFW profile
          rb=p1
          IF(rb==0.) THEN
             ! In this case there is no core
             win_norm=win_NFW(k,rmax,rs)
          ELSE
             ! Otherwise there is a core
             win_norm=win_cored_NFW(k,rmax,rs,rb)
          END IF
       ELSE IF(irho==28) THEN
          ! Shell
          win_norm=sinc(k*rv)
       ELSE
          ! Numerical integral over the density profile (slower)
          win_norm=winint(k,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       END IF

    END IF

  END FUNCTION win_norm

  FUNCTION winint(k,rmin,rmax,rv,rs,p1,p2,irho,imeth)

    ! Calculates W(k,M)
    IMPLICIT NONE
    REAL :: winint
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho, imeth

    ! Integration method
    ! imeth = 1 - normal integration
    ! imeth = 2 - bumps with normal integration
    ! imeth = 3 - storage integration
    ! imeth = 4 - bumps with storage integration
    ! imeth = 5 - linear bumps
    ! imeth = 6 - cubic bumps
    ! imeth = 7 - Hybrid with storage and cubic bumps

    ! Bump methods go crazy with some star profiles (those that drop too fast)
    ! You need to make sure that the rmax for the integration does not extend too far out

    ! The hybrid method seems not to be faster for practical calculations here

    IF(imeth==1) THEN
       winint=winint_normal(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,winint_order,acc_win)
    ELSE IF(imeth==2 .OR. imeth==4 .OR. imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
       winint=winint_bumps(k,rmin,rmax,rv,rs,p1,p2,irho,winint_order,acc_win,imeth)
    ELSE IF(imeth==3) THEN
       winint=winint_store(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,winint_order,acc_win)
    ELSE
       STOP 'WININT: Error, imeth not specified correctly'
    END IF

  END FUNCTION winint

  SUBROUTINE winint_diagnostics(rmin,rmax,rv,rs,p1,p2,irho,outfile)

    ! Write out the winint integrand as a function of k
    IMPLICIT NONE
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j
    REAL :: r, k
    REAL, ALLOCATABLE :: integrand(:)

    REAL, PARAMETER :: kmin=1e-1
    REAL, PARAMETER :: kmax=1e2
    INTEGER, PARAMETER :: nr=256 ! Number of points in r
    INTEGER, PARAMETER :: nk=16  ! Number of points in k

    WRITE(*,*) 'WININT_DIAGNOSTICS: Doing these'
    WRITE(*,*) 'WININT_DIAGNOSTICS: minimum r [Mpc/h]:', REAL(rmin)
    WRITE(*,*) 'WININT_DIAGNOSTICS: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'WININT_DIAGNOSTICS: virial radius [Mpc/h]:', REAL(rv)
    WRITE(*,*) 'WININT_DIAGNOSTICS: scale radius [Mpc/h]:', REAL(rs)
    WRITE(*,*) 'WININT_DIAGNOSTICS: concentration:', REAL(rv/rs)
    WRITE(*,*) 'WININT_DIAGNOSTICS: halo parameter 1:', p1
    WRITE(*,*) 'WININT_DIAGNOSTICS: halo parameter 2:', p2
    WRITE(*,*) 'WININT_DIAGNOSTICS: profile number:', irho
    WRITE(*,*) 'WININT_DIAGNOSTICS: outfile: ', TRIM(outfile)

    ALLOCATE(integrand(nk))
    
    OPEN(7,file=outfile)
    DO i=1,nr
       r=progression(0.,rmax,i,nr)
       DO j=1,nk
          k=exp(progression(log(kmin),log(kmax),j,nk))
          integrand(j)=winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)*sinc(r*k)
       END DO
       WRITE(7,*) r/rv, (integrand(j), j=1,nk)
    END DO
    CLOSE(7)

    WRITE(*,*) 'WININT_DIAGNOSTICS: Done'
    WRITE(*,*)
    
  END SUBROUTINE winint_diagnostics

  FUNCTION winint_normal(a,b,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)

    ! Integration routine using 'normal' method to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_normal
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder
    REAL, INTENT(IN) :: acc
    REAL, INTENT(IN) :: a, b
    DOUBLE PRECISION :: sum
    REAL :: r, dr, winold, weight
    INTEGER :: n, i, j

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=2

    winold=0.

    IF(a==b) THEN

       winint_normal=0.

    ELSE

       !Integrates to required accuracy!
       DO j=1,jmax

          !Increase the number of integration points each go until convergence
          n=ninit*(2**(j-1))

          !Set the integration sum variable to zero
          sum=0.

          DO i=1,n

             !Get the weights
             IF(iorder==1) THEN
                !Composite trapezium weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.5
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==2) THEN
                !Composite extended formula weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.416666666666
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.083333333333
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==3) THEN
                !Composite Simpson weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.375
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.166666666666
                ELSE IF(i==3 .OR. i==n-2) THEN
                   weight=0.958333333333
                ELSE
                   weight=1.
                END IF
             ELSE
                STOP 'WININT_NORMAL: Error, order specified incorrectly'
             END IF

             !Now get r and do the function evaluations
             r=progression(a,b,i,n)
             sum=sum+weight*winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)*sinc(r*k)

          END DO

          !The dr are all equally spaced
          dr=(b-a)/REAL(n-1)

          winint_normal=REAL(sum)*dr

          IF((j>jmin) .AND. winint_normal==0.) THEN
             EXIT
          ELSE IF((j>jmin) .AND. (ABS(-1.+winint_normal/winold)<acc)) THEN
             EXIT
          ELSE
             winold=winint_normal
          END IF

       END DO

    END IF

  END FUNCTION winint_normal

  FUNCTION winint_store(a,b,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: winint_store
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder, irho
    REAL, INTENT(IN) :: a, b
    INTEGER :: i, j, n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       winint_store=0.

    ELSE

       !Reset the sum variables for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=winint_integrand(a,rmin,rmax,rv,rs,p1,p2,irho)*sinc(a*k)
             f2=winint_integrand(b,rmin,rmax,rv,rs,p1,p2,irho)*sinc(b*k)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=winint_integrand(x,rmin,rmax,rv,rs,p1,p2,irho)*sinc(x*k)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'WININT_STORE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             winint_store=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             winint_store=0.d0
             STOP 'WININT_STORE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             winint_store=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION winint_store

  FUNCTION winint_bumps(k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc,imeth)

    !Integration routine to calculate the normalised halo FT
    IMPLICIT NONE
    REAL :: winint_bumps
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    INTEGER, INTENT(IN) :: iorder, imeth
    REAL, INTENT(IN) :: acc
    REAL :: sum, w, rn
    REAL :: r1, r2
    REAL :: a3, a2, a1, a0
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n

    INTEGER, PARAMETER :: nlim=3 ! Do the bumps approximation after this number of bumps

    !This MUST be set to zero for this routine
    IF(rmin .NE. 0.) STOP 'WININT_BUMPS: Error, rmin must be zero'

    !Calculate the number of nodes of sinc(k*rmax) for 0<=r<=rmax
    n=FLOOR(k*rmax/pi)

    !Set the sum variable to zero
    sum=0.

    !Integrate over each chunk between nodes separately
    DO i=0,n

       !Set the lower integration limit
       IF(k==0.) THEN
          !Special case when k=0 to avoid division by zero
          r1=0.
       ELSE
          r1=i*pi/k
       END IF

       !Set the upper integration limit
       IF(k==0. .OR. i==n) THEN
          !Special case when on last section because end is rmax, not a node!
          r2=rmax
       ELSE
          r2=(i+1)*pi/k
       END IF

       !WRITE(*,*) i, REAL(r1), REAL(r2)

       !Now do the integration along a section
       IF(imeth==2) THEN
          w=winint_normal(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE IF(k==0. .OR. imeth==4 .OR. (imeth==7 .AND. n<=nlim)) THEN
          w=winint_store(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
       ELSE IF(imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
          IF(i==0 .OR. i==n) THEN
             !First piece done 'normally' because otherwise /0 occurs in cubic
             !Last piece will not generally be over one full oscillation
             w=winint_store(r1,r2,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc)
          ELSE
             IF(imeth==5) THEN
                !Linear approximation to integral between nodes - see notes
                !All results from the analytic integral of a linear polynomial vs. one sine oscillation
                rn=pi*(2*i+1)/(2.*k)
                w=(2./k**2)*winint_integrand(rn,rmin,rmax,rv,rs,p1,p2,irho)*((-1.)**i)/rn !Note there is no sinc here
             ELSE IF(imeth==6 .OR. (imeth==7 .AND. n>nlim)) THEN
                !Cubic approximation to integral between nodes - see notes
                !All results from the analytic integral of a cubic polynomial vs. one sine oscillation
                x1=r1 !Beginning
                x2=r1+1.*(r2-r1)/3. !Middle
                x3=r1+2.*(r2-r1)/3. !Middle
                x4=r2 !End
                y1=winint_integrand(x1,rmin,rmax,rv,rs,p1,p2,irho)/x1 !Note there is no sinc here
                y2=winint_integrand(x2,rmin,rmax,rv,rs,p1,p2,irho)/x2 !Note there is no sinc here
                y3=winint_integrand(x3,rmin,rmax,rv,rs,p1,p2,irho)/x3 !Note there is no sinc here
                y4=winint_integrand(x4,rmin,rmax,rv,rs,p1,p2,irho)/x4 !Note there is no sinc here
                CALL fix_cubic(a3,a2,a1,a0,x1,y1,x2,y2,x3,y3,x4,y4)
                w=-6.*a3*(r2+r1)-4.*a2
                w=w+(k**2)*(a3*(r2**3+r1**3)+a2*(r2**2+r1**2)+a1*(r2+r1)+2.*a0)
                w=w*((-1)**i)/k**4
             ELSE
                STOP 'WININT_BUMPS: Error, imeth specified incorrectly'
             END IF
          END IF
       ELSE
          STOP 'WININT_BUMPS: Error, imeth specified incorrectly'
       END IF

       !WRITE(*,*) i, REAL(r1), REAL(r2), REAL(w)

       sum=sum+w

       !Exit if the contribution to the sum is very tiny
       !This seems to be necessary to prevent crashes
       IF(ABS(w)<acc*ABS(sum)) EXIT

    END DO

    winint_bumps=sum

  END FUNCTION winint_bumps

  FUNCTION winint_integrand(r,rmin,rmax,rv,rs,p1,p2,irho)

    !The integrand for the W(k) integral
    !Note that the sinc function is *not* included
    IMPLICIT NONE
    REAL :: winint_integrand
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho

    IF(r==0.) THEN
       winint_integrand=4.*pi*rhor2at0(irho)
    ELSE
       winint_integrand=4.*pi*(r**2)*rho(r,rmin,rmax,rv,rs,p1,p2,irho)
    END IF

  END FUNCTION winint_integrand

  REAL FUNCTION win_NFW(k,rv,rs)

    ! The analytic normalised (W(k=0)=1) Fourier Transform of the NFW profile
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, rs
    REAL :: c, ks
    REAL :: p1, p2, p3
    REAL :: rmin, rmax

    c=rv/rs
    ks=k*rv/c

    p1=cos(ks)*F_NFW(k,rv,c)
    p2=sin(ks)*G_NFW(k,rv,c)
    p3=sin(ks*c)/(ks*(1.+c))

    rmin=0.
    rmax=rv
    win_NFW=4.*pi*(rs**3)*(p1+p2-p3)/normalisation(rmin,rmax,rv,rs,zero,zero,4)

  END FUNCTION win_NFW

  REAL FUNCTION win_cored_NFW(k,rv,rs,rb)

    ! The analytic normalised (W(k=0)=1) Fourier Transform of the cored NFW profile
    ! Appendix A of Copeland, Taylor & Hall (1712.07112)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, rs, rb
    REAL :: b, c
    REAL :: rmin, rmax

    b=rv/rb
    c=rv/rs

    rmin=0.
    rmax=rv
    
    win_cored_NFW=NFW_factor(c)*win_NFW(k,rv,rs)
    win_cored_NFW=win_cored_NFW+(c/(b-c))*(1./(k*rs))*(G_NFW(k,rv,c)*cos(k*rs)-F_NFW(k,rv,c)*sin(k*rs))
    win_cored_NFW=win_cored_NFW-(c/(b-c))*(1./(k*rs))*(G_NFW(k,rv,b)*cos(k*rb)-F_NFW(k,rv,b)*sin(k*rb))
    win_cored_NFW=4.*pi*(rs**3)*(b/(b-c))*win_cored_NFW/normalisation(rmin,rmax,rv,rs,rb,zero,24)

  END FUNCTION win_cored_NFW

  REAL FUNCTION F_NFW(k,rv,c)

    ! Equation (A3; top) from Copeland, Taylor & Hall (1712.07112)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, c
    REAL :: ks

    ks=k*rv/c
    
    F_NFW=Ci(ks*(1.+c))-Ci(ks)

  END FUNCTION F_NFW

  REAL FUNCTION G_NFW(k,rv,c)

    ! Equation (A3; bottom) from Copeland, Taylor & Hall (1712.07112)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, rv, c
    REAL :: ks

    ks=k*rv/c

    G_NFW=Si(ks*(1.+c))-Si(ks)

  END FUNCTION G_NFW

  REAL FUNCTION b_nu(nu,hmod)

    ! Bias function selection
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       b_nu=b_ps(nu,hmod)
    ELSE IF(hmod%imf==2) THEN
       b_nu=b_st(nu,hmod)
    ELSE IF(hmod%imf==3) THEN
       b_nu=b_Tinker(nu,hmod)
    ELSE IF(hmod%imf==4) THEN
       b_nu=1.
    ELSE
       STOP 'B_NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b_nu

  REAL FUNCTION b_ps(nu,hmod)

    ! Press & Scheter (1974) halo bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc

    dc=hmod%dc

    b_ps=1.+(nu**2-1.)/dc

  END FUNCTION b_ps

  REAL FUNCTION b_st(nu,hmod)

    ! Sheth & Tormen (1999) halo bias (equation 12 in 9901122)
    ! Comes from peak-background split
    ! Haloes defined with SO relative to mean matter density with SC Delta_v relation
    ! A redshift dependent delta_c is used for barrier height, again from SC
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: dc

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707

    dc=hmod%dc
    !dc=1.686

    b_st=1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/dc

  END FUNCTION b_st

  REAL FUNCTION b_Tinker(nu,hmod)

    ! Tinker et al. (2010; 1001.3162) halo bias
    ! TODO: Common functions between this and mass function
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, Dv, dc
    REAL :: alpha, beta, gamma, phi, eta

    ! Parameter arrays from Tinker (2010)
    INTEGER, PARAMETER :: n=9 ! Number of entries in parameter lists
    REAL, PARAMETER :: Delta_v(n)=[200.,300.,400.,600.,800.,1200.,1600.,2400.,3200.]
    REAL, PARAMETER :: alpha0(n)=[0.368,0.363,0.385,0.389,0.393,0.365,0.379,0.355,0.327]
    REAL, PARAMETER :: beta0(n)=[0.589,0.585,0.544,0.543,0.564,0.623,0.637,0.673,0.702]
    REAL, PARAMETER :: gamma0(n)=[0.864,0.922,0.987,1.09,1.20,1.34,1.50,1.68,1.81]
    REAL, PARAMETER :: phi0(n)=[-0.729,-0.789,-0.910,-1.05,-1.20,-1.26,-1.45,-1.50,-1.49]
    REAL, PARAMETER :: eta0(n)=[-0.243,-0.261,-0.261,-0.273,-0.278,-0.301,-0.301,-0.319,-0.336]
    REAL, PARAMETER :: beta_z_exp=0.20
    REAL, PARAMETER :: gamma_z_exp=-0.01
    REAL, PARAMETER :: phi_z_exp=-0.08
    REAL, PARAMETER :: eta_z_exp=0.27

    ! Get these from the halo-model structure
    z=hmod%z
    IF(z>3.) z=3. ! Recommendation from Tinker et al. (2010)
    dc=hmod%dc
    Dv=hmod%Dv

    ! Delta_v dependence
    alpha=find(Dv,Delta_v,alpha0,n,3,3,2)
    beta=find(Dv,Delta_v,beta0,n,3,3,2)
    gamma=find(Dv,Delta_v,gamma0,n,3,3,2)
    phi=find(Dv,Delta_v,phi0,n,3,3,2)
    eta=find(Dv,Delta_v,eta0,n,3,3,2)

    ! Redshift dependence
    beta=beta*(1.+z)**beta_z_exp
    gamma=gamma**(1.+z)**gamma_z_exp
    phi=phi*(1.+z)**phi_z_exp
    eta=eta*(1.+z)**eta_z_exp
  
    b_Tinker=1.+(gamma*nu**2-(1.+2.*eta))/dc+(2.*phi/dc)/(1.+(beta*nu)**(2.*phi))
    
  END FUNCTION b_Tinker

  REAL FUNCTION b2_nu(nu,hmod)

    ! Bias function selection
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       b2_nu=b2_ps(nu,hmod)
    ELSE IF(hmod%imf==2) THEN
       b2_nu=b2_st(nu,hmod)
    ELSE IF(hmod%imf==3) THEN
       STOP 'B2_NU: Error, second-order bias not specified for Tinker mass function'
    ELSE
       STOP 'B2_NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b2_nu

  REAL FUNCTION b2_ps(nu,hmod)

    ! Press & Schechter (1974) second order bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: eps1, eps2, E1, E2, dc

    REAL, PARAMETER :: a2=-17./21.
    REAL, PARAMETER :: p=0.0
    REAL, PARAMETER :: q=1.0

    dc=hmod%dc

    STOP 'B2_PS: Check this very carefully'
    ! I just took the ST form and set p=0 and q=1

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2_ps=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2_ps

  REAL FUNCTION b2_st(nu,hmod)

    ! Sheth, Mo & Tormen (2001) second-order bias
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: eps1, eps2, E1, E2, dc

    ! Notation follows from Cooray & Sheth (2002) pp 25-26

    REAL, PARAMETER :: a2=-17./21.
    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707

    dc=hmod%dc

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2_st=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2_st

  REAL FUNCTION g_nu(nu,hmod)

    ! Mass function
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       g_nu=g_ps(nu,hmod)
    ELSE IF(hmod%imf==2) THEN
       g_nu=g_st(nu,hmod)
    ELSE IF(hmod%imf==3) THEN
       g_nu=g_Tinker(nu,hmod)
    ELSE
       STOP 'G_NU: Error, imf specified incorrectly'
    END IF

  END FUNCTION g_nu

  REAL FUNCTION g_ps(nu,hmod)

    ! Press & Scheter (1974) mass function!
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: crap

    REAL, PARAMETER :: A=sqrt(2./pi)

    ! Stop compile-time warnings
    crap=hmod%a

    g_ps=A*exp(-(nu**2)/2.)

  END FUNCTION g_ps

  REAL FUNCTION g_st(nu,hmod)

    ! Sheth & Tormen (1999) mass function, equation (10) in arXiv:9901122
    ! Note I use nu=dc/sigma(M) and this Sheth & Tormen (1999) use nu=(dc/sigma)^2, which accounts for some small differences
    ! Haloes defined with SO relative to mean matter density with SC Delta_v relation
    ! A redshift dependent delta_c is used for barrier height, again from SC
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: crap

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: A=0.21616

    ! Stop compile-time warnings
    crap=hmod%a

    g_st=A*(1.+((q*nu**2)**(-p)))*exp(-q*nu**2/2.)

  END FUNCTION g_st

  REAL FUNCTION g_Tinker(nu,hmod)

    ! Tinker et al. (2010; 1001.3162) mass function (also 2008; xxxx.xxxx)
    ! TODO: Common functions between this and mass function
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: alpha, beta, gamma, phi, eta
    REAL :: z, Dv

    ! Parameter arrays from Tinker (2010)
    INTEGER, PARAMETER :: n=9 ! Number of entries in parameter lists
    REAL, PARAMETER :: Delta_v(n)=[200.,300.,400.,600.,800.,1200.,1600.,2400.,3200.]
    REAL, PARAMETER :: alpha0(n)=[0.368,0.363,0.385,0.389,0.393,0.365,0.379,0.355,0.327]
    REAL, PARAMETER :: beta0(n)=[0.589,0.585,0.544,0.543,0.564,0.623,0.637,0.673,0.702]
    REAL, PARAMETER :: gamma0(n)=[0.864,0.922,0.987,1.09,1.20,1.34,1.50,1.68,1.81]
    REAL, PARAMETER :: phi0(n)=[-0.729,-0.789,-0.910,-1.05,-1.20,-1.26,-1.45,-1.50,-1.49]
    REAL, PARAMETER :: eta0(n)=[-0.243,-0.261,-0.261,-0.273,-0.278,-0.301,-0.301,-0.319,-0.336]
    REAL, PARAMETER :: beta_z_exp=0.20
    REAL, PARAMETER :: gamma_z_exp=-0.01
    REAL, PARAMETER :: phi_z_exp=-0.08
    REAL, PARAMETER :: eta_z_exp=0.27

    ! Get these from the halo-model structure
    z=hmod%z
    IF(z>3.) z=3. ! Recommendation from Tinker et al. (2010)
    Dv=hmod%Dv

    ! Delta_v dependence
    alpha=find(Dv,Delta_v,alpha0,n,3,3,2)
    beta=find(Dv,Delta_v,beta0,n,3,3,2)
    gamma=find(Dv,Delta_v,gamma0,n,3,3,2)
    phi=find(Dv,Delta_v,phi0,n,3,3,2)
    eta=find(Dv,Delta_v,eta0,n,3,3,2)

    ! Redshift dependence
    beta=beta*(1.+z)**beta_z_exp
    gamma=gamma**(1.+z)**gamma_z_exp
    phi=phi*(1.+z)**phi_z_exp
    eta=eta*(1.+z)**eta_z_exp

    ! The actual mass function
    g_Tinker=alpha*(1.+(beta*nu)**(-2.*phi))*nu**(2.*eta)*exp(-0.5*gamma*nu**2)
    
  END FUNCTION g_Tinker

  REAL FUNCTION gb_nu(nu,hmod)

    ! g(nu) times b(nu)
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    gb_nu=g_nu(nu,hmod)*b_nu(nu,hmod)

  END FUNCTION gb_nu

  FUNCTION wk_isothermal(x)

    ! The normlaised Fourier Transform of an isothermal profile
    IMPLICIT NONE
    REAL :: wk_isothermal
    REAL, INTENT(IN) :: x

    REAL, PARAMETER :: dx=1e-3

    ! Taylor expansion used for low |x| to avoid cancellation problems

    IF(ABS(x)<ABS(dx)) THEN
       ! Taylor series at low x
       wk_isothermal=1.-(x**2)/18.
    ELSE
       wk_isothermal=Si(x)/x
    END IF

  END FUNCTION wk_isothermal

  FUNCTION wk_isothermal_2(x,y)

    ! The normlaised Fourier Transform of an isothemral profile from x -> y
    IMPLICIT NONE
    REAL :: wk_isothermal_2
    REAL, INTENT(IN) :: x, y

    wk_isothermal_2=(Si(x)-Si(y))/(x-y)

  END FUNCTION wk_isothermal_2

  REAL FUNCTION halo_fraction(itype,m,hmod,cosm)

    ! Mass fraction of a type within a halo
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    If(itype==field_dmonly .OR. itype==field_matter) THEN
       halo_fraction=1.
    ELSE IF(itype==field_cdm) THEN
       halo_fraction=halo_CDM_fraction(m,hmod,cosm)
    ELSE IF(itype==field_gas) THEN
       halo_fraction=halo_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_star) THEN
       halo_fraction=halo_star_fraction(m,hmod,cosm)
    ELSE IF(itype==field_bound_gas) THEN
       halo_fraction=halo_boundgas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_free_gas) THEN
       halo_fraction=halo_freegas_fraction(m,hmod,cosm)
    ELSE IF(itype==field_cold_gas) THEN
       halo_fraction=halo_coldgas_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_FRACTION: Error, itype not specified correcntly'
    END IF

  END FUNCTION halo_fraction

  REAL FUNCTION halo_CDM_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in CDM
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap
    
    ! To prevent compile-time warning
    crap=m
    crap=hmod%a

    ! Always the universal value
    halo_CDM_fraction=cosm%om_c/cosm%om_m

  END FUNCTION halo_CDM_fraction

  REAL FUNCTION halo_gas_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in gas
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    halo_gas_fraction=halo_boundgas_fraction(m,hmod,cosm)+halo_freegas_fraction(m,hmod,cosm)+halo_coldgas_fraction(m,hmod,cosm)

  END FUNCTION halo_gas_fraction

  REAL FUNCTION halo_boundgas_fraction(m,hmod,cosm)

    ! Fraction of a halo in bound gas
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: m0, sigma, beta

    IF(hmod%frac_boundgas==1) THEN
       ! From Fedeli (2014a)
       m0=1e12
       sigma=3.
       IF(m<m0) THEN
          halo_boundgas_fraction=0.
       ELSE
          halo_boundgas_fraction=erf(log10(m/m0)/sigma)*cosm%om_b/cosm%om_m
       END IF
    ELSE IF(hmod%frac_boundgas==2) THEN
       ! From Schneider (2015)
       M0=HMx_M0(hmod)
       beta=0.6
       halo_boundgas_fraction=(cosm%om_b/cosm%om_m)/(1.+(M0/m)**beta)
    ELSE IF(hmod%frac_boundgas==3) THEN
       ! Universal baryon fraction model (account for stellar contribution)
       halo_boundgas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_BOUNDGAS_FRACTION: Error, frac_boundgas not specified correctly'
    END IF

    IF(hmod%frac_coldgas==1) THEN
       ! Cold gas is a constant fraction of halo gas
       halo_boundgas_fraction=(1.-hmod%fcold)*halo_boundgas_fraction
    ELSE
       STOP 'HALO_BOUNDGAS_FRACTION: Error, frac_coldgas not specified correctly'
    END IF

  END FUNCTION halo_boundgas_fraction

  REAL FUNCTION halo_coldgas_fraction(m,hmod,cosm)

    ! Fraction of cold gas in haloes
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%frac_coldgas==1) THEN
       ! Constant fraction of halo gas
       halo_coldgas_fraction=hmod%fcold*halo_boundgas_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_COLDGAS_FRACTION: Error, frac_coldgas not specified correctly'
    END IF
    
  END FUNCTION halo_coldgas_fraction

  REAL FUNCTION halo_freegas_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in free gas
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    ! This is necessarily all the gas that is not bound or in stars
    halo_freegas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,hmod,cosm)-halo_boundgas_fraction(m,hmod,cosm)-halo_coldgas_fraction(m,hmod,cosm)
    IF(halo_freegas_fraction<0.) halo_freegas_fraction=0.

  END FUNCTION halo_freegas_fraction

  REAL FUNCTION halo_star_fraction(m,hmod,cosm)

    ! Mass fraction of a halo in stars
    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m0, sigma, A
    REAL :: crap
    
    !REAL, PARAMETER :: fmin=0.01

    crap=cosm%A

    IF(hmod%frac_star==1 .OR. hmod%frac_star==3) THEN
       ! Fedeli (2014)
       A=HMx_Astar(hmod)
       m0=hmod%Mstar
       sigma=hmod%sstar
       halo_star_fraction=A*exp(-((log10(m/m0))**2)/(2.*sigma**2))
       IF(hmod%frac_star==3) THEN
          ! Suggested by Ian, the relation I have is for the central stellar mass
          ! in reality this saturates for high-mass haloes (due to satellite contribution)
          IF(halo_star_fraction<A/3. .AND. m>m0) halo_star_fraction=A/3.
          !IF(halo_star_fraction<fmin .AND. m>m0) halo_star_fraction=fmin
       END IF
    ELSE IF(hmod%frac_star==2) THEN
       ! Constant star fraction
       A=0.005
       halo_star_fraction=A
    ELSE IF(hmod%frac_star==4) THEN
       ! No stars (actually, things crash with exactly zero stars)
       halo_star_fraction=1e-4
    ELSE
       STOP 'HALO_STAR_FRACTION: Error, frac_star specified incorrectly'
    END IF

  END FUNCTION halo_star_fraction

  REAL FUNCTION halo_HI_fraction(M,hmod,cosm)

    ! Dimensionless M_HI/M_halo
    IMPLICIT NONE
    REAL, INTENT(IN) :: M ! Host-halo mass
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: alpha, M0, Mmin
    REAL :: crap

    ! Prevent compile-time warnings
    crap=cosm%A

    IF(hmod%frac_HI==1) THEN
       ! Simple model with hard truncation
       IF(m>=hmod%HImin .AND. m<=hmod%HImax) THEN
          halo_HI_fraction=1.
       ELSE
          halo_HI_fraction=0.
       END IF
    ELSE IF(hmod%frac_HI==2 .OR. hmod%frac_HI==3) THEN
       ! From Villaescusa-Navarro et al. (2018; 1804.09180)
       ! 2 - Just z=0 results
       ! 3 - Using z evolution
       IF(hmod%frac_HI==2 .OR. hmod%z==0.) THEN
          ! FoF values from Table 1
          !alpha=0.24
          !M0=4.3e10
          !Mmin=2e12
          ! FoF-SO values from Table 1
          alpha=0.16
          M0=4.1e10
          Mmin=2.4e12          
       ELSE IF(hmod%z==1.) THEN
          ! FoF values from Table 1
          !alpha=0.53
          !M0=1.5e10
          !Mmin=6e11
          ! FoF-SO values from Table 1
          alpha=0.43
          M0=1.8e10
          Mmin=8.6e11          
       ELSE IF(hmod%z==2.) THEN
          ! FoF values from Table 1
          !alpha=0.60
          !M0=1.3e10
          !Mmin=3.6e11
          ! FoF-SO values from Table 1
          alpha=0.51
          M0=1.5e10
          Mmin=4.6e11        
       ELSE IF(hmod%z==3.) THEN
          ! FoF values from Table 1
          !alpha=0.76
          !M0=2.9e9
          !Mmin=6.7e10
          ! FoF-SO values from Table 1
          alpha=0.69
          M0=3.7e9
          Mmin=9.6e10
       ELSE IF(hmod%z==4.) THEN
          ! FoF values from Table 1
          !alpha=0.79
          !M0=1.4e9
          !Mmin=2.1e10
          ! FoF-SO values from Table 1
          alpha=0.61
          M0=4.5e9
          Mmin=7.6e10         
       ELSE IF(hmod%z==5.) THEN
          ! FoF values from Table 1
          !alpha=0.74
          !M0=1.9e9
          !Mmin=2e10
          ! FoF-SO values from Table 1
          alpha=0.59
          M0=4.1e9
          Mmin=5.4e10          
       ELSE
          STOP 'HALO_HI_FRACTION: Error, redshift not supported here'
       END IF
       halo_HI_fraction=(M0/M)*((M/Mmin)**alpha)*exp(-(M/Mmin)**(-0.35))
    ELSE
       STOP 'HALO_HI_FRACTION: Error, frac_HI specified incorrectly'
    END IF
       
  END FUNCTION halo_HI_fraction

  FUNCTION integrate_hmod(a,b,f,hmod,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_hmod
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(halomod), INTENT(INOUT) :: hmod
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
    
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(nu,hmod)
         IMPORT :: halomod
         REAL, INTENT(IN) :: nu
         TYPE(halomod), INTENT(INOUT) :: hmod
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_hmod=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax
          
          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             ! The first go is just the trapezium of the end points
             f1=f(a,hmod)
             f2=f(b,hmod)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,hmod)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod

  FUNCTION integrate_hmod_cosm(a,b,f,hmod,cosm,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_hmod_cosm
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
    
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(nu,hmod,cosm)
         IMPORT :: halomod
         IMPORT :: cosmology
         REAL, INTENT(IN) :: nu
         TYPE(halomod), INTENT(INOUT) :: hmod
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_hmod_cosm=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax
          
          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             ! The first go is just the trapezium of the end points
             f1=f(a,hmod,cosm)
             f2=f(b,hmod,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,hmod,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod_cosm=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod_cosm

  FUNCTION integrate_hmod_cosm_exp(a,b,f,hmod,cosm,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_hmod_cosm_exp
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
    
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(nu,hmod,cosm)
         IMPORT :: halomod
         IMPORT :: cosmology
         REAL, INTENT(IN) :: nu
         TYPE(halomod), INTENT(INOUT) :: hmod
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_hmod_cosm_exp=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax
          
          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             ! The first go is just the trapezium of the end points
             f1=f(exp(a),hmod,cosm)*exp(a)
             f2=f(exp(b),hmod,cosm)*exp(b)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(exp(x),hmod,cosm)*exp(x)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod_cosm_exp=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod_cosm_exp

  REAL FUNCTION integrate_scatter(c,dc,ih,k,m,rv,hmod,cosm,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL, INTENT(IN) :: c
    REAL, INTENT(IN) :: dc    
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
    
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    REAL, PARAMETER :: nsig=5

    a=c/(1.+nsig*dc)
    b=c*(1.+nsig*dc)
    
    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_scatter=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             ! The first go is just the trapezium of the end points
             f1=scatter_integrand(a,c,dc,ih,k,m,rv,hmod,cosm)
             f2=scatter_integrand(b,c,dc,ih,k,m,rv,hmod,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=scatter_integrand(x,c,dc,ih,k,m,rv,hmod,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_SCATTER: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_SCATTER: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_scatter=REAL(sum_new)

    END IF

  END FUNCTION integrate_scatter

  REAL FUNCTION scatter_integrand(c,mean_c,sigma_lnc,ih,k,m,rv,hmod,cosm)

    ! Integrand for computing halo profiles with scatter
    IMPLICIT NONE
    REAL, INTENT(IN) :: c
    REAL, INTENT(IN) :: mean_c
    REAL, INTENT(IN) :: sigma_lnc
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: m
    REAL, INTENT(IN) :: rv
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: wk(2), pc, rs
    INTEGER :: j

    !Halo profiles
    DO j=1,2
       rs=rv/c
       wk(j)=win_type(.FALSE.,ih(j),k,m,rv,rs,hmod,cosm)
    END DO

    !Probability distribution
    pc=lognormal(c,mean_c,sigma_lnc)

    !The full integrand
    scatter_integrand=wk(1)*wk(2)*pc
    
  END FUNCTION scatter_integrand

  SUBROUTINE winint_speed_tests(k,nk,rmin,rmax,rv,rs,p1,p2,irho)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k(nk), rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: nk, irho
    CHARACTER(len=256) :: base, ext, outfile
    INTEGER :: j, i, ii, imeth, n, ntime
    LOGICAL :: timing
    REAL :: t1, t2, w

    base='winint/results_'
    ext='.dat'

    !j=1 - Time calculation
    !j=2 - Write calculation out
    DO j=1,2

       timing=.FALSE.
       IF(j==2) THEN
          timing=.TRUE.
          WRITE(*,*) 'WININT_SPEED_TESTS: Doing this many evaluations for timing test:', ntime
          WRITE(*,*)
       END IF

       DO imeth=1,7

          WRITE(*,*) 'WININT_SPEED_TESTS: Method:', imeth
          IF(imeth==1) WRITE(*,*) 'WININT_SPEED_TESTS: winint_normal'
          IF(imeth==2) WRITE(*,*) 'WININT_SPEED_TESTS: winint_normal - bumps'
          IF(imeth==3) WRITE(*,*) 'WININT_SPEED_TESTS: winint_store'
          IF(imeth==4) WRITE(*,*) 'WININT_SPEED_TESTS: winint_store - bumps'
          IF(imeth==5) WRITE(*,*) 'WININT_SPEED_TESTS: linear - bumps'
          IF(imeth==6) WRITE(*,*) 'WININT_SPEED_TESTS: cubic - bumps'
          IF(imeth==7) WRITE(*,*) 'WININT_SPEED_TESTS: hybrid'

          IF(timing .EQV. .FALSE.) THEN
             outfile=number_file(base,imeth,ext)
             WRITE(*,*) 'WININT_SPEED_TESTS: Writing data: ', TRIM(outfile)
             OPEN(7,file=outfile)            
             n=1
             IF(imeth==1) CALL cpu_time(t1)
          ELSE
             CALL cpu_time(t1)
             n=ntime
          END IF

          ! Loop over number of iterations
          DO ii=1,n

             ! Loop over wave number and do integration
             DO i=1,nk
                w=winint(k(i),rmin,rmax,rv,rs,p1,p2,irho,imeth)/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
                IF(.NOT. timing) THEN
                   WRITE(7,*) k(i), w
                END IF
             END DO

          END DO

          IF(.NOT. timing) THEN
             CLOSE(7)
             IF(imeth==1) THEN
                CALL cpu_time(t2)
                ntime=CEILING(winint_test_seconds/(t2-t1))
             END IF
          ELSE
             CALL cpu_time(t2)
             WRITE(*,*) 'WININT_SPEED_TESTS: Time [s]:', t2-t1
          END IF

          WRITE(*,*)

       END DO

    END DO

  END SUBROUTINE winint_speed_tests

END MODULE HMx
