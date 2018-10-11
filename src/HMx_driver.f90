PROGRAM HMx_driver

  USE HMx
  USE Limber
  USE file_info
  USE random_numbers
  USE cosmic_emu_stuff
 
  IMPLICIT NONE

  ! Parameter definitions
  REAL, ALLOCATABLE :: k(:), a(:)
  REAL, ALLOCATABLE :: pow(:), pow_lin(:), pow_2h(:), pow_1h(:), pow_full(:)
  REAL, ALLOCATABLE :: powa(:,:), powa_lin(:,:), powa_2h(:,:), powa_1h(:,:), powa_full(:,:), pows_full(:,:,:)
  REAL, ALLOCATABLE :: ell(:), Cell(:), theta(:), xi(:,:), zs(:)
  REAL, ALLOCATABLE :: z_tab(:), HI_frac(:)
  INTEGER :: i, j, ii, l, ik, nk, na, j1, j2, iz, itriad
  INTEGER :: n, nl, nz, nth, nnz, m, ipa, npa, ncos, ncore, ntriad
  INTEGER :: ip(2), ix(2), ixx(2), ihalo
  REAL :: kmin, kmax, amin, amax, lmin, lmax, thmin, thmax, zmin, zmax
  REAL :: rcore_min, rcore_max
  REAL :: z, z1, z2, r1, r2, a1, a2
  TYPE(cosmology) :: cosm
  TYPE(cosmology), ALLOCATABLE :: cosms(:)
  TYPE(halomod) :: hmod
  TYPE(halomod), ALLOCATABLE :: hmods(:)
  TYPE(projection) :: proj(2)
  TYPE(lensing) :: lens
  CHARACTER(len=256) :: infile, outfile, base, mid, ext, dir, name, fname, inbase, outbase, inext, outext
  CHARACTER(len=256) :: mode, halomodel, red, cosmo
  INTEGER :: imode, icosmo, iowl, ihm, irho, itest, jtest
  REAL :: sig8min, sig8max
  REAL :: mass, m1, m2, nu, nu_min, nu_max, mf
  REAL :: c, rmin, rmax, rv, rs, p1, p2
  REAL :: spam
  CHARACTER(len=1) :: crap
  
  ! Baryon stuff
  REAL :: param_min, param_max, param, param_neat
  LOGICAL :: ilog

  ! Halo-model Parameters
  LOGICAL, PARAMETER :: verbose=.TRUE. ! Verbosity
  REAL, PARAMETER :: mmin=1e7 ! Minimum halo mass for the calculation
  REAL, PARAMETER :: mmax=1e17 ! Maximum halo mass for the calculation

  ! Test parameters
  REAL, PARAMETER :: tolerance=3e-3
  REAL :: error, error_max
  LOGICAL, PARAMETER :: verbose_tests=.FALSE.
  LOGICAl :: ifail=.FALSE.

  ! Benchmark parameters
  LOGICAL, PARAMETER :: Alonso_k=.TRUE.

  ! Output choices
  LOGICAL, PARAMETER :: icumulative=.TRUE. ! Do cumlative distributions for breakdown
  LOGICAL, PARAMETER :: ixi=.TRUE. ! Do correlation functions from C(l)
  LOGICAL, PARAMETER :: ifull=.FALSE. ! Do only full halo model C(l), xi(theta) calculations (quicker, no breakdown ...)

  ! Fitting
  REAL, ALLOCATABLE :: pow_sim(:), pows_sim(:,:,:), ks_sim(:,:,:), k_sim(:)
  REAL, ALLOCATABLE :: pmin(:), pmax(:), pbest(:), pnow(:), prange(:)
  REAL, ALLOCATABLE :: pnew(:), pold(:), porig(:)
  LOGICAL, ALLOCATABLE :: plog(:)
  INTEGER, ALLOCATABLE :: nrange(:), iii(:), ipbest(:)
  REAL :: fom, fom_best, fom_new, fom_old, fom_orig
  INTEGER :: i1, i2, i3, ifit, nref, np, ibest!, ipbest(3), iii(3)
  LOGICAL :: refine, jump
  CHARACTER(len=256) :: model

  CALL get_command_argument(1,mode)
  IF(mode=='') THEN
     imode=-1
  ELSE
     READ(mode,*) imode
  END IF

  CALL get_command_argument(2,cosmo)
  IF(cosmo=='') THEN
     icosmo=-1
  ELSE
     READ(cosmo,*) icosmo
  END IF

  CALL get_command_argument(3,halomodel)
  IF(halomodel=='') THEN
     ihm=-1
  ELSE
     READ(halomodel,*) ihm
  END IF

  ! Initial white space
  WRITE(*,*)

  ! Choose mode
  IF(imode==-1) THEN
     WRITE(*,*) 'HMx_DRIVER: Choose what to do'
     WRITE(*,*) '============================='
     WRITE(*,*) ' 0 - Gravity-only power spectrum at z=0'
     WRITE(*,*) ' 1 - 3D Matter power spectrum over multiple z'
     WRITE(*,*) ' 2 - Produce all halo components 3D cross and auto spectra'
     WRITE(*,*) ' 3 - Run diagnostics for haloes'
     WRITE(*,*) ' 4 - Do random cosmologies for bug testing'
     WRITE(*,*) ' 5 - NOT SUPPORTED: Pressure field comparison'
     WRITE(*,*) ' 6 - n(z) check'
     WRITE(*,*) ' 7 - Do general angular cross correlation'
     WRITE(*,*) ' 8 - Angular cross correlation as a function of cosmology'
     WRITE(*,*) ' 9 - Breakdown angular correlations in halo mass'
     WRITE(*,*) '10 - Breakdown angular correlations in redshift'
     WRITE(*,*) '11 - Breakdown angular correlations in halo radius'
     WRITE(*,*) '12 - Project triad'
     WRITE(*,*) '13 - Cross-correlation coefficient'
     WRITE(*,*) '14 - 3D spectra for variations in baryon parameters'
     WRITE(*,*) '15 - 3D spectra for cosmo-OWLS models (k values)'
     WRITE(*,*) '16 - 3D spectra for BAHAMAS models (k values)'
     WRITE(*,*) '17 - 3D spectra for user choice of fields'
     WRITE(*,*) '18 - 3D bias'
     WRITE(*,*) '19 - Make CCL benchmark data'
     WRITE(*,*) '20 - Make data for Ma et al. (2015) Fig. 1'
     WRITE(*,*) '21 - W(k) integrand diagnostics'
     WRITE(*,*) '22 - Time W(k) integration methods'
     WRITE(*,*) '23 - Produce DE response results from Mead (2017)'
     WRITE(*,*) '24 - Grid fitting'
     WRITE(*,*) '25 - MCMC fitting for Mira Titan'
     WRITE(*,*) '26 - TESTS: DMONLY spectra; HMcode'
     WRITE(*,*) '27 - Comparison with Mira Titan nodes'
     WRITE(*,*) '28 - Comparison with FrankenEmu nodes'
     WRITE(*,*) '29 - Comparison with random Mira Titan cosmology'
     WRITE(*,*) '30 - Comparison with random FrankenEmu cosmology'
     WRITE(*,*) '31 - PAPER: Breakdown 3D hydro power in halo mass'
     WRITE(*,*) '32 - PAPER: Hydro power for baseline model'
     WRITE(*,*) '33 - PAPER: Effect of parameter variations on baseline model hydro power'
     WRITE(*,*) '34 - Write HMx hydro parameter variations with T_AGN and z'
     WRITE(*,*) '35 - Power spectra of cored halo profiles'
     WRITE(*,*) '36 - TESTS: hydro spectra'
     WRITE(*,*) '37 - Produce CFHTLenS correlation functions'
     WRITE(*,*) '38 - MCMC fitting for FrankenEmu nodes'
     WRITE(*,*) '39 - MCMC fitting for hydro'
     WRITE(*,*) '40 - MCMC fitting for random Mira Titan'
     WRITE(*,*) '41 - MCMC fitting for random FrankenEmu'
     WRITE(*,*) '42 - PAPER: Contributions to k-k C(l) integral'
     WRITE(*,*) '43 - PAPER: Contributions to k-y C(l) integral'
     WRITE(*,*) '44 - PAPER: Project triad'
     WRITE(*,*) '45 - Comparison of Sheth-Tormen vs. Tinker mass function'
     WRITE(*,*) '46 - Mass function and bias plots'
     WRITE(*,*) '47 - Make CMB lensing to compare with CAMB'
     WRITE(*,*) '48 - HI bias'
     WRITE(*,*) '49 - HI mass fractions'
     READ(*,*) imode
     WRITE(*,*) '============================'
     WRITE(*,*)
  END IF

  IF(imode==0) THEN

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Sets the redshift
     !WRITE(*,*) 'HMx_DRIVER: Choose redshift:'
     !READ(*,*) z
     !WRITE(*,*)
     z=0.

     !Initiliasation for the halomodel calcualtion
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     ! Do the halo-model calculation
     ip(1)=-1
     ip(2)=-1
     CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)

     ! Write out the results
     outfile='data/power.dat'
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     ! Write the one-void term if necessary
     IF(hmod%voids) THEN
        OPEN(8,file='data/power_1void.dat')
        DO i=1,nk     
           WRITE(8,*) k(i), p_1void(k(i),hmod)
        END DO
        CLOSE(8)
     END IF

  ELSE IF(imode==1) THEN

     ! Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     CALL assign_halomod(ihm,hmod,verbose)
     
     ! Set number of k points and k range (log spaced)
     ! The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Set the number of redshifts and range (linearly spaced) and convert z -> a
     nz=16
     zmin=0.
     zmax=4.
     CALL fill_array(zmin,zmax,a,nz)
     a=1./(1.+a)
     na=nz

     ip=-1 ! Set DMONLY profiles 
     CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response=.FALSE.)

     base='data/power'
     CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

  ELSE IF(imode==2 .OR. imode==15 .OR. imode==16 .OR. imode==32) THEN

     ! Make cross/auto power spectra of all different components of haloes as well as pressure

     ! Generic hydro
     IF(imode==2 .OR. imode==32) THEN     
        
        ! Only do one 'model' here
        n=1

        ! Set the redshift
        nz=4
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.0
        z_tab(2)=0.5
        z_tab(3)=1.0
        z_tab(4)=2.0
     
        ! Set number of k points and k range (log spaced)
        nk=128
        kmin=1e-3
        kmax=1e2
        CALL fill_array(log(kmin),log(kmax),k,nk)
        k=exp(k)

     ! cosmo-OWLS
     ELSE IF(imode==15) THEN

        ! Do from REF, NOCOOL, AGN, AGN 8.5, AGN 8.7
        n=5

        ! Set the redshift
        nz=1
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.

        ! Get the k values from the simulation measured P(k)
        infile='/Users/Mead/Physics/cosmo-OWLS/power/N800/DMONLY_all_all_power.dat'
        CALL get_k_values(infile,k,nk)

     ! BAHAMAS
     ELSE IF(imode==16) THEN
       
        ! Do AGN, AGN-lo and AGN-hi
        n=3

        ! Set the redshift
        nz=4
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.0
        z_tab(2)=0.5
        z_tab(3)=1.0
        z_tab(4)=2.0

        ! Get the k values from the simulation measured P(k)
        infile='/Users/Mead/Physics/BAHAMAS/power/M1024/DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'
        CALL get_k_values(infile,k,nk)
        
     END IF

     !Allocate the arrays for P(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     DO iowl=1,n

        DO j=1,nz

           z=z_tab(j)

           !Assigns the cosmological model
           IF(imode==2)  THEN
              icosmo=4
           ELSE IF(imode==15) THEN
              icosmo=2
           ELSE IF(imode==16) THEN
              icosmo=4
           ELSE IF(imode==32) THEN
              icosmo=4
              ihm=3
           END IF
           
           CALL assign_cosmology(icosmo,cosm,verbose)
           CALL init_cosmology(cosm)
           CALL print_cosmology(cosm)
           CALL assign_halomod(ihm,hmod,verbose)
           
           IF(imode==15) THEN
           
              ! cosmo-OWLS
              IF(imode==15 .AND. iowl==1) THEN
                 name='REF'
                 fname=name
                 ! From my fitting by eye
                 hmod%alpha=2.
                 hmod%eps=1.
                 hmod%Gamma=1.24
                 hmod%M0=1e13
                 hmod%Astar=0.055
              ELSE IF(imode==15 .AND. iowl==2) THEN
                 name='NOCOOL'
                 fname=name
                 ! From my fitting by eye
                 hmod%alpha=2.
                 hmod%eps=1.
                 hmod%Gamma=1.1
                 hmod%M0=0.
                 hmod%Astar=0.
              ELSE IF(imode==15 .AND. iowl==3) THEN
                 name='AGN'
                 fname=name
                 ! From Tilman's preliminary results
                 hmod%alpha=0.52
                 hmod%eps=1.
                 hmod%Gamma=1.17
                 hmod%M0=1.047e14
                 hmod%Astar=0.02
              ELSE IF(imode==15 .AND. iowl==4) THEN
                 name='AGN 8.5'
                 fname='AGN8p5'
                 ! From Tilman's preliminary results
                 hmod%alpha=0.56
                 hmod%eps=1.
                 hmod%Gamma=1.19
                 hmod%M0=3.548e14
                 hmod%Astar=0.01
              ELSE IF(imode==15 .AND. iowl==5) THEN
                 name='AGN 8.7'
                 fname='AGN8p7'
                 ! From Tilman's preliminary results
                 hmod%alpha=0.53
                 hmod%eps=1.
                 hmod%Gamma=1.21
                 hmod%M0=7.586e14
                 hmod%Astar=0.01
              END IF

           END IF

           !BAHAMAS
           IF(imode==16) THEN

              IF(iowl==1) THEN

                 ! Simulation name and file name
                 name='AGN'
                 fname='AGN'
  
                 ! Best z=0 fit on 21/06/2018
                 IF(ihm==4) THEN
                    
                    IF(z==0.) THEN
                       hmod%alpha=0.379
                       hmod%eps=10**(-0.061)
                       hmod%Gamma=1.205
                       hmod%M0=10**(13.823)
                       hmod%Astar=0.029
                       hmod%Twhim=10**(5.754)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.537
                       hmod%eps=10**(-0.209)
                       hmod%Gamma=1.163
                       hmod%M0=10**(13.964)
                       hmod%Astar=0.024
                       hmod%Twhim=10**(5.750)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.755
                       hmod%eps=10**(-0.985)
                       hmod%Gamma=1.162
                       hmod%M0=10**(13.673)
                       hmod%Astar=0.019
                       hmod%Twhim=10**(5.057)
                    END IF
                    
                 ELSE IF(ihm==6) THEN
                    
                    IF(z==0.) THEN
                       hmod%alpha=0.195
                       hmod%eps=10**(-0.484)
                       hmod%Gamma=1.399
                       hmod%M0=10**(13.807)
                       hmod%Astar=0.020
                       hmod%Twhim=10**(6.049)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.619
                       hmod%eps=10**(-0.309)
                       hmod%Gamma=1.507
                       hmod%M0=10**(14.937)
                       hmod%Astar=0.021
                       hmod%Twhim=10**(5.987)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.384
                       hmod%eps=10**(-0.475)
                       hmod%Gamma=1.183
                       hmod%M0=10**(14.562)
                       hmod%Astar=0.017
                       hmod%Twhim=10**(5.796)
                    END IF
                    
                 ELSE IF(ihm==3 .OR. ihm==14) THEN
                    
                    IF(z==0.) THEN
                       hmod%alpha=0.428
                       hmod%eps=10**(0.015)
                       hmod%Gamma=1.287
                       hmod%M0=10**(13.233)
                       hmod%Astar=0.030
                       hmod%Twhim=10**(5.404)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.742
                       hmod%eps=10**(0.148)
                       hmod%Gamma=1.516
                       hmod%M0=10**(12.688)
                       hmod%Astar=0.026
                       hmod%Twhim=10**(5.531)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.830
                       hmod%eps=10**(0.169)
                       hmod%Gamma=1.487
                       hmod%M0=10**(12.004)
                       hmod%Astar=0.024
                       hmod%Twhim=10**(5.643)
                    END IF

                 ELSE IF(ihm==17 .OR. ihm==18 .OR. ihm==19) THEN

                    hmod%Theat=10**7.8
                    
                 END IF
                 
              ELSE IF(iowl==3) THEN

                 ! Simulation name and file name
                 name='AGN high'
                 fname='AGN-hi'
                 hmod%Theat=10**8.0

                 ! Best z=0 fit on 21/06/2018
                 IF(ihm==4) THEN

                    IF(z==0.) THEN
                       hmod%alpha=0.421
                       hmod%eps=10**(-0.154)
                       hmod%Gamma=1.211
                       hmod%M0=10**(14.316)
                       hmod%Astar=0.026
                       hmod%Twhim=10**(5.801)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.599
                       hmod%eps=10**(-0.686)
                       hmod%Gamma=1.158
                       hmod%M0=10**(14.455)
                       hmod%Astar=0.022
                       hmod%Twhim=10**(5.849)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.659
                       hmod%eps=10**(-1.993)
                       hmod%Gamma=1.151
                       hmod%M0=10**(14.386)
                       hmod%Astar=0.018
                       hmod%Twhim=10**(5.829)
                    END IF

                 ELSE IF(ihm==6) THEN

                    IF(z==0.) THEN
                       hmod%alpha=0.302
                       hmod%eps=10**(-0.429)
                       hmod%Gamma=1.339
                       hmod%M0=10**(14.662)
                       hmod%Astar=0.026
                       hmod%Twhim=10**(6.080)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.751
                       hmod%eps=10**(-0.013)
                       hmod%Gamma=1.502
                       hmod%M0=10**(14.958)
                       hmod%Astar=0.014
                       hmod%Twhim=10**(5.959)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.371
                       hmod%eps=10**(-0.127)
                       hmod%Gamma=1.101
                       hmod%M0=10**(14.966)
                       hmod%Astar=0.018
                       hmod%Twhim=10**(6.040)
                    END IF

                 ELSE IF(ihm==3 .OR. ihm==14) THEN

                    IF(z==0.) THEN
                       hmod%alpha=0.528
                       hmod%eps=10**(0.038)
                       hmod%Gamma=1.505
                       hmod%M0=10**(13.638)
                       hmod%Astar=0.027
                       hmod%Twhim=10**(5.078)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.742
                       hmod%eps=10**(0.125)
                       hmod%Gamma=1.547
                       hmod%M0=10**(13.481)
                       hmod%Astar=0.024
                       hmod%Twhim=10**(5.786)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.918
                       hmod%eps=10**(0.289)
                       hmod%Gamma=1.996
                       hmod%M0=10**(13.022)
                       hmod%Astar=0.022
                       hmod%Twhim=10**(5.849)
                    END IF

                 ELSE IF(ihm==17 .OR. ihm==18 .OR. ihm==19) THEN

                    hmod%Theat=10**8.0

                 END IF
                 
              ELSE IF(iowl==2) THEN

                 ! Simulation name and file name
                 name='AGN low'
                 fname='AGN-lo'
                 hmod%Theat=10**7.6
                 
                 ! Best z=0 21/06/2018
                 IF(ihm==4) THEN

                    IF(z==0.) THEN
                       hmod%alpha=0.353
                       hmod%eps=10**(-0.017)
                       hmod%Gamma=1.199
                       hmod%M0=10**(13.517)
                       hmod%Astar=0.031
                       hmod%Twhim=10**(5.767)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.491
                       hmod%eps=10**(-0.104)
                       hmod%Gamma=1.158
                       hmod%M0=10**(13.663)
                       hmod%Astar=0.025
                       hmod%Twhim=10**(5.721)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.532
                       hmod%eps=10**(-0.990)
                       hmod%Gamma=1.079
                       hmod%M0=10**(13.815)
                       hmod%Astar=0.021
                       hmod%Twhim=10**(5.601)
                    END IF

                 ELSE IF(ihm==6) THEN

                    IF(z==0.) THEN
                       hmod%alpha=0.181
                       hmod%eps=10**(-0.489)
                       hmod%Gamma=1.432
                       hmod%M0=10**(13.632)
                       hmod%Astar=0.023
                       hmod%Twhim=10**(6.099)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.420
                       hmod%eps=10**(-0.413)
                       hmod%Gamma=1.467
                       hmod%M0=10**(14.606)
                       hmod%Astar=0.021
                       hmod%Twhim=10**(5.874)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.715
                       hmod%eps=10**(-0.405)
                       hmod%Gamma=1.154
                       hmod%M0=10**(14.786)
                       hmod%Astar=0.019
                       hmod%Twhim=10**(5.434)
                    END IF

                 ELSE IF(ihm==3 .OR. ihm==14) THEN

                    IF(z==0.) THEN
                       hmod%alpha=0.409
                       hmod%eps=10**(0.045)
                       hmod%Gamma=1.275
                       hmod%M0=10**(12.737)
                       hmod%Astar=0.032
                       hmod%Twhim=10**(5.357)
                    ELSE IF(z==0.5) THEN
                       hmod%alpha=0.698
                       hmod%eps=10**(0.159)
                       hmod%Gamma=1.393
                       hmod%M0=10**(12.012)
                       hmod%Astar=0.028
                       hmod%Twhim=10**(5.491)
                    ELSE IF(z==1. .OR. z==2.) THEN
                       hmod%alpha=0.832
                       hmod%eps=10**(0.312)
                       hmod%Gamma=1.344
                       hmod%M0=10**(12.020)
                       hmod%Astar=0.025
                       hmod%Twhim=10**(5.723)
                    END IF

                 ELSE IF(ihm==17 .OR. ihm==18 .OR. ihm==19) THEN

                    hmod%Theat=10**7.6

                 END IF

              END IF

           END IF

           IF(imode==15) WRITE(*,*) 'Comparing to OWLS model: ', TRIM(name)
           IF(imode==16) WRITE(*,*) 'Comparing to BAHAMAS model: ', TRIM(name)
            
           !Initiliasation for the halomodel calcualtion after variables changed
           CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
           CALL print_halomod(hmod,cosm,verbose)

           !Runs the diagnostics
           IF(imode==2) THEN
              dir='data'
              CALL halo_diagnostics(hmod,cosm,dir)
              CALL halo_definitions(hmod,cosm,dir)
              CALL halo_properties(hmod,dir)
           END IF

           IF(imode==2 .OR. imode==32) THEN
              !File base and extension
              IF(j==1) base='data/power_z0.0_'
              IF(j==2) base='data/power_z0.5_'
              IF(j==3) base='data/power_z1.0_'
              IF(j==4) base='data/power_z2.0_'
              mid=''
              ext='.dat'
           ELSE IF(imode==15) THEN
              base='cosmo-OWLS/power_'//TRIM(fname)//'_'
              mid=''
              ext='.dat'
           ELSE IF(imode==16) THEN
              IF(j==1) base='BAHAMAS/power_'//TRIM(fname)//'_z0.0_'
              IF(j==2) base='BAHAMAS/power_'//TRIM(fname)//'_z0.5_'
              IF(j==3) base='BAHAMAS/power_'//TRIM(fname)//'_z1.0_'
              IF(j==4) base='BAHAMAS/power_'//TRIM(fname)//'_z2.0_'
              mid=''
              ext='.dat'
           END IF

           !Dark-matter only
           IF(imode==2 .OR. imode==32) THEN
              IF(j==1) outfile='data/power_z0.0.dat'
              IF(j==2) outfile='data/power_z0.5.dat'
              IF(j==3) outfile='data/power_z1.0.dat'
              IF(j==4) outfile='data/power_z2.0.dat'
           ELSE IF(imode==15) THEN
              outfile='cosmo-OWLS/power_DMONLY_00.dat'
           ELSE IF(imode==16) THEN
              IF(j==1) outfile='BAHAMAS/power_DMONLY_z0.0_00.dat'
              IF(j==2) outfile='BAHAMAS/power_DMONLY_z0.5_00.dat'
              IF(j==3) outfile='BAHAMAS/power_DMONLY_z1.0_00.dat'
              IF(j==4) outfile='BAHAMAS/power_DMONLY_z2.0_00.dat'
           END IF
           ip=-1
           WRITE(*,*) ip(1), ip(2), TRIM(outfile)
           CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

           ! Loop over matter types and do auto- and cross-spectra
           DO j1=0,6
              DO j2=j1,6

                 ! Skip for the bound- and free-gas spectra, fuck 'em
                 IF(j1==4 .OR. j1==5) CYCLE
                 IF(j2==4 .OR. j2==5) CYCLE

                 ! Fix output file and write to screen
                 outfile=number_file2(base,j1,mid,j2,ext)

                 ! Set the halo types and write to screen
                 ip(1)=j1
                 ip(2)=j2
                 WRITE(*,*) ip(1), ip(2), TRIM(outfile)

                 ! Do the calculation and write P(k) to disk
                 CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose=.FALSE.,response=.FALSE.)
                 CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose=.FALSE.)

              END DO
           END DO

        END DO

     END DO

  ELSE IF(imode==3) THEN

     !Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Loop over redshifts
     DO j=1,4

        !Set the redshift
        IF(j==1) z=0.0
        IF(j==2) z=0.5
        IF(j==3) z=1.0
        IF(j==4) z=2.0

        !Initiliasation for the halomodel calcualtion
        CALL assign_halomod(ihm,hmod,verbose)
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        !Runs the diagnostics
        dir='data'
        CALL halo_diagnostics(hmod,cosm,dir)
        CALL halo_definitions(hmod,cosm,dir)
        CALL halo_properties(hmod,dir)

     END DO

  ELSE IF(imode==4) THEN

     ! Random baryon parameters

     ! Ignore this, only useful for bug tests
     CALL RNG_set(0)

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Assign the default halo model
     CALL assign_halomod(ihm,hmod,verbose)

     ! Set the redshift
     z=0.

     DO

        ! Initiliasation for the halomodel calcualtion        
        CALL random_baryon_parameters(hmod)
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Do the halo-model calculation for a range of halo types
        DO j1=0,6
           DO j2=j1,6
              IF(j1==1 .OR. j1==4 .OR. j1==5) CYCLE
              IF(j2==1 .OR. j2==4 .OR. j2==5) CYCLE
              ip(1)=j1
              ip(2)=j2
              WRITE(*,*) 'Halo types:', ip(1), ip(2)
              CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose=.FALSE.,response=.FALSE.)
              DO i=1,nk
                 IF(ISNAN(pow_full(i))) STOP 'HMX_DRIVER: Error, NaN found in pow_full array'
              END DO
           END DO
        END DO

        WRITE(*,*)

     END DO

  ELSE IF(imode==5) THEN

     STOP 'HMX_DRIVER: Error, imode=5 not supported any more'

     !Create spectra for 'all matter' and 'electron pressure' and their cross

     !Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Assigns the cosmological model
     icosmo=2
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     !File base and extension
     base='data/power_'
     ext='.dat'

     DO i=1,4

        !Set the redshift
        IF(i==1) THEN
           z=0.0
           red='z0.0'
        ELSE IF(i==2) THEN
           z=0.5
           red='z0.5'
        ELSE IF(i==3) THEN
           z=1.0
           red='z1.0'
        ELSE IF(i==4) THEN
           z=2.0
           red='z2.0'
        END IF

        !Initiliasation for the halomodel calcualtion
        CALL assign_halomod(ihm,hmod,verbose)
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        !Runs the diagnostics
        dir='data'
        CALL halo_diagnostics(hmod,cosm,dir)
        CALL halo_definitions(hmod,cosm,dir)
        CALL halo_properties(hmod,dir)

        !Do the calculation
        DO j=0,3

           IF(j==0) THEN
              !DMONLY
              ip(1)=-1
              ip(2)=-1
              outfile='data/power_'//TRIM(red)//TRIM(ext)
           ELSE IF(j==1) THEN
              !matter - matter
              ip(1)=0
              ip(2)=0
              outfile='dd'
           ELSE IF(j==2) THEN
              !matter - electron pressure
              ip(1)=0
              ip(2)=6
              outfile='dp'
           ELSE IF(j==3) THEN
              !electron pressure - electron pressure
              ip(1)=6
              ip(2)=6
              outfile='pp'
           END IF

           IF(j .NE. 0) outfile=TRIM(base)//TRIM(outfile)//'_'//TRIM(red)//TRIM(ext)

           WRITE(*,fmt='(3I5,A30)') j, ip(1), ip(2), TRIM(outfile)

           CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose=.FALSE.,response=.FALSE.)
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose=.FALSE.)

        END DO

        WRITE(*,*)

     END DO

  ELSE IF(imode==6) THEN

     !n(z) normalisation check

     WRITE(*,*) 'HMx_DRIVER: Checking n(z) functions'
     WRITE(*,*)

     ! Number of n(z) to check
     nnz=11
     DO i=1,nnz
        IF(i==1)  nz=1
        IF(i==2)  nz=4
        IF(i==3)  nz=5
        IF(i==4)  nz=6
        IF(i==5)  nz=7
        IF(i==6)  nz=8
        IF(i==7)  nz=9
        IF(i==8)  nz=11
        IF(i==9)  nz=12
        IF(i==10) nz=13
        IF(i==11) nz=14
        WRITE(*,*) 'HMx_DRIVER: n(z) number:', nz
        WRITE(*,*)
        CALL get_nz(nz,lens)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (linear):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,1)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (quadratic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,2)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (cubic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,2,lens%nnz,3)
        WRITE(*,*)
     END DO

  ELSE IF(imode==7 .OR. imode==8 .OR. imode==9 .OR. imode==10 .OR. imode==11 .OR. imode==37 .OR. imode==42 .OR. imode==43 .OR. imode==47) THEN

     ! General stuff for all cross correlations
     !  7 - Do general angular cross correlation
     !  8 - Angular cross correlation as a function of cosmology
     !  9 - Breakdown angular correlations in halo mass
     ! 10 - Breakdown angular correlations in redshift
     ! 11 - Breakdown angular correlations in halo radius
     ! 36 - CFHTLenS angular correlations
     ! 42 - PAPER: breakdown lensing-lensing per ell
     ! 43 - PAPER: breakdown lensing-y per ell
     ! 47 - Make CMB-lensing data to compare with CAMB
     
     ! Set the fields
     ix=-1
     IF(imode==37) ix=4  ! CFHTLenS autospectrum
     IF(imode==42) ix=11 ! KiDS-450 z=0.1->0.9 autospectrum
     IF(imode==47) ix=3  ! CMB lensing autospectrum
     IF(imode==43) THEN
        ix(1)=11 ! KiDS-450 z=0.1->0.9
        ix(2)=2 ! Compton y
     END IF
     CALL set_xcorr_type(ix,ip)
     IF(imode==37 .OR. imode==47) ip=-1 ! Set DMONLY haloes

     ! Assign the cosmological model
     IF(imode==37 .OR. imode==42 .OR. imode==43) icosmo=4 ! WMAP9
     IF(imode==47) icosmo=26 ! Boring with CAMB linear spectrum
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)

     ! Set the k range
     kmin=1e-3
     kmax=1e1
     nk=128
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Set the a range
     amin=0.1 ! Problems with one-halo term if amin is less than 0.1
     amax=1.
     na=16
     CALL fill_array(amin,amax,a,na)

     ! Need to call 'comoving_distance' at least once first so as to stop
     ! a write trying to happen while printing to screen
     spam=comoving_distance(1.,cosm) ! CARE: This needs to be called before the write-to-screen below
     lmax=5000.
     WRITE(*,*) 'HMx_DRIVER: lmax:', lmax
     WRITE(*,*) '======================================================='
     WRITE(*,*) '            a             z     r [Mpc/h]     k [h/MPc]'
     WRITE(*,*) '======================================================='
     DO i=1,na
        WRITE(*,fmt='(4F14.4)') a(i), redshift_a(a(i)), comoving_distance(a(i),cosm), k_ell(lmax,a(i),cosm)
     END DO
     WRITE(*,*) '======================================================='
     WRITE(*,*)

     ! Assign arrays for power spectra
     ALLOCATE(powa_lin(nk,na),powa_1h(nk,na),powa_2h(nk,na),powa_full(nk,na),powa(nk,na))

     ! Set the ell range
     lmin=1
     lmax=1e4 ! Problems if this is pushed up to 10^5
     nl=128

     ! Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     ! Set the angular arrays in degrees and allocate arrays for theta and xi(theta)
     IF(ixi) THEN
        thmin=0.01
        thmax=10.
        nth=128
        CALL fill_array(log(thmin),log(thmax),theta,nth)
        theta=exp(theta)
        ALLOCATE(xi(3,nth))
     END IF

     ! Output directory
     dir='data/'

     ! Write to screen
     WRITE(*,*) 'HMx_DRIVER: Cross-correlation information'
     WRITE(*,*) 'HMx_DRIVER: output directiory: ', TRIM(dir)
     WRITE(*,*) 'HMx_DRIVER: Profile type 1: ', TRIM(halo_type(ip(1)))
     WRITE(*,*) 'HMx_DRIVER: Profile type 2: ', TRIM(halo_type(ip(2)))
     WRITE(*,*) 'HMx_DRIVER: cross-correlation type 1: ', TRIM(xcorr_type(ix(1)))
     WRITE(*,*) 'HMx_DRIVER: cross-correlation type 2: ', TRIM(xcorr_type(ix(2)))
     WRITE(*,*) 'HMx_DRIVER: P(k) minimum k [h/Mpc]:', REAL(kmin)
     WRITE(*,*) 'HMx_DRIVER: P(k) maximum k [h/Mpc]:', REAL(kmax)
     WRITE(*,*) 'HMx_DRIVER: minimum a:', REAL(amin)
     WRITE(*,*) 'HMx_DRIVER: maximum a:', REAL(amax)
     WRITE(*,*) 'HMx_DRIVER: minimum ell:', REAL(lmin)
     WRITE(*,*) 'HMx_DRIVER: maximum ell:', REAL(lmax)
     WRITE(*,*)     

     IF(imode==7 .OR. imode==37 .OR. imode==42 .OR. imode==43 .OR. imode==47) THEN

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        ! Initialise the lensing part of the calculation
        !CALL write_distances(cosm)

        ! Write out diagnostics
        IF(imode==37 .OR. imode==47) ihm=1 ! HMcode (2016)
        IF(imode==42 .OR. imode==43) ihm=20 ! Standard halo-model in response
        CALL assign_halomod(ihm,hmod,verbose)
        CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response=.FALSE.)

        base=TRIM(dir)//'power'
        CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

        ! Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        ! Set the distance range for the Limber integral
        !r1=100.
        !r2=proj%rs
        r1=0.
        r2=maxdist(proj)

        ! Write information to screen
        WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
        WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
        WRITE(*,*) 'HMx_DRIVER: lower limit of Limber integral [Mpc/h]:', REAL(r1)
        WRITE(*,*) 'HMx_DRIVER: upper limit of Limber integral [Mpc/h]:', REAL(r2)
        WRITE(*,*)

        ! Loop over all types of C(l) to create
        DO j=4,1,-1

           IF(ifull .AND. (j .NE. 4)) CYCLE
           !IF(j==3) CYCLE ! Skip the fucking one-halo term

           ! Write information to screen
           IF(j==1) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing linear'
              outfile=TRIM(dir)//'cl_linear.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_linear.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_linear.dat'
              powa=powa_lin
           ELSE IF(j==2) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing 2-halo'
              outfile=TRIM(dir)//'cl_2halo.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_2halo.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_2halo.dat'
              powa=powa_2h
           ELSE IF(j==3) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing 1-halo'
              outfile=TRIM(dir)//'cl_1halo.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_1halo.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_1halo.dat'
              powa=powa_1h
           ELSE IF(j==4) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing full'
              outfile=TRIM(dir)//'cl_full.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_full.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_full.dat'
              powa=powa_full
           ELSE
              STOP 'HMx_DRIVER: Something went wrong'
           END IF

           WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

           ! Actually calculate the C(ell)
           CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
           CALL write_Cell(ell,Cell,nl,outfile)

           IF(j==4) CALL Cell_contribution(r1,r2,k,a,powa,nk,na,proj,cosm)

           IF(ixi) THEN

              ! Set xi output files
              IF(j==1) THEN
                 outfile=TRIM(dir)//'xi_linear.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_linear.dat'
              ELSE IF(j==2) THEN
                 outfile=TRIM(dir)//'xi_2halo.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_2halo.dat'                 
              ELSE IF(j==3) THEN
                 outfile=TRIM(dir)//'xi_1halo.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_1halo.dat'
              ELSE IF(j==4) THEN
                 outfile=TRIM(dir)//'xi_full.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_full.dat'
              ELSE
                 STOP 'HMx_DRIVER: Something went wrong'
              END IF
              
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

              ! Actually calculate the xi(theta)
              CALL calculate_xi(theta,xi,nth,ell,Cell,nl,NINT(lmax))
              CALL write_xi(theta,xi,nth,outfile)

           END IF

        END DO
        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     ELSE IF(imode==8) THEN

        !Assess cross-correlation as a function of cosmology

        !Loop over cosmology
        sig8min=0.7
        sig8max=0.9
        ncos=5! I may have changed this number inadvertantly 
        DO i=1,ncos

           cosm%sig8=progression(sig8min,sig8max,i,ncos)
           CALL init_cosmology(cosm)
           CALL print_cosmology(cosm)

           CALL assign_halomod(ihm,hmod,verbose)
           CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response=.FALSE.)

           !Fill out the projection kernels
           CALL fill_projection_kernels(ix,proj,cosm)
           CALL write_projection_kernels(proj,cosm)

           !Now do the C(l) calculations
           !Set l range, note that using Limber and flat-sky for sensible results lmin to ~10
           CALL fill_array(log(lmin),log(lmax),ell,nl)
           ell=exp(ell)
           IF(ALLOCATED(Cell)) DEALLOCATE(Cell)
           ALLOCATE(Cell(nl))

           !Write to screen
           WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
           WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
           WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
           WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
           WRITE(*,*)

           !Loop over all types of C(l) to create
           dir='data/'
           base=TRIM(dir)//'cosmology_'
           DO j=1,4
              IF(j==1) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) linear'
                 ext='_cl_linear.dat'
                 powa=powa_lin
              ELSE IF(j==2) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) 2-halo'
                 ext='_cl_2halo.dat'
                 powa=powa_2h
              ELSE IF(j==3) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) 1-halo'
                 ext='_cl_1halo.dat'
                 powa=powa_1h
              ELSE IF(j==4) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) full'
                 ext='_cl_full.dat'
                 powa=powa_full
              END IF           
              outfile=number_file(base,i,ext)
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)
              !Actually calculate the C(l)
              CALL calculate_Cell(0.,maxdist(proj),ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)
           END DO
           WRITE(*,*) 'HMx_DRIVER: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==9) THEN

        ! Breakdown cross-correlation in terms of mass

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        !CALL initialise_cosmology(verbose,cosm)
        CALL print_cosmology(cosm)

        ! Initialise the lensing part of the calculation
        !CALL initialise_distances(verbose,cosm)
        !CALL write_distances(cosm)

        !Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        DO i=0,6
           IF(icumulative .EQV. .FALSE.) THEN
              ! Set the mass intervals
              IF(i==0) THEN
                 m1=mmin
                 m2=mmax
              ELSE IF(i==1) THEN
                 m1=mmin
                 m2=1e11
              ELSE IF(i==2) THEN
                 m1=1e11
                 m2=1e12
              ELSE IF(i==3) THEN
                 m1=1e12
                 m2=1e13
              ELSE IF(i==4) THEN
                 m1=1e13
                 m2=1e14
              ELSE IF(i==5) THEN
                 m1=1e14
                 m2=1e15
              ELSE IF(i==6) THEN
                 m1=1e15
                 m2=1e16
              END IF
           ELSE
              !Set the mass intervals
              IF(i==0) THEN
                 m1=mmin
                 m2=mmax
              ELSE IF(i==1) THEN
                 m1=mmin
                 m2=1e11
              ELSE IF(i==2) THEN
                 m1=mmin
                 m2=1e12
              ELSE IF(i==3) THEN
                 m1=mmin
                 m2=1e13
              ELSE IF(i==4) THEN
                 m1=mmin
                 m2=1e14
              ELSE IF(i==5) THEN
                 m1=mmin
                 m2=1e15
              ELSE IF(i==6) THEN
                 m1=mmin
                 m2=1e16
              END IF
           END IF

           !Set the code to not 'correct' the two-halo power for missing
           !mass when doing the calcultion binned in halo mass
           !STOP 'HMx_DRIVER: Extreme caution here, need to set ip2h=0, but it is defined as parameter in HMx.f90'
           !IF((icumulative .EQV. .FALSE.) .AND. i>1) hmod%ip2h=0
           IF((icumulative .EQV. .TRUE.) .AND. i>0) hmod%ip2h=0
           
           WRITE(*,fmt='(A16)') 'HMx_DRIVER: Mass range'
           WRITE(*,fmt='(A16,I5)') 'HMx_DRIVER: Iteration:', i
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_min [Msun/h]:', m1
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_max [Msun/h]:', m2
           WRITE(*,*)

           !Loop over redshifts
           DO j=1,na

              z=redshift_a(a(j))

              !Initiliasation for the halomodel calcualtion
              CALL assign_halomod(ihm,hmod,verbose)
              CALL init_halomod(m1,m2,z,hmod,cosm,verbose)
              CALL print_halomod(hmod,cosm,verbose)
              CALL calculate_halomod(ip,k,nk,z,powa_lin(:,j),powa_2h(:,j),powa_1h(:,j),powa_full(:,j),hmod,cosm,verbose,response=.FALSE.)

              !Write progress to screen
              IF(j==1) THEN
                 WRITE(*,fmt='(A5,A7)') 'i', 'a'
                 WRITE(*,fmt='(A13)') '   ============'
              END IF
              WRITE(*,fmt='(I5,F8.3)') j, a(j)

           END DO
           WRITE(*,fmt='(A13)') '   ============'
           WRITE(*,*)

           dir='data/'
           IF(i==0) THEN
              outfile=TRIM(dir)//'power'
           ELSE
              base=TRIM(dir)//'mass_'
              mid='_'
              ext='_power'
              outfile=number_file2(base,NINT(log10(m1)),mid,NINT(log10(m2)),ext)
           END IF
           WRITE(*,*) 'HMx_DRIVER: File: ', TRIM(outfile)
           !CALL write_power_a(k,a,powa,nk,na,output)

           !Loop over all types of C(l) to create
           base=TRIM(dir)//'mass_'
           mid='_' 
           DO j=1,4

              !Skip the 1-halo C(l) because it takes ages (2017/02/06)
              IF(j==3) CYCLE

              !Set output files
              IF(j==1) THEN
                 powa=powa_lin
                 outfile=TRIM(dir)//'cl_linear.dat'
                 ext='_cl_linear.dat'
              ELSE IF(j==2) THEN
                 powa=powa_2h
                 outfile=TRIM(dir)//'cl_2halo.dat'
                 ext='_cl_2halo.dat'
              ELSE IF(j==3) THEN
                 powa=powa_1h
                 outfile=TRIM(dir)//'cl_1halo.dat'
                 ext='_cl_1halo.dat'
              ELSE IF(j==4) THEN
                 powa=powa_full
                 outfile=TRIM(dir)//'cl_full.dat'
                 ext='_cl_full.dat'
              END IF

              IF(i>0) outfile=number_file2(base,NINT(log10(m1)),mid,NINT(log10(m2)),ext)

              WRITE(*,*) 'HMx_DRIVER: File: ', TRIM(outfile)

              CALL calculate_Cell(0.,maxdist(proj),ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)

           END DO
           WRITE(*,*) 'HMx_DRIVER: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==10) THEN

        ! Break down cross-correlation in terms of redshift

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        IF(imode==42) ihm=20 ! Standard halo-model but response
        CALL assign_halomod(ihm,hmod,verbose)

        CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response=.FALSE.)

        ! Write distances
        !CALL write_distances(cosm)

        ! Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        CALL write_projection_kernels(proj,cosm)

        ! Write to screen
        WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
        WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
        WRITE(*,*)

        zmin=0.
        zmax=1.
        nz=8

        DO i=0,nz

           IF(i==0) THEN
              r1=0.
              r2=maxdist(proj)
           ELSE
              IF(icumulative .EQV. .FALSE.) THEN
                 z1=progression(zmin,zmax,i,nz)
              ELSE
                 z1=zmin
              END IF
              z2=zmin+(zmax-zmin)*float(i)/float(nz)
              a1=scale_factor_z(z1)
              a2=scale_factor_z(z2)
              r1=comoving_distance(a1,cosm)
              r2=comoving_distance(a2,cosm)
           END IF

           WRITE(*,*) 'HMx_DRIVER:', i
           IF(i>0) THEN
              WRITE(*,*) 'HMx_DRIVER: z1:', REAL(z1)
              WRITE(*,*) 'HMx_DRIVER: z2:', REAL(z2)
           END IF
           WRITE(*,*) 'HMx_DRIVER: r1 [Mpc/h]:', REAL(r1)
           WRITE(*,*) 'HMx_DRIVER: r2 [Mpc/h]:', REAL(r2)

           ! Loop over all types of C(l) to create
           dir='data/'
           base=TRIM(dir)//'redshift_'
           mid='_'
           DO j=1,4

              !Set output files
              IF(j==1) THEN
                 ext='_cl_linear.dat'
                 outfile=TRIM(dir)//'cl_linear.dat'
                 powa=powa_lin
              ELSE IF(j==2) THEN
                 ext='_cl_2halo.dat'
                 outfile=TRIM(dir)//'cl_2halo.dat'
                 powa=powa_2h
              ELSE IF(j==3) THEN
                 ext='_cl_1halo.dat'
                 outfile=TRIM(dir)//'cl_1halo.dat'
                 powa=powa_1h
              ELSE IF(j==4) THEN
                 ext='_cl_full.dat'
                 outfile=TRIM(dir)//'cl_full.dat'
                 powa=powa_full
              END IF

              IF(i>0 .AND. (icumulative .EQV. .FALSE.)) THEN
                 outfile=number_file2(base,i-1,mid,i,ext)
              ELSE IF(i>0 .AND. icumulative) THEN
                 outfile=number_file2(base,0,mid,i,ext)
              END IF
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

              !This crashes for the low r2 values for some reason
              !Only a problem if lmax ~ 10^5
              !STOP 'This crashes for the low r2 values for high ell for some reason - should debug'
              CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)
              CALL write_Cell(ell,Cell,nl,outfile)

           END DO
           WRITE(*,*)

        END DO

        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     ELSE IF(imode==11) THEN

        STOP 'HMx_DRIVER: Error, breakdown in radius is not supported yet'

     ELSE

        STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

     END IF

  ELSE IF(imode==12 .OR. imode==44) THEN

     ! Triad stuff (uses xcorr)
     ! 12 - Project triad
     ! 44 - Project triad for paper

     ! Directory for data output
     dir='data'

     ! Assign the cosmological model
     IF(imode==44) icosmo=4 ! WMAP 9 - BAHAMAS
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     IF(imode==44) ihm=18 ! AGN 7.8
     CALL assign_halomod(ihm,hmod,verbose)

     ! Initialise the lensing part of the calculation
     !CALL write_distances(cosm)

     ! Set the ell range
     lmin=100.
     lmax=4000.
     nl=64

     ! Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     WRITE(*,*) 'HMx_DRIVER: Cross-correlation information'
     WRITE(*,*) 'HMx_DRIVER: output directiory: ', TRIM(dir)
     WRITE(*,*) 'HMx_DRIVER: minimum ell:', REAL(lmin)
     WRITE(*,*) 'HMx_DRIVER: maximum ell:', REAL(lmax)
     WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
     WRITE(*,*)

     ! Set to do the new or old triad
     ! 1 - Old triad
     ! 2 - New triad
     ! 3 - Newer triad
     itriad=3

     IF(itriad==1) THEN
        ntriad=3
     ELSE IF(itriad==2) THEN
        ntriad=7
     ELSE IF(itriad==3) THEN
        ntriad=13
     ELSE
        STOP 'HMX_DRIVER: Error, itriad specified incorrectly'
     END IF

     ! Loop over the triad (septad?)
     ! Original triad is 1,2,3 extended is tomographic lensing
     DO i=1,ntriad

        IF(i==1) THEN           
           IF(itriad==1) THEN
              ix(1)=5 ! KiDS (z = 0.1 -> 0.9)
              ix(2)=3 ! CMB  
              outfile=TRIM(dir)//'/triad_Cl_gal-CMB.dat'
           ELSE IF(itriad==2 .OR. itriad==3) THEN
              ix(1)=11 ! KiDS 450 (z = 0.1 -> 0.9)
              ix(2)=3 ! CMB
              outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.9-CMB.dat'
           ELSE
              STOP 'HMX_DRIVER: Error, itriad specified incorrectly'
           END IF
        ELSE IF(i==2) THEN
           ix(1)=3 ! CMB
           ix(2)=2 ! y
           outfile=TRIM(dir)//'/triad_Cl_CMB-y.dat'
        ELSE IF(i==3) THEN
           IF(itriad==1) THEN
              ix(1)=2 ! y
              ix(2)=5 ! KiDS (z = 0.1 -> 0.9)
              outfile=TRIM(dir)//'/triad_Cl_y-gal.dat'
           ELSE IF(itriad==2 .OR. itriad==3) THEN
              ix(1)=2 ! y
              ix(2)=11 ! KiDS 450 (z = 0.1 -> 0.9)
              outfile=TRIM(dir)//'/triad_Cl_y-gal_z0.1-0.9.dat'
           ELSE
              STOP 'HMX_DRIVER: Error, itriad specified incorrectly'
           END IF
        ELSE IF(i==4) THEN
           ix(1)=12 ! KiDS 450 (z = 0.1 -> 0.5)
           ix(2)=3 ! CMB
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.5-CMB.dat'
        ELSE IF(i==5) THEN
           ix(1)=13 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=3 ! CMB
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-CMB.dat'
        ELSE IF(i==6) THEN
           ix(1)=2 ! y
           ix(2)=12 ! KiDS 450 (z = 0.1 -> 0.5)
           outfile=TRIM(dir)//'/triad_Cl_y-gal_z0.1-0.5.dat'
        ELSE IF(i==7) THEN
           ix(1)=2 ! y
           ix(2)=13 ! KiDS 450 (z = 0.5 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_y-gal_z0.5-0.9.dat'
        ELSE IF(i==8) THEN
           ix(1)=12 ! KiDS 450 (z = 0.1 -> 0.5)
           ix(2)=12 ! KiDS 450 (z = 0.1 -> 0.5)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.5-gal_z0.1-0.5.dat'
        ELSE IF(i==9) THEN
           ix(1)=12 ! KiDS 450 (z =  0.1 -> 0.5)
           ix(2)=11 ! KiDS 450 (z =  0.1 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.5-gal_z0.1-0.9.dat'
        ELSE IF(i==10) THEN
           ix(1)=11 ! KiDS 450 (z = 0.1 -> 0.9)
           ix(2)=11 ! KiDS 450 (z = 0.1 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.9-gal_z0.1-0.9.dat'
        ELSE IF(i==11) THEN
           ix(1)=13 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=12 ! KiDS 450 (z = 0.1 -> 0.5)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-gal_z0.1-0.5.dat'
        ELSE IF(i==12) THEN
           ix(1)=13 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=11 ! KiDS 450 (z = 0.1 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-gal_z0.1-0.9.dat'
        ELSE IF(i==13) THEN
           ix(1)=13 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=13 ! KiDS 450 (z = 0.5 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-gal_z0.5-0.9.dat'
        END IF

        ! Do the cross correlation
        CALL xcorr(ix,mmin,mmax,ell,Cell,nl,hmod,cosm,verbose)
        CALL write_Cell(ell,Cell,nl,outfile)

        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)
        
     END DO

  ELSE IF(imode==13) THEN

     ! 13 - Calculate a cross-correlation coefficient (uses xcorr)

     ! Assign the cosmology
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)

     CALL assign_halomod(ihm,hmod,verbose)

     ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL print_cosmology(cosm)

     ! Initialise the lensing part of the calculation
     !CALL write_distances(cosm)

     ! Set the ell range and allocate arrays for l and C(l)
     lmin=1e0
     lmax=1e4 ! Errors if this is increased to 10^5
     nl=64 
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cell(nl))

     dir='data'

     ixx=-1
     CALL set_xcorr_type(ixx,ip)

     DO i=1,3
        
        IF(i==1) THEN
           ix(1)=ixx(1)
           ix(2)=ixx(1)
           outfile=TRIM(dir)//'/cl_first.dat'
        ELSE IF(i==2) THEN
           ix(1)=ixx(2)
           ix(2)=ixx(2)
           outfile=TRIM(dir)//'/cl_second.dat'
        ELSE IF(i==3) THEN
           ix(1)=ixx(1)
           ix(2)=ixx(2)
           outfile=TRIM(dir)//'/cl_full.dat'
        END IF

        ! Do the cross correlation
        CALL xcorr(ix,mmin,mmax,ell,Cell,nl,hmod,cosm,verbose)
        CALL write_Cell(ell,Cell,nl,outfile)
        
     END DO

  ELSE IF(imode==14 .OR. imode==33) THEN

     ! Make power spectra variations as a function of baryon parameter variations

     !Number of values to try for each parameter
     m=9

     !Set number of k points and k range (log spaced)
     kmin=1e-3
     kmax=1e1
     nk=128
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     !Set the redshift
     z=0.

     !Assigns the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)

     !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     !CALL initialise_cosmology(verbose,cosm)
     CALL print_cosmology(cosm)

     !Initiliasation for the halo-model calcualtion
     IF(imode==33) ihm=3
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     !DMONLY
     ip=-1
     outfile='variations/DMONLY.dat'
     CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)
     CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     !Prevents warning
     ilog=.FALSE.

     !Number of parameters
     npa=6

     !Loop over parameters     
     DO ipa=1,npa
        
        !Set maximum and minimum parameter values and linear or log range
        IF(ipa==1) THEN
           !alpha - virial temperature pre factor
           param_min=0.05
           param_max=0.65
           ilog=.FALSE.
        ELSE IF(ipa==2) THEN
           !epsilon - concentration change due to gas
           param_min=0.5
           param_max=2.0
           ilog=.TRUE.
        ELSE IF(ipa==3) THEN
           !Gamma - KS polytropic index
           param_min=1.12
           param_max=1.22
           ilog=.FALSE.
        ELSE IF(ipa==4) THEN
           !M0 - bound gas transition in Msun/h
           param_min=1e13
           param_max=1e15
           ilog=.TRUE.
        ELSE IF(ipa==5) THEN
           !A* - Stellar mass fraction
           param_min=0.02
           param_max=0.04
           ilog=.FALSE.
        ELSE IF(ipa==6) THEN
           !WHIM temperature in K
           param_min=1e5
           param_max=1e7
           ilog=.TRUE.
        END IF

        !Loop over parameter values
        DO i=1,m

           !Set the parameter value that is being varied
           IF(ilog) THEN
              param=progression_log(param_min,param_max,i,m)
           ELSE
              param=progression(param_min,param_max,i,m)
           END IF

           CALL assign_halomod(ihm,hmod,verbose)

           IF(hmod%fixed_HMx) THEN
              IF(ipa==1) hmod%alpha=param
              IF(ipa==2) hmod%eps=param
              IF(ipa==3) hmod%Gamma=param
              IF(ipa==4) hmod%M0=param
              IF(ipa==5) hmod%Astar=param
              IF(ipa==6) hmod%Twhim=param
           ELSE
              IF(ipa==1) THEN
                 hmod%A_alpha=0.
                 hmod%B_alpha=0.
                 hmod%C_alpha=0.
                 hmod%D_alpha=param
              ELSE IF(ipa==2) THEN
                 hmod%A_eps=0.
                 hmod%B_eps=0.
                 hmod%C_eps=0.
                 hmod%D_eps=log10(param)
              ELSE IF(ipa==3) THEN
                 hmod%A_Gamma=0.
                 hmod%B_Gamma=0.
                 hmod%C_Gamma=0.
                 hmod%D_Gamma=param
              ELSE IF(ipa==4) THEN
                 hmod%A_M0=0.
                 hmod%B_M0=0.
                 hmod%C_M0=0.
                 hmod%D_M0=log10(param)
              ELSE IF(ipa==5) THEN
                 hmod%A_Astar=0.
                 hmod%B_Astar=0.
                 hmod%C_Astar=0.
                 hmod%D_Astar=param
              ELSE IF(ipa==6) THEN
                 hmod%A_Twhim=0.
                 hmod%B_Twhim=0.
                 hmod%C_Twhim=0.
                 hmod%D_Twhim=log10(param)
              END IF
           END IF
           
           CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
           CALL print_halomod(hmod,cosm,verbose)

           !DO NOT DELETE THIS
           !It is only used to print values to the screen later
           !For example, mass, which is inconvenient if written out in full
           IF(ilog) THEN
              param_neat=log10(param)
           ELSE
              param_neat=param
           END IF

           !Write out halo matter and electron-pressure profile information
           !All the string crap is in the loop for a reason
           DO j=10,16
              base='variations/profile_mass_'
              ext='_param_'
              base=number_file(base,j,ext)
              mid='_value_'
              ext='.dat'
              outfile=number_file2(base,ipa,mid,i,ext)
              mass=10.**j
              CALL write_halo_profiles(mass,hmod,cosm,outfile)
           END DO

           !Write out halo mass fraction information
           base='variations/mass_fractions_param_'
           outfile=number_file2(base,ipa,mid,i,ext)
           CALL write_mass_fractions(hmod,cosm,outfile)

           !File base and extension
           base='variations/power_param_'
           mid='_value_'

           !Do the calculation
           DO j=1,9

              IF(j==1) THEN
                 ! matter - matter
                 ip(1)=0
                 ip(2)=0
                 ext='_dd.dat'
              ELSE IF(j==2) THEN
                 ! matter - electron pressure
                 ip(1)=0
                 ip(2)=6
                 ext='_dp.dat'
              ELSE IF(j==3) THEN
                 ! electron pressure - electron pressure
                 ip(1)=6
                 ip(2)=6
                 ext='_pp.dat'
              ELSE IF(j==4) THEN
                 ! matter - CDM
                 ip(1)=0
                 ip(2)=1
                 ext='_dc.dat'
              ELSE IF(j==5) THEN
                 ! CDM - CDM
                 ip(1)=1
                 ip(2)=1
                 ext='_cc.dat'
              ELSE IF(j==6) THEN
                 ! matter - gas
                 ip(1)=0
                 ip(2)=2
                 ext='_dg.dat'
              ELSE IF(j==7) THEN
                 ! gas - gas
                 ip(1)=2
                 ip(2)=2
                 ext='_gg.dat'
              ELSE IF(j==8) THEN
                 ! matter - star
                 ip(1)=0
                 ip(2)=3
                 ext='_ds.dat'
              ELSE IF(j==9) THEN
                 ! star - star
                 ip(1)=3
                 ip(2)=3
                 ext='_ss.dat'
              END IF

              !Set output file
              outfile=number_file2(base,ipa,mid,i,ext)

              !Write progress to screen
              WRITE(*,fmt='(4I5,F14.7,A50)') ipa, i, ip(1), ip(2), param_neat, TRIM(outfile)

              !Do the halo-model calculation and write to file
              CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose=.FALSE.,response=.FALSE.)
              CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose=.FALSE.)
              !Loop over k values
              !Fill the tables with one- and two-halo terms as well as total
              !OPEN(7,file=outfile)
              !DO ii=1,nk       
              !   WRITE(7,fmt='(6ES20.10)') k(ii), pow_lin(ii), pow_2h(ii), pow_1h(ii), pow_full(ii), param
              !END DO
              !CLOSE(7)
              
           END DO

        END DO

     END DO

  ELSE IF(imode==17) THEN

     ! 3D spectra for a user choice of fields

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     CALL assign_halomod(ihm,hmod,verbose)

     ! Set number of k points and k range (log spaced)
     ! The range kmin=1e-3 to kmax=1e4 is necessary to compare to HMcode
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

      !Set the number of redshifts and range (linearly spaced) and convert z -> a
     nz=16
     zmin=0.
     zmax=4.
     CALL fill_array(zmin,zmax,a,nz)
     a=1./(1.+a)
     na=nz

     ! Choose the field types
     DO i=1,2
        WRITE(*,*) 'HMx_driver: Choose halo', i
        CALL set_halo_type(ip(i))
     END DO

     ! User chooses halo model   
     CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response=.FALSE.)

     base='data/power'
     CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

  ELSE IF(imode==18 .OR. imode==48) THEN

     ! Create 3D bias function
     ! 18 - User choice
     ! 48 - HI bias

     ! Assign the cosmological model
     IF(imode==48) icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Assign the halo model
     !IF(imode==48) ihm=3 ! This is what I sent Kiyo & Richard
     IF(imode==48) ihm=25 ! Villaescua-Navarro halo model
     CALL assign_halomod(ihm,hmod,verbose)

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Set the number of redshifts and range (linearly spaced) and convert z -> a
     IF(imode==18) THEN
        nz=16
        zmin=0.
        zmax=4.
     ELSE IF(imode==48) THEN
        nz=6
        zmin=0.
        zmax=5.
     END IF
     CALL fill_array(zmin,zmax,a,nz)
     a=1./(1.+a)
     na=nz

     ! Select field type for bias study
     IF(imode==18) CALL set_halo_type(ihalo)
     IF(imode==48) ihalo=12 ! HI 

     DO j=1,3

        IF(j==1) THEN
           ! DMONLY-DMONLY
           ip(1)=-1
           ip(2)=-1
           base='data/power_mm'
        ELSE IF(j==2) THEN
           ! DMONLY-field
           ip(1)=ihalo
           ip(2)=-1
           base='data/power_mf'
        ELSE IF(j==3) THEN
           ! field-field
           ip(1)=ihalo
           ip(2)=ihalo
           base='data/power_ff'
        END IF
        
        CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response=.FALSE.)
        CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose)

     END DO

  ELSE IF(imode==19) THEN

     ! Create CCL benchmark data

     ! Set number of k points and k range (log spaced)
     IF(Alonso_k) THEN
        ! This is *almost* kmin=1e-4, kmax=1e2, but minus the end points
        nk=256
        kmin=1.027350768179302566e-04
        kmax=9.733773809039202263e+01
     ELSE
        nk=128
        kmin=1e-3
        kmax=1e1
     END IF
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Set the halo model to that for CCL tests
     ihm=9

     ! Loop over tests
     DO j=1,3

        ! Assigns the cosmological model
        IF(j==1) icosmo=1
        IF(j==2) icosmo=2
        IF(j==3) icosmo=3
        CALL assign_cosmology(icosmo,cosm,verbose)
        CALL init_cosmology(cosm)
        CALL print_cosmology(cosm)

        ! Loop over redshifts
        DO i=1,2

           ! Sets the redshift and file names
           IF(i==1) THEN
              z=0.
              IF(j==1) outfile='/Users/Mead/Physics/HMx/CCL/HMx_power_model1_z0.txt'
              IF(j==2) outfile='/Users/Mead/Physics/HMx/CCL/HMx_power_model2_z0.txt'
              IF(j==3) outfile='/Users/Mead/Physics/HMx/CCL/HMx_power_model3_z0.txt'
           ELSE IF(i==2) THEN
              z=1.
              IF(j==1) outfile='/Users/Mead/Physics/HMx/CCL/HMx_power_model1_z1.txt'
              IF(j==2) outfile='/Users/Mead/Physics/HMx/CCL/HMx_power_model2_z1.txt'
              IF(j==3) outfile='/Users/Mead/Physics/HMx/CCL/HMx_power_model3_z1.txt'
           END IF

           ! Initialise the halo-model calculation
           CALL assign_halomod(ihm,hmod,verbose)
           CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
           CALL print_halomod(hmod,cosm,verbose)

           ! Do the halo-model calculation
           ip=-1
           CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)

           ! Write out the results
           CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

        END DO

     END DO


  ELSE IF(imode==20) THEN

     ! Ma et al. Fig. 1

     ! Set the cosmology
     icosmo=3
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the halo model
     z=0.
     ihm=4
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     ! Make the Figure
     CALL YinZhe_Fig1(hmod,cosm)

  ELSE IF(imode==21) THEN

     ! Stuff for diagnosing problems with the window function integrand
     
     outfile='winint/integrand.dat'
     irho=11
     rv=1.
     c=4.
     rs=rv/c
     p1=1.18
     p2=0.
     rmin=0.
     rmax=rv
     CALL winint_diagnostics(rmin,rmax,rv,rs,p1,p2,irho,outfile)

  ELSE IF(imode==22) THEN

     ! Speed tests for W(M,k) integrals

     !k range
     kmin=1e-2
     kmax=1e3
     nk=512
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     !Halo parameters
     rv=1.
     c=4.
     rs=rv/c
     p1=1.5
     p2=0.
     irho=11
     rmin=0.
     rmax=rv

     CALL winint_speed_tests(k,nk,rmin,rmax,rv,rs,p1,p2,irho)

  ELSE IF(imode==23) THEN

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e1
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Set the redshift
     z=0.

     ! Directory for output
     dir='Mead2017'

     ! Set the halo model
     ihm=12

     ! Loop over cosmologies
     DO i=1,9

        IF(i==1) THEN
           ! LCDM
           icosmo=1
           outfile='LCDM'
        ELSE IF(i==2) THEN
           ! OCDM
           icosmo=5
           outfile='OCDM'
        ELSE IF(i==3) THEN
           ! EdS
           icosmo=15
           outfile='SCDM'
        ELSE IF(i==4) THEN
           ! w = -0.7
           icosmo=16
           outfile='w-0.7'
        ELSE IF(i==5) THEN
           ! w = -1.3
           icosmo=17
           outfile='w-1.3'
        ELSE IF(i==6) THEN
           ! wa = 0.5
           icosmo=18
           outfile='wa0.5'
        ELSE IF(i==7) THEN
           ! wa = -0.5
           icosmo=19
           outfile='wa-0.5'
        ELSE IF(i==8) THEN
           ! w = -0.7; wa = -1.5
           icosmo=20
           outfile='w0-0.7wa-1.5'
        ELSE IF(i==9) THEN
           ! w = -1.3; wa = 0.5
           icosmo=21
           outfile='w0-1.3wa0.5'
        END IF

        ! Assigns the cosmological model
        CALL assign_cosmology(icosmo,cosm,verbose)
        CALL init_cosmology(cosm)

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        ! Initiliasation for the halomodel calcualtion
        CALL assign_halomod(ihm,hmod,verbose)
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Do the halo-model calculation
        ip=-1 ! DMONLY
        CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)

        ! Write out the results
        outfile=TRIM(dir)//'/power_'//TRIM(outfile)//'.dat'
        CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     END DO

  ELSE IF(imode==24) THEN

     ! Assign the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Assign the base halo model
     CALL assign_halomod(ihm,hmod,verbose)

     ! Choose fitting to do
     ! 1 - Star autospectrum
     ! 2 - Gas autospectrum
     ifit=1

     ! Default number of refinements
     nref=3

     np=3
     ALLOCATE(pmin(np),pmax(np),pbest(np),pnow(np),prange(np),plog(np))
     ALLOCATE(nrange(np),iii(np),ipbest(np))

     !Set the redshift
     z=0.0

     ! Set the model
     model='AGN_TUNED_nu0'
     !model='AGN_7p6_nu0'
     !model='AGN_8p0_nu0'

     ! Default number of parameters
     nrange=11
    
     IF(ifit==1) THEN

        ! Star power spectra
        
        ! Halo type
        ip=3

        ! Astar (height of f*(M) distribution)
        pmin(1)=1e-3
        pmax(1)=1e-1
        plog(1)=.TRUE.

        ! Mstar (location of max of f*(M) distribution)
        pmin(2)=1e11
        pmax(2)=1e13
        plog(2)=.TRUE.

        ! sigma (width of f*(M) distribution)
        pmin(3)=1.2
        pmax(3)=2.0
        plog(3)=.TRUE.
        nrange(3)=1

     ELSE IF(ifit==2) THEN

        ! Gas power spectra

        ! Halo type
        ip=2

        ! Gamma
        pmin(1)=1.05
        pmax(1)=1.55
        plog(1)=.FALSE.

        ! M0
        pmin(2)=1e13
        pmax(2)=1e15
        plog(2)=.TRUE.

        ! epsilon
        pmin(3)=1e-1
        pmax(3)=1e1
        plog(3)=.TRUE.
        nrange(3)=1

     ELSE

        STOP 'FITTING: Error, ifit specified incorrectly'
        
     END IF

     ! Read in the data
     infile=BAHAMAS_power_file_name(model,z,ip)
     CALL read_simulation_power_spectrum(k,pow_sim,nk,infile,cut_nyquist=.TRUE.,subtract_shot=.TRUE.,verbose=.TRUE.)
     
     ! Fix the paramter range
     DO j=1,np
        IF(plog(j)) THEN
           prange(j)=log(pmax(j)/pmin(j))
        ELSE
           prange(j)=pmax(j)-pmin(j)
        END IF
     END DO

     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Loop over refinements
     i=0
     refine=.FALSE.
     DO

        ! Always set this to be a HUGE number at the beginning of any refinement level
        fom_best=HUGE(fom)

        ! Update the range to be centred on the best fit
        ! Maybe also decrease the parameter space
        IF(refine) THEN
           DO j=1,np
              IF(nrange(j)>1) THEN
                 IF(plog(j)) THEN
                    pmin(j)=exp(log(pbest(j))-prange(j)/(2.**(i+1)))
                    pmax(j)=exp(log(pbest(j))+prange(j)/(2.**(i+1)))
                 ELSE
                    pmin(j)=pbest(j)-prange(j)/(2.**(i+1))
                    pmax(j)=pbest(j)+prange(j)/(2.**(i+1))
                 END IF
              END IF
           END DO
        END IF

        ! Write data to screen
        WRITE(*,*) 'FITTING: Refinement', i
        WRITE(*,*) '======================================================================'
        WRITE(*,*) ' Refinement   Parameter              minimum                   maximum'
        WRITE(*,*) '======================================================================'
        DO j=1,np
           IF(nrange(j)>1) THEN
              IF(plog(j)) THEN
                 WRITE(*,*) i, j, log10(pmin(j)), log10(pmax(j))
              ELSE
                 WRITE(*,*) i, j, pmin(j), pmax(j)
              END IF
           END IF
        END DO
        WRITE(*,*) '======================================================================'     
        WRITE(*,*)

        ! Loop over parameters
        DO i1=1,nrange(1)
           DO i2=1,nrange(2)
              DO i3=1,nrange(3)

                 ! Update the parameters
                 iii(1)=i1
                 iii(2)=i2
                 iii(3)=i3
                 DO j=1,np
                    IF(plog(j)) THEN
                       pnow(j)=progression_log(pmin(j),pmax(j),iii(j),nrange(j))
                    ELSE
                       pnow(j)=progression(pmin(j),pmax(j),iii(j),nrange(j))
                    END IF
                 END DO

                 ! Update fitted parameters
                 IF(ifit==1) THEN
                    hmod%Astar=pnow(1)
                    hmod%Mstar=pnow(2)
                    hmod%sstar=pnow(3)
                 ELSE IF(ifit==2) THEN
                    hmod%gamma=pnow(1)
                    hmod%M0=pnow(2)
                    hmod%eps=pnow(3)
                 ELSE
                    STOP 'FITTING: Error, ifit specified incorrectly'
                 END IF

                  ! Initiliasation for the halomodel calcualtion
                 CALL init_halomod(mmin,mmax,z,hmod,cosm,.FALSE.)

                 ! Do the halo-model calculation
                 CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

                 ! Calculate the figure-of-merit and update best fit
                 fom=figure_of_merit(pow_full,pow_sim,nk)
                 IF(fom<fom_best) THEN
                    ipbest=iii
                    pbest=pnow
                    fom_best=fom
                 END IF

              END DO
           END DO
        END DO

        ! Write best-fit to screen
        WRITE(*,*) 'FITTING: Best fit'
        WRITE(*,*) '============================================'
        WRITE(*,*) '  Parameter    location                value'
        WRITE(*,*) '============================================'
        DO j=1,np
           IF(nrange(j)>1) THEN
              IF(plog(j)) THEN
                 WRITE(*,*) j, ipbest(j), log10(pbest(j)), pbest(j)
              ELSE
                 WRITE(*,*) j, ipbest(j), pbest(j)
              END IF
           END IF
        END DO
        WRITE(*,*) '============================================'
        WRITE(*,*)

        ! Check that the best fit is not on the edge of the range
        refine=.TRUE.
        DO j=1,np
           IF(nrange(j)>1) THEN
              IF(ipbest(j)==1 .OR. ipbest(j)==nrange(j)) THEN
                 WRITE(*,*) 'FITTING: Parameter on edge of range'
                 WRITE(*,*) 'FITTING: Parameter:', j
                 IF(plog(j)) THEN
                    WRITE(*,*) 'FITTING: Value:', log10(pbest(j)), pbest(j)
                 ELSE
                    WRITE(*,*) 'FITTING: Value:', pbest(j)
                 END IF
                 WRITE(*,*) 'FITTING: Number of points:', nrange(j)
                 WRITE(*,*) 'FITTING: Best fit location:', ipbest(j)
                 WRITE(*,*)
                 refine=.FALSE.
              END IF
           END IF
        END DO

        ! Update refinement level or exit
        IF(refine .AND. i==nref) THEN
           EXIT
        ELSE IF(refine) THEN
           i=i+1
        END IF
        refine=.TRUE.

     END DO

     ! Initiliasation for the halomodel calcualtion
     CALL assign_halomod(ihm,hmod,verbose)

     ! Set to the best fit
     IF(ifit==1) THEN
        hmod%Astar=pbest(1)
        hmod%Mstar=pbest(2)
        hmod%sstar=pbest(3)
     ELSE IF(ifit==2) THEN
        hmod%Gamma=pbest(1)
        hmod%M0=pbest(2)
        hmod%eps=pbest(3)
     ELSE
        STOP 'FITTING: Error, ifit specified incorrectly'
     END IF
        
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)
     CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)
     
     ! Write out the results
     !outfile='data/power.dat'
     !CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)
     outfile='fitting/power.dat'
     OPEN(7,file=outfile)
     DO i=1,nk
        WRITE(7,*) k(i), pow_full(i), pow_sim(i)
     END DO
     CLOSE(7)

  ELSE IF(imode==25 .OR. imode==38 .OR. imode==39 .OR. imode==40 .OR. imode==41) THEN

     ! MCMC-like fitting
     ! 25 - Mira Titan
     ! 38 - FrankenEmu
     ! 39 - Hydro
     ! 40 - Random Mira Titan
     ! 41 - Random FrankenEmu

     ! Set the random-number generator
     CALL RNG_set(1)

     ! Number of MCMC points
     n=1000

     ! Allocate array for cosmological models
     nz=4
     IF(imode==25 .OR. imode==40) THEN
        ncos=10 ! Number of Mita Titan nodes - only 10 have Omega_nu = 0.
     ELSE IF(imode==38 .OR. imode==41) THEN
        ncos=37 ! Number of FrankenEmu nodes
     ELSE IF(imode==39) THEN
        ncos=1 ! Only one cosmology for hydro
     ELSE
        STOP 'HMX_DRIVER: Error, something went wrong'
     END IF
     ALLOCATE(cosms(ncos))

     ! MCMC: Set redshifts and halo-profiles here
     
     IF(imode==25 .OR. imode==38 .OR. imode==40 .OR. imode==41) THEN
        ! Mira Titan or FrankenEmu
        nz=4
        ALLOCATE(zs(nz))
        zs(1)=0.0
        zs(2)=0.5
        zs(3)=1.0
        zs(4)=2.0
        ip=-1 ! DMONLY
     ELSE IF(imode==39) THEN
        ! Hydro simulations
        nz=1
        ALLOCATE(zs(nz))
        zs(1)=2.0 ! Set the redshift
        !ip(1)=3   ! Stars
        !ip(2)=3   ! Stars
        !ip(1)=0   ! Matter
        !ip(2)=6   ! Pressure
        ip(1)=2   ! Gas
        ip(2)=2   ! Gas
        name='AGN_TUNED_nu0'
     ELSE
        STOP 'HMX_DRIVER: Error, something went wrong'
     END IF

     ! Assign the cosmological models
     DO i=1,ncos
        IF(imode==40) icosmo=24 ! Random Mira Titan cosmology
        IF(imode==41) icosmo=25 ! Random FrankenEmu cosmology
        IF(imode==25) icosmo=100+i ! Set set Mira Titan node
        IF(imode==38) icosmo=200+i ! Set set FrankenEmu node
        IF(imode==39) icosmo=4 ! WMAP9
        CALL assign_cosmology(icosmo,cosms(i),verbose=.TRUE.)
        CALL init_cosmology(cosms(i))
        CALL print_cosmology(cosms(i))
     END DO

     ! MCMC: Set halo-model here
     ALLOCATE(hmods(ncos))
     IF(imode==25 .OR. imode==38 .OR. imode==40 .OR. imode==41) ihm=15 ! HMcode 2018
     IF(imode==39) ihm=20 ! Standard halo-model response
     DO i=1,ncos
        CALL assign_halomod(ihm,hmods(i),verbose=.FALSE.)
     END DO

     ! Just print one halo model to screen to see 
     CALL init_halomod(mmin,mmax,0.,hmods(1),cosms(1),verbose=.TRUE.)
     CALL print_halomod(hmods(1),cosms(1),verbose=.TRUE.)

     ! Get the simulation power spectrum for the models
     DO i=1,ncos
        DO j=1,nz

           ! Read in power spectra
           IF(imode==25 .OR. imode==40) THEN
              CALL get_Mira_Titan_power(k_sim,pow_sim,nk,zs(j),cosms(i),rebin=.TRUE.)
           ELSE IF(imode==38 .OR. imode==41) THEN
              CALL get_FrankenEmu_power(k_sim,pow_sim,nk,zs(j),cosms(i),rebin=.TRUE.)
           ELSE IF(imode==39) THEN
              CALL get_BAHAMAS_power(k_sim,pow_sim,nk,zs(j),name,ip,cosms(i))
           ELSE
              STOP 'HMX_DRIVER: Error, something went wrong'
           END IF
           
           ! Allocate big arrays for P(k,z,cosm)
           IF(i==1 .AND. j==1) THEN
              ALLOCATE(ks_sim(nk,nz,ncos),pows_sim(nk,nz,ncos))
           END IF
           ks_sim(:,j,i)=k_sim
           pows_sim(:,j,i)=pow_sim
           
        END DO
     END DO

     ! Allocate arrays for halo-model power
     ALLOCATE(pows_full(nk,nz,ncos))

     ! MCMC: Set varying parameters, number of them, and initial values here

     ! Set initial parameters
     IF(imode==25 .OR. imode==38 .OR. imode==40 .OR. imode==41) THEN

        ! Mira Titan or FrankenEmu fitting for HMcode
        np=12
        ALLOCATE(pbest(np),pnew(np),pold(np),prange(np),porig(np),plog(np))
        plog=.FALSE. ! None of these parameters are explored in log
                
        porig(1)=418. ! Dv0
        porig(2)=-0.352 ! Dvp
        porig(3)=1.59 ! dc0
        porig(4)=0.0314 !dcp
        porig(5)=0.603 !eta0
        porig(6)=0.300 !eta1
        !porig(7)=0.0095 ! Mead (2016) damping
        !porig(8)=1.37 ! Mead (2016) damping
        porig(7)=0.188 ! Mead (2015) damping
        porig(8)=4.29 ! Mead (2015) damping
        porig(9)=0.584 ! kstar
        porig(10)=3.13 ! As 
        porig(11)=3.24 ! alpha0
        porig(12)=1.85 ! alpha1
        
     ELSE IF(imode==39) THEN

        ! Hydro fitting
        np=3
        ALLOCATE(pbest(np),pnew(np),pold(np),prange(np),porig(np),plog(np))
        plog=.TRUE. ! All of these parameters are explored in log

        ! stars-stars - complex
        !porig(1)=0.03 ! A_star
        !porig(2)=1e12 ! M_star        
        !porig(3)=10.  ! c_star
        !porig(4)=1.2  ! sigma_star

        ! stars-stars - simple
        !porig(1)=0.03 ! A_star
        !porig(2)=10.  ! c_star

        ! matter-pressure 
        !porig(1)=1e6 ! T_whim

        ! gas-gas
        porig(1)=1.1    ! eps (not one because log(1)=0 and this messes things up in setting the ranges
        porig(2)=1.17   ! Gamma
        plog(2)=.FALSE. ! Gamma not explored in log
        porig(3)=1e14   ! M0      
        
     END IF

     ! Take the log of those parameters that are explored in log
     DO i=1,np
        IF(plog(i)) porig(i)=log10(porig(i))
     END DO

     ! Set the new parameters
     pold=porig
     pnew=porig

     ! Set the ranges (sigma) for the parameters
     CALL set_mcmc_parameter_ranges(imode,ip,prange,pold,plog,np,ks_sim,nk,zs,nz,pows_sim,hmods,cosms,ncos,verbose=.TRUE.)

     ! Set the best figures-of-merit to some huge value
     fom_old=HUGE(fom)
     fom_new=HUGE(fom)
     fom_best=HUGE(fom)

      ! Set the random-number generator
     CALL RNG_set(0)

     ! Loop over number of runs
     WRITE(*,*) 'MCMC: Starting MCMC'
     ii=0
     OPEN(10,file='data/mcmc.dat')
     DO l=1,n

        IF(l==1) THEN
           ! Do nothing
        ELSE IF(l==n) THEN
           ! Set to best-fitting parameters on last go
           pnew=pbest           
        ELSE
           ! Randomly jump parameters
           !IF(jump) CALL set_ranges(prange,delta,pold,np,ks_sim,nk,zs,nz,pows_sim,hmods,cosms,ncos,verbose=.FALSE.)
           DO i=1,np
              pnew(i)=random_Gaussian(pold(i),prange(i))
           END DO
        END IF

        ! Calculate the figure-of-merit
        CALL fom_multiple(imode,ip,fom_new,pnew,plog,np,ks_sim,nk,zs,nz,pows_full,pows_sim,hmods,cosms,ncos)

        ! Write original power spectra to disk
        IF(l==1) THEN

           ! Set the 
           fom_orig=fom_new

           base='fitting/original_cosmo'
           mid='_z'
           ext='.dat'

           ! Write data to disk
           DO i=1,ncos
              DO j=1,nz
                 outfile=number_file2(base,i,mid,j,ext)
                 OPEN(7,file=outfile)
                 DO ik=1,nk
                    WRITE(7,*) ks_sim(ik,j,i), pows_full(ik,j,i), pows_sim(ik,j,i)
                 END DO
                 CLOSE(7)
              END DO
           END DO

        END IF

        WRITE(*,fmt='(I10,F14.7)') l, fom_new
           
        IF(fom_new<fom_best) THEN
           ! If fom is best then always accept...
           pbest=pnew
           ibest=l
           fom_best=fom_new
           jump=.TRUE.
        ELSE IF(fom_new <= fom_old) THEN
           ! ... also accept if fom is better
           jump=.TRUE.
        ELSE IF(fom_old/fom_new<random_uniform(0.,1.)) THEN
           ! ...otherwise accept poorer fom with some probability...
           jump=.TRUE.
        ELSE
           ! ...otherwise, do nothing!
           jump=.FALSE.
        END IF

        IF(jump) THEN
           ii=ii+1
           pold=pnew
           fom_old=fom_new
           WRITE(10,*) fom_old, (pold(j), j=1,np)
        END IF

     END DO
     CLOSE(10)
     WRITE(*,*) 'MCMC: Done'
     WRITE(*,*)

     WRITE(*,*) 'MCMC: Best location:', ibest
     WRITE(*,*) 'MCMC: Total attempts:', n
     WRITE(*,*) 'MCMC: Accepted steps:', ii
     WRITE(*,*) 'MCMC: Fraction accepted:', REAL(ii)/REAL(n)
     WRITE(*,*) 'MCMC: Original figure-of-merit:', fom_orig
     WRITE(*,*) 'MCMC: Final figure-of-merit:', fom_new
     WRITE(*,*)

     WRITE(*,*) 'MCMC: Best-fitting parameters'
     WRITE(*,*) '=========================================='
     WRITE(*,*) 'Parameter            Original         Best'
     WRITE(*,*) '=========================================='
     DO i=1,np
        IF(plog(i)) THEN
           WRITE(*,fmt='(I10,A5,2F14.7)') i, 'lin', 10**porig(i), 10**pbest(i)
        END IF
        WRITE(*,fmt='(I10,A5,2F14.7)') i, 'log', porig(i), pbest(i)
     END DO
     WRITE(*,*) '=========================================='
     WRITE(*,*)

     base='fitting/cosmo'
     mid='_z'
     ext='.dat'

     ! Output data
     DO i=1,ncos
        DO j=1,nz
           outfile=number_file2(base,i,mid,j,ext)
           OPEN(7,file=outfile)
           DO ik=1,nk
              WRITE(7,*) ks_sim(ik,j,i), pows_full(ik,j,i), pows_sim(ik,j,i)
           END DO
           CLOSE(7)
        END DO
     END DO

  ELSE IF(imode==26) THEN

     ! Automated testing

     ! Loop over different tests
     DO itest=1,3

        IF(itest==1) THEN
           infile='benchmarks/power_HMcode_Mead.txt'
           ihm=1
           base='data/power_HMx_Mead'
           icosmo=1
        ELSE IF(itest==2) THEN
           infile='benchmarks/power_HMcode_basic.txt'
           ihm=2
           base='data/power_HMx_basic'
           icosmo=1
        ELSE IF(itest==3) THEN
           infile='benchmarks/power_HMcode_standard.txt'
           ihm=5
           base='data/power_HMx_standard'
           icosmo=1
        END IF

        ! Allocate arrays of k and a
        nk=128
        na=16
        ALLOCATE(k(nk),a(na),powa(nk,na))

        ! Read-in test data
        OPEN(7,file=infile,status='old')
        DO i=0,nk
           IF(i==0) THEN
              READ(7,*) crap, (a(j), j=1,na)
           ELSE
              READ(7,*) k(i), (powa(i,j), j=1,na)
           END IF
        END DO
        CLOSE(7)

        ! Convert read-in z to a; note that 'a' is the correct argument
        ! for the function here because actually read-in a is z.
        DO i=1,na
           a(i)=scale_factor_z(a(i))
        END DO
           
        ! Assigns the cosmological model
        CALL assign_cosmology(icosmo,cosm,verbose_tests)
        CALL init_cosmology(cosm)
        CALL print_cosmology(cosm)

        ! Assign the halo model
        CALL assign_halomod(ihm,hmod,verbose_tests)
        !CALL print_halomod(hmod)
        
        ip=-1 ! Set DMONLY profiles 
        CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose_tests,response=.FALSE.)

        ! Loop over k and a
        error_max=0.
        DO j=1,na
           DO i=1,nk
              error=ABS(-1.+powa_full(i,j)/powa(i,j))
              IF(error>error_max) error_max=error
              IF(error>tolerance) THEN
                 WRITE(*,*) 'HMx_DRIVER: Test:', itest
                 WRITE(*,*) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                 WRITE(*,*) 'HMx_DRIVER: Scale-factor:', a(j)
                 WRITE(*,*) 'HMx_DRIVER: Expected power:', powa(i,j)
                 WRITE(*,*) 'HMx_DRIVER: Model power:', powa_full(i,j)
                 WRITE(*,*) 'HMx_DRIVER: Tolerance:', tolerance
                 WRITE(*,*) 'HMx_DRIVER: Error:', error
                 WRITE(*,*)
                 ifail=.TRUE.
              END IF
           END DO
        END DO

        WRITE(*,*) 'HMx_DRIVER: Test:', itest
        WRITE(*,*) 'HMx_DRIVER: Tolerance:', tolerance
        WRITE(*,*) 'HMx_DRIVER: Max error:', error_max
        WRITE(*,*)

        ! Write data to file
        !base='data/power'
        CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose_tests)

        DEALLOCATE(k,a,powa)
        
     END DO

     IF(ifail) THEN
        STOP 'HMx_DRIVER: Error, tests failed'
     ELSE
        WRITE(*,*) 'HMx_DRIVER: Tests should take around 0.80 seconds to run'
        WRITE(*,*) 'HMx_DRIVER: Tests passed'
        WRITE(*,*)
     END IF

  ELSE IF(imode==27 .OR. imode==28 .OR. imode==29 .OR. imode==30)  THEN

     ! Comparison with FrankenEmu or Mira Titan

     ! Set DMONLY
     ip=-1

     ! Number of cosmological models (+1)
     IF(imode==28 .OR. imode==30) n=37
     IF(imode==27 .OR. imode==29) n=10

     ! Allocate arrays
     nz=4
     na=nz
     ALLOCATE(zs(nz))

     ! Set redshifts/scale factors
     zs(1)=0.0
     zs(2)=0.5
     zs(3)=1.0
     zs(4)=2.0

     base='emulator/cosmo'
     mid='_z'
     ext='.dat'

     ! Initiliasation for the halomodel calcualtion
     CALL assign_halomod(ihm,hmod,verbose)

     ! Loop over cosmologies
     DO i=0,n

        IF(imode==27) icosmo=100+i ! Mira Titan nodes 
        IF(imode==28) icosmo=200+i ! FrankenEmu nodes
        IF(imode==29) icosmo=24 ! Random Mira Titan
        IF(imode==30) icosmo=25 ! Random FrankenEmu
        CALL assign_cosmology(icosmo,cosm,verbose)
        CALL init_cosmology(cosm)
        CALL print_cosmology(cosm)

        ! Loop over redshift
        DO j=1,nz

           IF(imode==27 .OR. imode==29) CALL get_Mira_Titan_power(k_sim,pow_sim,nk,zs(j),cosm,rebin=.FALSE.)
           IF(imode==28 .OR. imode==30) CALL get_FrankenEmu_power(k_sim,pow_sim,nk,zs(j),cosm,rebin=.FALSE.)

           ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow(nk))

           CALL init_halomod(mmin,mmax,zs(j),hmod,cosm,verbose=.FALSE.)
           CALL print_halomod(hmod,cosm,verbose=.FALSE.)
           CALL calculate_halomod(ip,k_sim,nk,zs(j),pow_lin,pow_2h,pow_1h,pow,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

           ! Write data to disk
           outfile=number_file2(base,i,mid,j,ext)
           OPEN(7,file=outfile)
           DO ii=1,nk
              WRITE(7,*) k_sim(ii), pow(ii), pow_sim(ii)
           END DO
           CLOSE(7)

           DEALLOCATE(pow_lin,pow_2h,pow_1h,pow)
           DEALLOCATE(k_sim,pow_sim)
           
        END DO

     END DO

  ELSE IF(imode==31) THEN

     ! Power breakdown as a function of mass (paper)

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Assign the cosmological model
     icosmo=4 ! BAHAMAS cosmology
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the redshift
     z=0.

     ! Set the halo model
     ihm=3
     CALL assign_halomod(ihm,hmod,verbose)

     ! Loop over upper limit of mass integral
     DO i=10,16

        ! Set the upper limit for the mass integration
        ! Needs to be 10. to enforce that m2 is real
        m2=10.**i

        ! Initiliasation for the halomodel calcualtion
        CALL init_halomod(mmin,m2,z,hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Loop over fields
        DO j1=1,2
           DO j2=j1,2

              ! Set the fields
              IF(j1==1) ip(1)=0
              IF(j1==2) ip(1)=6
              IF(j2==1) ip(2)=0
              IF(j2==2) ip(2)=6

              ! Do the halo-model calculation
              CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)

              ! Write out the results
              base='data/power_'
              mid=''
              ext='_m'
              base=number_file2(base,ip(1),mid,ip(2),ext)
              ext='.dat'
              outfile=number_file(base,i,ext)
              CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

           END DO
        END DO

     END DO

  ELSE IF(imode==34) THEN

     ! Write out the variation of HMx parameters with T_AGN and z

     !Assign the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     !Sets the redshift
     z=0.

     !Initiliasation for the halomodel calcualtion
     ihm=18
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     zmin=0.
     zmax=4.
     nz=101

     DO j=1,3
        
        IF(j==1) THEN
           outfile='data/HMx_params_AGN-lo.dat'
           hmod%Theat=10**7.6
        ELSE IF(j==2) THEN
           outfile='data/HMx_params_AGN.dat'
           hmod%Theat=10**7.8
        ELSE IF(j==3) THEN
           outfile='data/HMx_params_AGN-hi.dat'
           hmod%Theat=10**8.0
        END IF

        OPEN(7,file=outfile)
        DO i=1,nz
           hmod%z=progression(zmin,zmax,i,nz)
           WRITE(7,*) hmod%z, HMx_alpha(hmod), HMx_eps(hmod), HMx_Gamma(hmod), HMx_M0(hmod), HMx_Astar(hmod), HMx_Twhim(hmod)
        END DO
        CLOSE(7)
        
     END DO

  ELSE IF(imode==35) THEN

     ! Look at the effect of cores on halo profiles

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Assign the cosmological model
     icosmo=1 ! Boring
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the redshift
     z=0.

     ! Initiliasation for the halomodel calcualtion
     ihm=3
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     ! Do the halo-model calculation
     ip(1)=-1
     ip(2)=-1

      ! Initiliasation for the halomodel calcualtion
     ihm=21
     CALL assign_halomod(ihm,hmod,verbose)     
     CALL print_halomod(hmod,cosm,verbose)

     ! Range of cores to explore
     rcore_min=0.
     rcore_max=0.1
     ncore=16
     
     DO i=1,ncore

        ! Set the core radius
        hmod%rcore=progression(rcore_min,rcore_max,i,ncore)

        ! Initialise the halo-model calculation
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)

        ! Do the halo-model calculation
        CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)

        ! Write out the results
        !outfile='data/power_cored.dat'
        base='data/power_cored_'
        ext='.dat'
        outfile=number_file(base,i,ext)
        CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     END DO

  ELSE IF(imode==36) THEN

     ! Automated testing of hydro models

     ! Assigns the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Assign the halo model
     ihm=3
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)
     
     ! Initially assume all the tests will pass
     ifail=.FALSE.

     ! File naming things
     mid=''
     inext='.txt'
     outext='.dat'

     ! Loop over redshifts
     DO iz=1,4

        ! Set the redshift and input and output file bases
        IF(iz==1) THEN
           z=0.0
           inbase='benchmarks/power_z0.0_'
           outbase='data/power_z0.0_'
        ELSE IF(iz==2) THEN
           z=0.5
           inbase='benchmarks/power_z0.5_'
           outbase='data/power_z0.5_'
        ELSE IF(iz==3) THEN
           z=1.0
           inbase='benchmarks/power_z1.0_'
           outbase='data/power_z1.0_'
        ELSE IF(iz==4) THEN
           z=2.0
           inbase='benchmarks/power_z2.0_'
           outbase='data/power_z2.0_'
        ELSE
           STOP 'HMX_DRIVER: Error, iz specified incorrectly'
        END IF

        ! Initialise the halo-model calculation for this redshift
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose_tests)
        CALL print_halomod(hmod,cosm,verbose_tests)

        ! Loop over different fields
        DO itest=0,6
           DO jtest=itest,6

              ! Cycle for these fields as they are not used
              IF(itest==4 .OR. itest==5) CYCLE
              IF(jtest==4 .OR. jtest==5) CYCLE

              ! Set the input and output files
              infile=number_file2(inbase,itest,mid,jtest,inext)
              outfile=number_file2(outbase,itest,mid,jtest,outext)

              ! Allocate arrays of k and power
              nk=file_length(infile,verbose=.FALSE.)
              ALLOCATE(k(nk),pow(nk))

              ! Read in benchmark data
              OPEN(7,file=infile,status='old')
              DO i=1,nk
                 READ(7,*) k(i), spam, spam, spam, pow(i)
              END DO
              CLOSE(7)

              ! Allocate arrays for the halo-model calculation
              ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

              ! Do the actual halo-model calculation
              ip(1)=itest
              ip(2)=jtest
              CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose_tests,response=.FALSE.)

              ! Loop over k and check for accuracy
              error_max=0.
              DO i=1,nk
                 error=ABS(-1.+pow_full(i)/pow(i))
                 IF(error>tolerance) THEN
                    WRITE(*,*) 'HMx_DRIVER: Test failing'
                    WRITE(*,*) 'HMx_DRIVER: Test fields:', itest, jtest
                    WRITE(*,*) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                    WRITE(*,*) 'HMx_DRIVER: Redshift:', z
                    WRITE(*,*) 'HMx_DRIVER: Expected power:', pow(i)
                    WRITE(*,*) 'HMx_DRIVER: Model power:', pow_full(i)
                    WRITE(*,*) 'HMx_DRIVER: Tolerance:', tolerance
                    WRITE(*,*) 'HMx_DRIVER: Error:', error
                    WRITE(*,*)
                    ifail=.TRUE.
                 END IF
              END DO

              ! Write test information to the screen
              WRITE(*,*) 'HMx_DRIVER: Test fields:', itest, jtest
              WRITE(*,*) 'HMx_DRIVER: Tolerance:', tolerance
              WRITE(*,*) 'HMx_DRIVER: Max error:', error_max
              WRITE(*,*)

              ! Write data to file
              CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

              ! Deallocate arrays ready for next field combination
              DEALLOCATE(k,pow,pow_lin,pow_2h,pow_1h,pow_full)

           END DO
        END DO

     END DO

     ! Write pass/fail information to the screen
     IF(ifail) THEN
        WRITE(*,*) 'HMx_DRIVER: Hydro tests failed'
     ELSE
        WRITE(*,*) 'HMx_DRIVER: Hydro tests passed'        
     END IF
     WRITE(*,*)

  ELSE IF(imode==45) THEN

     ! 45 - Comparison of Sheth-Tormen vs. Tinker mass function

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_lin(nk),pow_2h(nk),pow_1h(nk),pow_full(nk))

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Sets the redshift
     z=0.

     DO j=1,2

        IF(j==1) THEN
           ! Sheth-Tormen mass function and bias
           ihm=3  
           outfile='data/power_ShethTormen.dat'
        ELSEIF(j==2) THEN
           ihm=23 ! Tinker mass function and bias
           outfile='data/power_Tinker.dat'
        END IF
        
        !Initiliasation for the halomodel calcualtion
        CALL assign_halomod(ihm,hmod,verbose)
        CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Do the halo-model calculation
        ip(1)=-1
        ip(2)=-1
        CALL calculate_halomod(ip,k,nk,z,pow_lin,pow_2h,pow_1h,pow_full,hmod,cosm,verbose,response=.FALSE.)

        ! Write out data to disk
        CALL write_power(k,pow_lin,pow_2h,pow_1h,pow_full,nk,outfile,verbose)

     END DO

  ELSE IF(imode==46) THEN

     ! 46 - Mass function plots

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Sets the redshift
     z=0.

     !Initiliasation for the halomodel calcualtion
     ihm=3
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     nu_min=0.1
     nu_max=6.
     n=256

     OPEN(7,file='data/bnu_functions.dat')
     OPEN(8,file='data/gnu_functions.dat')
     OPEN(9,file='data/mass_functions.dat')
     DO i=1,n
        nu=progression(nu_min,nu_max,i,n)
        mass=M_nu(nu,hmod)
        mf=mass_function(mass,hmod,cosm)
        WRITE(7,*) nu, mass, b_ps(nu,hmod), b_st(nu,hmod), b_Tinker(nu,hmod)
        WRITE(8,*) nu, mass, g_ps(nu,hmod), g_st(nu,hmod), g_Tinker(nu,hmod)
        WRITE(9,*) nu, mass, mf, mf*mass/comoving_matter_density(cosm), mf*mass**2/comoving_matter_density(cosm)
     END DO
     CLOSE(7)
     CLOSE(8)
     CLOSE(9)

  ELSE IF(imode==49) THEN

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the range of redshifts
     z1=0.
     z2=5.
     nz=6
     ALLOCATE(hmods(nz),HI_frac(nz),z_tab(nz))
     
     ! Initiliasation for the halomodel calcualtion at all zs
     ihm=25
     DO j=1,nz
        z_tab(j)=progression(z1,z2,j,nz)
        CALL assign_halomod(ihm,hmods(j),verbose)
        CALL init_halomod(mmin,mmax,z_tab(j),hmods(j),cosm,verbose)
        CALL print_halomod(hmods(j),cosm,verbose)
     END DO

     ! Set the range of halo masses
     m1=mmin
     m2=mmax
     n=256

     ! Loop over mass and z and do the calculation
     OPEN(7,file='data/HI_mass_fraction.dat')
     DO i=1,n
        mass=exp(progression(log(m1),log(m2),i,n))
        DO j=1,nz
           HI_frac(j)=halo_HI_fraction(mass,hmods(j),cosm)
        END DO
        WRITE(7,*) mass, (HI_frac(j), j=1,nz)
     END DO
     CLOSE(7)
     CLOSE(8)
     
  ELSE
        
     STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

  END IF

CONTAINS

  REAL FUNCTION figure_of_merit(a,b,n)

    ! A figure-of-merit or cost function
    IMPLICIT NONE
    REAL, INTENT(IN) :: a(n), b(n)
    INTEGER, INTENT(IN) :: n

    !figure_of_merit=(SUM(a/b)-REAL(n))**2
    figure_of_merit=sqrt(SUM((a/b-1.)**2)/REAL(n))
    !figure_of_merit=SUM(log(a/b)**2)/REAL(n)
    
  END FUNCTION figure_of_merit

  SUBROUTINE fom_multiple(imode,ip,fom,p,plog,np,k,nk,z,nz,pow,pow_sim,hmod,cosm,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: imode ! Mode for the driver
    INTEGER, INTENT(IN) :: ip(2) ! Field types
    REAL, INTENT(OUT) :: fom ! Output figure of merit
    REAL, INTENT(IN) :: p(np) ! Array of varying parameters
    LOGICAL, INTENT(IN) :: plog(np) ! Array of which parameters are explored in log
    INTEGER, INTENT(IN) :: np ! Number of varying parameters
    REAL, INTENT(IN) :: k(nk,nz,n) ! Array of k values for comparison data
    INTEGER, INTENT(IN) :: nk ! Number of k values for comparison data
    REAL, INTENT(IN) :: z(nz) ! Array of z values for comparison data
    INTEGER, INTENT(IN) :: nz ! Number of z values
    REAL, INTENT(OUT) :: pow(nk,nz,n) ! Output array for power as a function of k, z, cosmology
    REAL, INTENT(IN) :: pow_sim(nk,nz,n) ! Comparison power as a function of k, z, cosmology 
    TYPE(halomod), INTENT(INOUT) :: hmod(n) ! Array of halo models for each comparison
    TYPE(cosmology), INTENT(INOUT) :: cosm(n) ! Array of cosmological models for each comparison
    INTEGER, INTENT(IN) :: n ! Number of cosmological models being compared
    INTEGER :: i, j   
    REAL :: pow_lin(nk,nz,n), pow_2h(nk,nz,n), pow_1h(nk,nz,n), pexp(np)
    
    ! Set this counting output variable to zero
    fom=0.

    ! Set this to zero too, for the banter
    pow=0.

    ! Loop over cosmologies
    DO i=1,n

       ! MCMC: Set parameters here

       ! Exponentiate those parameters that need it
       DO j=1,np
          IF(plog(j)) THEN
             pexp(j)=10**p(j)
          ELSE
             pexp(j)=p(j)
          END IF
       END DO
       
       If(imode==25 .OR. imode==38 .OR. imode==40 .OR. imode==41) THEN

          ! Set HMcode parameters to those being varied
          hmod(i)%Dv0=pexp(1)
          hmod(i)%Dv1=pexp(2)
          hmod(i)%dc0=pexp(3)
          hmod(i)%dc1=pexp(4)
          hmod(i)%eta0=pexp(5)
          hmod(i)%eta1=pexp(6)
          hmod(i)%f0=pexp(7)
          hmod(i)%f1=pexp(8)
          hmod(i)%ks=pexp(9)
          hmod(i)%A=pexp(10)
          hmod(i)%alp0=pexp(11)
          hmod(i)%alp1=pexp(12)

       ELSE IF(imode==39) THEN

          ! Set hydro parameters to those being varied

          ! stars-stars complex
          !hmod(i)%Astar=pexp(1)
          !hmod(i)%Mstar=pexp(2)          
          !hmod(i)%cstar=pexp(3)
          !hmod(i)%sstar=pexp(4)

          ! stars-stars simple
          !hmod(i)%Astar=pexp(1)
          !hmod(i)%cstar=pexp(2)

          ! matter-pressure
          !hmod(i)%Twhim=10**pexp(1)

          ! gas-gas
          hmod(i)%eps=pexp(1)
          hmod(i)%Gamma=pexp(2)
          hmod(i)%M0=pexp(3)
          
       ELSE

          STOP 'FOM_MULTIPLE: Error, imode not specified correctly'
          
       END IF

       ! Loop over redshifts
       DO j=1,nz

          ! Initialise the halo-model calculation
          CALL init_halomod(mmin,mmax,z(j),hmod(i),cosm(i),verbose=.FALSE.)
          CALL print_halomod(hmod(i),cosm(i),verbose=.FALSE.)

          ! Calculate the halo-model power spectrum
          CALL calculate_halomod(ip,k(:,j,i),nk,z(j),pow_lin(:,j,i),pow_2h(:,j,i),pow_1h(:,j,i),pow(:,j,i),hmod(i),cosm(i),verbose=.FALSE.,response=.FALSE.)

          ! Calculate figure of merit and add to total
          fom=fom+figure_of_merit(pow(:,j,i),pow_sim(:,j,i),nk)**2

       END DO

    END DO

    ! Divide the figure-of-merit by the number of redshifts and cosmologies
    ! This is then the rms error per log-k, per z, per cosmology
    fom=sqrt(fom/REAL(nz*n))
    
  END SUBROUTINE fom_multiple

  SUBROUTINE set_mcmc_parameter_ranges(imode,ip,sigma,p,plog,np,k,nk,z,nz,pow_sim,hmod,cosm,n,verbose)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: imode
    INTEGER, INTENT(IN) :: ip(2)
    REAL, INTENT(OUT) :: sigma(n)
    REAL, INTENT(IN) :: p(np)
    LOGICAL, INTENT(IN) :: plog(np)
    INTEGER, INTENT(IN) :: np
    REAL, INTENT(IN) :: k(nk,nz,n)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: z(nz)
    INTEGER, INTENT(IN) :: nz
    REAL, INTENT(IN) :: pow_sim(nk,nz,n)
    TYPE(halomod), INTENT(INOUT) :: hmod(n)
    TYPE(cosmology), INTENT(INOUT) :: cosm(n)
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i
    REAL :: fom_base, fom, df, p2(np), pow(nk,nz,n)
    REAL :: delta, dp
    LOGICAL, PARAMETER :: check=.TRUE.

    ! MCMC: Set step-size stuff here

    IF(imode==25 .OR. imode==38 .OR. imode==40 .OR. imode==41) THEN
       ! Mira Titan or FrankenEmu
       delta=1e-4 ! Required change in figure-of-merit
       dp=1e-10 ! Used for derivative
    ELSE IF(imode==39) THEN
       ! Hydro
       delta=1e-2 ! Required change in figure-of-merit
       dp=1e-10 ! Used for derivative
    ELSE
       STOP 'SET_MCMC_PARAMETER_RANGES: Error imode specified incorrectly'
    END IF
       
    ! Get the figure of merit for the base set of parameters
    CALL fom_multiple(imode,ip,fom_base,p,plog,np,k,nk,z,nz,pow,pow_sim,hmod,cosm,n)

!!$    ! Calculate the change in parameter
!!$    DO i=1,np
!!$       IF(plog(i)) THEN
!!$          sigma(i)=p(i)+dp
!!$       ELSE
!!$          sigma(i)=p(i)*dp
!!$       END IF
!!$    END DO

    ! Is this the correct thing to do in log?
    sigma=p*dp

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Setting parameter jump sizes'
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Number of parameters:', np
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Number of cosmologies:', n
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Number of wavenumbers:', nk
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Number of redshifts:', nz
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Derivatives being calculated with:', dp
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Fixing sigma to give change in fom:', delta
       WRITE(*,*) '========================================================'
       WRITE(*,*) 'Parameter              Value         Sigma         Ratio'
       WRITE(*,*) '========================================================'
    END IF

    ! Loop over parameters
    DO i=1,np

       ! Set the range of p to take the derivative over
       p2=p ! Reset all
       p2(i)=p(i)+sigma(i) ! Perturb parameter i
       
       ! Get the figure of merit for the updated parameter
       CALL fom_multiple(imode,ip,fom,p2,plog,np,k,nk,z,nz,pow,pow_sim,hmod,cosm,n)

       ! Calculate the change in the figure of merit for this parameter
       df=fom-fom_base

       ! Check that perturbing the parameter actually changes the figure of merit
       IF(df==0.) THEN
          WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Parameter:', i
          WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: sigma:', sigma(i)
          IF(plog(i)) THEN
             WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Original value:', 10**p(i)
             WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Perturbed value:', 10**p2(i)
          ELSE
             WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Original value:', p(i)
             WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Perturbed value:', p2(i)
          END IF
          WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Original figure-of-merit:', fom_base
          WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Perturbed figure-of-merit:', fom
          WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Change in figure-of-merit:', df
          STOP 'SET_MCMC_PARAMETER_RANGES: Error, changing parameter does not change power spectra'
       END IF

       ! Set sigma so that it gives a change of 'delta' in fom
       sigma(i)=ABS(sigma(i)/df)*delta

       ! Write parameters to screen
       IF(verbose) THEN
          IF(plog(i)) WRITE(*,fmt='(I10,A5,3F14.7)') i, 'lin', 10**p(i), sigma(i), sigma(i)/ABS(p(i))
          WRITE(*,fmt='(I10,A5,3F14.7)') i, 'log', p(i), sigma(i), sigma(i)/ABS(p(i))
       END IF

       ! Do an extra check if necessary
       IF(check) THEN
          p2(i)=p(i)+sigma(i)
          CALL fom_multiple(imode,ip,fom,p2,plog,np,k,nk,z,nz,pow,pow_sim,hmod,cosm,n)
          WRITE(*,*) i, fom_base, fom, fom-fom_base
       END IF
     
    END DO

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) '========================================================'
       WRITE(*,*) 'SET_MCMC_PARAMETER_RANGES: Done'
       WRITE(*,*)
    END IF
    
  END SUBROUTINE set_mcmc_parameter_ranges

  SUBROUTINE read_simulation_power_spectrum(k,Pk,n,infile,cut_nyquist,subtract_shot,verbose)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), Pk(:) ! Output simulation k and power
    INTEGER, INTENT(OUT) :: n ! Number of output k values
    CHARACTER(len=*), INTENT(IN) :: infile ! Input file location
    LOGICAL, INTENT(IN) :: cut_nyquist ! Logical to cut Nyquist or not
    LOGICAL, INTENT(IN) :: subtract_shot ! Logical to subtract shot noise or not
    LOGICAL, INTENT(IN) :: verbose ! Logical verbose
    INTEGER :: i
    REAL :: shot
    LOGICAL :: lexist

    ! Deallocate arrays if they are already allocated
    IF(ALLOCATED(k))  DEALLOCATE(k)
    IF(ALLOCATED(Pk)) DEALLOCATE(Pk)

    ! Check file exists
    INQUIRE(file=infile,exist=lexist)
    IF(.NOT. lexist) THEN
       WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: File: ', TRIM(infile)
       STOP 'READ_SIMULATION_POWER_SPECTRUM: File does not exist'
    END IF

    ! Get file length and allocate arrays for output
    n=file_length(infile,verbose=.FALSE.)
    ALLOCATE(k(n),Pk(n))

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: Reading in data'
       WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: File: ', TRIM(infile)
       WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: nk:', nk
    END IF

    ! Read in data from file
    OPEN(7,file=infile,status='old')
    DO i=1,n
       READ(7,*) k(i), Pk(i), shot
       IF(subtract_shot) Pk(i)=Pk(i)-shot
    END DO
    CLOSE(7)  

    IF(cut_nyquist) THEN

       ! Find position in array of half-Nyquist
       kmax=k(n)
       DO i=1,n
          IF(k(i)>kmax/2.) EXIT
       END DO

       ! Cut arrays down to half-Nyquist
       CALL amputate(k,n,i)
       CALL amputate(Pk,n,i)
       n=i

       ! Write to screen
       IF(verbose) THEN
          WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: Trimmed to Nyquist frequency'
       END IF

    END IF

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'READ_SIMULATION_POWER_SPECTRUM: Done'
       WRITE(*,*)
    END IF
    
  END SUBROUTINE read_simulation_power_spectrum

  SUBROUTINE get_k_values(infile,k,nk)

    ! Get k values from a data file, assumes they are the first column of a data file
    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile ! Input file location
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:) ! Output array of k values
    INTEGER, INTENT(OUT) :: nk ! Output number of k values

    ! Get the number of k values
    nk=file_length(infile,verbose=.FALSE.)

    ! Allocate the array in k
    ALLOCATE(k(nk))

    ! Read in the k values
    OPEN(7,file=infile)
    DO i=1,nk
       READ(7,*) k(i)
    END DO
    CLOSE(7)
    
  END SUBROUTINE get_k_values

  SUBROUTINE YinZhe_Fig1(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmological model
    INTEGER :: i
    REAL :: r, rs, rv, c, Mh, rh, r500c, m500c

    REAL, PARAMETER :: M=1e15 !Halo virial? mass [Msun]
    REAL, PARAMETER :: rmin=1e-3!Minimum radius [Mpc]
    REAL, PARAMETER :: rmax=8 !Maximum radius [Mpc] 
    INTEGER, PARAMETER :: nr=512 !Number of points in radius

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod,cosm)

    Mh=M*cosm%h !This virial mass is now [Msun/h]

    rv=exp(find(log(Mh),hmod%log_m,log(hmod%rv),hmod%n,3,3,2)) ![Mpc/h]
    c=find(log(Mh),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c ![Mpc/h]

    m500c=exp(find(log(Mh),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ![Mpc/h]
    r500c=exp(find(log(Mh),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ![Mpc/h]

    WRITE(*,*) 'YINZHE_FIG1: Making data for this figure'
    WRITE(*,*) 'YINZHE_FIG1: Redshift:', z
    WRITE(*,*) 'YINZHE_FIG1: Virial radius [Mpc]:', rv/cosm%h
    WRITE(*,*) 'YINZHE_FIG1: Virial radius [Mpc/h]:', rv
    WRITE(*,*) 'YINZHE_FIG1: r_500,c [Mpc]:', r500c/cosm%h
    WRITE(*,*) 'YINZHE_FIG1: r_500,c [Mpc/h]:', r500c
    WRITE(*,*) 'YINZHE_FIG1: r_500,c / r_v:', r500c/rv
    WRITE(*,*) 'YINZHE_FIG1: Virial halo mass [log10 Msun]:', log10(M)
    WRITE(*,*) 'YINZHE_FIG1: Virial halo mass [log10 Msun/h]:', log10(Mh)    
    WRITE(*,*) 'YINZHE_FIG1: M_500,c [log10 Msun]:', log10(M500c/cosm%h)
    WRITE(*,*) 'YINZHE_FIG1: M_500,c [log10 Msun/h]:', log10(M500c)
    WRITE(*,*) 'YINZHE_FIG1: M_500,c / M_v:', M500c/Mh
    WRITE(*,*) 'YINZHE_FIG1: Halo concentraiton:', c

    OPEN(7,file='data/YinZhe_Fig1.dat')
    DO i=1,nr
       r=progression(rmin,rmax,i,nr) !Radius [Mpc]
       rh=r*cosm%h !Convert [Mpc/h]
       WRITE(7,*) r, UPP(.TRUE.,rh,Mh,rv,rs,hmod,cosm)*r**2, win_electron_pressure(.TRUE.,1,rh,Mh,rv,rs,hmod,cosm)*r**2
    END DO
    CLOSE(7)

    WRITE(*,*) 'YINZHE_FIG1: Done'
    WRITE(*,*)
    
  END SUBROUTINE YinZhe_Fig1

   SUBROUTINE write_power(k,pow_lin,pow_2h,pow_1h,pow,nk,output,verbose)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    IF(verbose) WRITE(*,*) 'WRITE_POWER: Writing power to ', TRIM(output)

    ! Loop over k values
    ! Fill the tables with one- and two-halo terms as well as total
    OPEN(7,file=output)
    DO i=1,nk       
       WRITE(7,fmt='(5ES20.10)') k(i), pow_lin(i), pow_2h(i), pow_1h(i), pow(i)
    END DO
    CLOSE(7)

    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_power

  SUBROUTINE write_power_a_multiple(k,a,pow_lin,pow_2h,pow_1h,pow_full,nk,na,base,verbose)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: base
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow_lin(nk,na), pow_2h(nk,na), pow_1h(nk,na), pow_full(nk,na)
    LOGICAL, INTENT(IN) :: verbose
    REAL :: pow(nk,na)
    INTEGER :: i
    CHARACTER(len=512) :: output
    LOGICAL :: verbose2

    DO i=1,4
       IF(i==1) THEN
          output=TRIM(base)//'_linear.dat'
          pow=pow_lin
       ELSE IF(i==2) THEN
          output=TRIM(base)//'_2halo.dat'
          pow=pow_2h
       ELSE IF(i==3) THEN
          output=TRIM(base)//'_1halo.dat'
          pow=pow_1h
       ELSE IF(i==4) THEN
          output=TRIM(base)//'_full.dat'
          pow=pow_full
       ELSE
          STOP 'WRITE_POWER_A_MULTIPLE: Error, something went FUBAR'
       END IF
       IF(i==1) THEN
          verbose2=verbose
       ELSE
          verbose2=.FALSE.
       END IF
       CALL write_power_a(k,a,pow,nk,na,output,verbose2)
    END DO

  END SUBROUTINE write_power_a_multiple

  SUBROUTINE write_power_a(k,a,pow,nk,na,output,verbose)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: output
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i, j

    ! Print to screen
    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER_A: The first entry of the file is hashes - #####'
       WRITE(*,*) 'WRITE_POWER_A: The remainder of the first row are the scale factors - a'
       WRITE(*,*) 'WRITE_POWER_A: The remainder of the first column are the wave numbers - k'
       WRITE(*,*) 'WRITE_POWER_A: Each row then gives the power at that k and a'
       WRITE(*,*) 'WRITE_POWER_A: Output:', TRIM(output)
    END IF

    ! Write out data to files
    OPEN(7,file=output)
    DO i=0,nk
       IF(i==0) THEN
          WRITE(7,fmt='(A20,40F20.10)') '#####', (a(j), j=1,na)
       ELSE
          WRITE(7,fmt='(F20.10,40E20.10)') k(i), (pow(i,j), j=1,na)
       END IF
    END DO
    CLOSE(7)

    ! Print to screen
    IF(verbose) THEN
       WRITE(*,*) 'WRITE_POWER_A: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_power_a

!!$  SUBROUTINE write_distances(cosm)
!!$
!!$    ! Write file of z vs. r(z)
!!$    IMPLICIT NONE
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    CHARACTER(len=256) :: output
!!$    INTEGER :: i
!!$    REAL :: z
!!$
!!$    ! Now write the results of r(z) calculation
!!$    output='data/distance.dat'
!!$    WRITE(*,*) 'WRITE_DISTANCE: Writing r(a): ', TRIM(output)
!!$    OPEN(7,file=output)
!!$    DO i=1,cosm%n_r
!!$       z=redshift_a(cosm%a_r(i))
!!$       WRITE(7,*) z, cosm%r(i), f_k(cosm%r(i),cosm)
!!$    END DO
!!$    CLOSE(7)
!!$    WRITE(*,*) 'WRITE_DISTANCE: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE write_distances

  SUBROUTINE random_baryon_parameters(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod

    REAL, PARAMETER :: alpha_min=0.05
    REAL, PARAMETER :: alpha_max=2.5

    REAL, PARAMETER :: eps_min=0.5
    REAL, PARAMETER :: eps_max=2.0

    REAL, PARAMETER :: Gamma_min=1.05
    REAL, PARAMETER :: Gamma_max=2.00

    REAL, PARAMETER :: M0_min=10**(12.)
    REAL, PARAMETER :: M0_max=10**(15.)

    REAL, PARAMETER :: Astar_min=0.002
    REAL, PARAMETER :: Astar_max=0.2

    REAL, PARAMETER :: Twhim_min=10**(5.)
    REAL, PARAMETER :: Twhim_max=10**(7.)
    
    hmod%alpha=random_uniform(alpha_min,alpha_max)

    hmod%eps=random_uniform(log(eps_min),log(eps_max))
    hmod%eps=exp(hmod%eps)

    hmod%Gamma=random_uniform(Gamma_min,Gamma_max)

    hmod%M0=random_uniform(log(M0_min),log(M0_max))
    hmod%M0=exp(hmod%M0)
    
    hmod%Astar=random_uniform(Astar_min,Astar_max)

    hmod%Twhim=random_uniform(log(Twhim_min),log(Twhim_max))
    hmod%Twhim=exp(hmod%Twhim)
    
  END SUBROUTINE random_baryon_parameters
  
  !SUBROUTINE BAHAMAS_power_file_name(name,model,z,ip)
  CHARACTER(len=256) FUNCTION BAHAMAS_power_file_name(model,z,ip)

    IMPLICIT NONE
    !CHARACTER(len=*), INTENT(OUT) :: name
    CHARACTER(len=*), INTENT(IN) :: model
    REAL, INTENT(IN) :: z
    INTEGER, INTENT(IN) :: ip(2)
    CHARACTER(len=64) :: dir
    CHARACTER(len=32) :: snap, field(2)
    LOGICAL :: lexist
    INTEGER :: j
    
    ! Directory containing everything
    dir='/Users/Mead/Physics/BAHAMAS/power/M1024'

    ! Set the redshift
    IF(z==0.0) THEN
       snap='snap32'
    ELSE IF(z==0.5) THEN
       snap='snap28'
    ELSE IF(z==1.0) THEN
       snap='snap26'
    ELSE IF(z==2.0) THEN
       snap='snap22'
    ELSE
       STOP 'BAHAMAS_POWER_FILE_NAME: Error, redshift specified incorrectly'
    END IF

    ! Set the fields
    DO j=1,2
       IF(ip(j)==0) THEN
          field(j)='all'
       ELSE IF(ip(j)==1) THEN
          field(j)='dm'
       ELSE IF(ip(j)==2) THEN
          field(j)='gas'
       ELSE IF(ip(j)==3) THEN
          field(j)='stars'
       ELSE IF(ip(j)==6) THEN
          field(j)='epressure'
       ELSE
          WRITE(*,*) 'BAHAMAS_POWER_FILE_NAME: Field number', j
          WRITE(*,*) 'BAHAMAS_POWER_FILE_NAME: Halo type', ip(j)
          STOP 'BAHAMAS_POWER_FILE_NAME: Error, ip specified incorrectly'
       END IF
    END DO

    ! File name
    BAHAMAS_power_file_name=TRIM(dir)//'/'//TRIM(model)//'_L400N1024_WMAP9_'//TRIM(snap)//'_'//TRIM(field(1))//'_'//TRIM(field(2))//'_power.dat'

    INQUIRE(file=BAHAMAS_power_file_name,exist=lexist)
    IF(.NOT. lexist) THEN
       WRITE(*,*) 'BAHAMAS_POWER_FILE_NAME: ', TRIM(BAHAMAS_power_file_name)
       STOP 'BAHAMAS_POWER_FILE_NAME: Error, file does not exist'
    END IF
    
  END FUNCTION BAHAMAS_power_file_name

  SUBROUTINE get_BAHAMAS_power(k,Pk,nk,z,name,field,cosm)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
    INTEGER, INTENT(OUT) :: nk
    REAL, INTENT(IN) :: z
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: field(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, ALLOCATABLE :: Pk_DM(:), Pk_HMcode(:)
    CHARACTER(len=256) :: infile, dmonly
    
    INTEGER, PARAMETER :: field_all_matter(2)=0
    REAL, PARAMETER :: mmin=1e7
    REAL, PARAMETER :: mmax=1e17

    dmonly=BAHAMAS_power_file_name('DMONLY_2fluid_nu0',z,field_all_matter)
    infile=BAHAMAS_power_file_name(name,z,field)

    CALL read_simulation_power_spectrum(k,Pk_DM,nk,dmonly,cut_nyquist=.FALSE.,subtract_shot=.TRUE.,verbose=.TRUE.)
    CALL read_simulation_power_spectrum(k,Pk,nk,infile,cut_nyquist=.FALSE.,subtract_shot=.TRUE.,verbose=.TRUE.)
    Pk=Pk/Pk_DM
    
    ALLOCATE(Pk_HMcode(nk))
    CALL calculate_HMcode(k,z,Pk_HMcode,nk,cosm)
    Pk=Pk*Pk_HMcode
    
  END SUBROUTINE get_BAHAMAS_power

END PROGRAM HMx_driver
