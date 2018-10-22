PROGRAM HMx_driver

  USE HMx
  USE Limber
  USE file_info
  USE random_numbers
  USE cosmic_emu_stuff
  USE calculus_table
  USE string_operations

  ! TODO: Many pow(1,1,nk) could be pow(nk)
  ! TODO: Too many different pow(:,:), pow(:,:,:), pow(:,:,:,:), ...
  ! TODO: Maybe could use pow1(:), pow2(:,:), ...
  
  IMPLICIT NONE

  ! Parameter definitions
  REAL, ALLOCATABLE :: k(:), a(:)
  REAL, ALLOCATABLE :: k_sim(:), pow_sim(:)
  REAL, ALLOCATABLE :: pow_li(:), pow_2h(:,:,:), pow_1h(:,:,:), pow_hm(:,:,:)
  REAL, ALLOCATABLE :: pow_ka(:,:)
  REAL, ALLOCATABLE :: powd_li(:), powd_2h(:), powd_1h(:), powd_hm(:)
  REAL, ALLOCATABLE :: pows_li(:,:), pows_2h(:,:,:,:), pows_1h(:,:,:,:), pows_hm(:,:,:,:)
  REAL, ALLOCATABLE :: powb_hm(:,:,:,:)
  REAL, ALLOCATABLE :: ell(:), Cl(:), theta(:), xi(:,:), zs(:), masses(:)
  REAL, ALLOCATABLE :: z_tab(:), HI_frac(:)
  INTEGER :: i, j, ii, nk, na, j1, j2, itriad, nf
  INTEGER :: n, nl, nz, nth, nnz, m, ipa, npa, ncos, ncore, ntriad
  INTEGER :: ip(2), ix(2), ixx(2), field(1)
  INTEGER, ALLOCATABLE :: fields(:)
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
  CHARACTER(len=256) :: mode, halomodel, cosmo
  INTEGER :: imode, icosmo, iowl, ihm, irho, itest, jtest
  REAL :: sig8min, sig8max
  REAL :: mass, m1, m2, nu, nu_min, nu_max, mf
  REAL :: c, rmin, rmax, rv, rs, p1, p2
  REAL :: spam
  CHARACTER(len=1) :: crap
  LOGICAL :: verbose2
  
  ! Baryon stuff
  REAL :: param_min, param_max, param, param_neat
  LOGICAL :: ilog

  ! Halo-model Parameters
  Logical, PARAMETER :: verbose=.TRUE. ! Verbosity
  REAL, PARAMETER :: mmin=1e7        ! Minimum halo mass for the calculation
  REAL, PARAMETER :: mmax=1e17       ! Maximum halo mass for the calculation

  ! Test parameters
  REAL, PARAMETER :: tolerance=3e-3
  REAL :: error, error_max
  LOGICAL :: verbose_tests=.FALSE.
  LOGICAl :: ifail=.FALSE.

  ! Benchmark parameters
  LOGICAL, PARAMETER :: Alonso_k=.TRUE.

  ! Output choices
  LOGICAL, PARAMETER :: icumulative=.TRUE. ! Do cumlative distributions for breakdown
  LOGICAL, PARAMETER :: ifull=.FALSE.      ! Do only full halo model C(l), xi(theta) calculations (quicker, no breakdown ...)

  ! Cross correlation
  REAL, PARAMETER :: kmin_xcorr=1e-3 ! Minimum k
  REAL, PARAMETER :: kmax_xcorr=1e1  ! Maximum k (some halo-model things go hectic for k>10)
  INTEGER, PARAMETER :: nk_xcorr=64  ! Number of k values (used to be 32)
  REAL, PARAMETER :: amin_xcorr=0.1 ! Minimum scale factor (problems with one-halo term if amin is less than 0.1 (CMB lensing?))
  REAL, PARAMETER :: amax_xcorr=1.0 ! Maximum scale factor
  INTEGER, PARAMETER :: na_xcorr=16 ! Number of scale factores

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
     WRITE(*,*) ' 5 - Lensing diagnostics'
     WRITE(*,*) ' 6 - n(z) check'
     WRITE(*,*) ' 7 - Do general angular cross correlation'
     WRITE(*,*) ' 8 - Angular cross correlation as a function of cosmology'
     WRITE(*,*) ' 9 - Breakdown angular correlations in halo mass'
     WRITE(*,*) '10 - Breakdown angular correlations in redshift'
     WRITE(*,*) '11 - Do general angular cross correlation and correlation functions'
     WRITE(*,*) '12 - Project triad'
     WRITE(*,*) '13 - Cross-correlation coefficient'
     WRITE(*,*) '14 - 3D spectra for variations in baryon parameters'
     WRITE(*,*) '15 - 3D spectra for cosmo-OWLS models (k values) NOT SUPPORTED'
     WRITE(*,*) '16 - 3D spectra for BAHAMAS models (k values)'
     WRITE(*,*) '17 - 3D spectra for user choice of fields'
     WRITE(*,*) '18 - 3D bias'
     WRITE(*,*) '19 - Make CCL benchmark data'
     WRITE(*,*) '20 - Make data for Ma et al. (2015) Fig. 1'
     WRITE(*,*) '21 - W(k) integrand diagnostics'
     WRITE(*,*) '22 - Time W(k) integration methods'
     WRITE(*,*) '23 - Produce DE response results from Mead (2017)'
     WRITE(*,*) '24 - '
     WRITE(*,*) '25 - '
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
     WRITE(*,*) '38 - '
     WRITE(*,*) '39 - '
     WRITE(*,*) '40 - '
     WRITE(*,*) '41 - '
     WRITE(*,*) '42 - PAPER: Contributions to k-k C(l) integral'
     WRITE(*,*) '43 - PAPER: Contributions to k-y C(l) integral'
     WRITE(*,*) '44 - PAPER: Project triad'
     WRITE(*,*) '45 - Comparison of Sheth-Tormen vs. Tinker mass function'
     WRITE(*,*) '46 - Mass function and bias plots'
     WRITE(*,*) '47 - Make CMB lensing to compare with CAMB'
     WRITE(*,*) '48 - HI bias'
     WRITE(*,*) '49 - HI mass fractions'
     WRITE(*,*) '50 - Mass function changes with Lbox'
     WRITE(*,*) '51 - Compare power with and without scatter'
     READ(*,*) imode
     WRITE(*,*) '============================'
     WRITE(*,*)
  END IF
          
  IF(imode==0) THEN

     ! Calculate halo model at one z
     ! TODO: Change to calculate_HMx_a

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Sets the redshift
     z=0.

     !Initiliasation for the halomodel calcualtion
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     ! Allocate arrays
     ALLOCATE(powd_li(nk),powd_2h(nk),powd_1h(nk),powd_hm(nk))

     ! Do the halo-model calculation
     field=field_dmonly
     CALL calculate_HMx_a(field,1,k,nk,powd_li,powd_2h,powd_1h,powd_hm,hmod,cosm,verbose,response=.FALSE.)

     ! Write out the results
     outfile='data/power.dat'
     CALL write_power(k,powd_li,powd_2h,powd_1h,powd_hm,nk,outfile,verbose)

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
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Set the number of redshifts and range (linearly spaced) and convert z -> a
     na=16
     amin=0.2
     amax=1.0
     CALL fill_array(amin,amax,a,na)
     !a=1./(1.+a)
     !na=nz

     field=field_dmonly
     CALL calculate_HMx(field,1,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)
     
     base='data/power'
     CALL write_power_a_multiple(k,a,pows_li,pows_2h(1,1,:,:),pows_1h(1,1,:,:),pows_hm(1,1,:,:),nk,na,base,verbose)

  ELSE IF(imode==2 .OR. imode==15 .OR. imode==16 .OR. imode==32) THEN

     ! Make cross power spectra of all different components of haloes as well as pressure

     IF(imode==2 .OR. imode==32) THEN

        ! Generic hydro
        
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
        
     ELSE IF(imode==15) THEN

        ! cosmo-OWLS

        STOP 'HMx_DRIVER: not tested in ages, be very careful'

        ! Do from REF, NOCOOL, AGN, AGN 8.5, AGN 8.7
        n=5

        ! Set the redshift
        nz=1
        ALLOCATE(z_tab(nz))
        z_tab(1)=0.

        ! Get the k values from the simulation measured P(k)
        infile='/Users/Mead/Physics/cosmo-OWLS/power/N800/DMONLY_all_all_power.dat'
        CALL read_k_values(infile,k,nk)

     ELSE IF(imode==16) THEN

        ! BAHAMAS
       
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
        CALL read_k_values(infile,k,nk)
        
     END IF

     ! Field types
     nf=5
     ALLOCATE(fields(nf))
     fields(1)=field_matter
     fields(2)=field_cdm
     fields(3)=field_gas
     fields(4)=field_star
     fields(5)=field_electron_pressure

     ! Allocate the arrays for P(k)
     ALLOCATE(powd_li(nk),powd_2h(nk),powd_1h(nk),powd_hm(nk))
     ALLOCATE(pow_li(nk),pow_2h(nf,nf,nk),pow_1h(nf,nf,nk),pow_hm(nf,nf,nk))

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
            
           ! Initiliasation for the halomodel calcualtion after variables changed
           CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
           CALL print_halomod(hmod,cosm,verbose)

           ! Runs the diagnostics
           IF(imode==2) THEN
              dir='data'
              CALL halo_diagnostics(hmod,cosm,dir)
              CALL halo_definitions(hmod,cosm,dir)
              CALL halo_properties(hmod,dir)
           END IF

           IF(imode==2 .OR. imode==32) THEN
              ! File base and extension
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

           ! Dark-matter only
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

           ! Write some things to the screen
           field=field_dmonly
           WRITE(*,*) field(1), field(1), TRIM(outfile)

           ! Do the calculation for DMonly
           CALL calculate_HMx_a(field,1,k,nk,powd_li,powd_2h,powd_1h,powd_hm,hmod,cosm,verbose,response=.FALSE.)
           CALL write_power(k,powd_li,powd_2h,powd_1h,powd_hm,nk,outfile,verbose)

           ! Do the calculation for the rest of the fields
           CALL calculate_HMx_a(fields,nf,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response=.FALSE.)

           ! Loop over fields and write data to disk
           DO j1=1,nf
              DO j2=j1,nf

                 ! Fix output file and write to screen
                 outfile=number_file2(base,fields(j1),mid,fields(j2),ext)

                 ! Set the halo types and write to screen
                 WRITE(*,*) fields(j1), fields(j2), TRIM(outfile)

                 ! Write P(k) to disk
                 CALL write_power(k,pow_li,pow_2h(j1,j2,:),pow_1h(j1,j2,:),pow_hm(j1,j2,:),nk,outfile,verbose=.FALSE.)

              END DO
           END DO

        END DO

     END DO

  ELSE IF(imode==3) THEN

     ! Diagnostics
     
     ! Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Loop over redshifts
     DO j=1,4

        ! Set the redshift
        IF(j==1) z=0.0
        IF(j==2) z=0.5
        IF(j==3) z=1.0
        IF(j==4) z=2.0

        ! Initiliasation for the halomodel calcualtion
        CALL assign_halomod(ihm,hmod,verbose)
        CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Runs the diagnostics
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

     ! Set the number of fields
     nf=5
     ALLOCATE(fields(nf))
     fields(1)=field_matter
     fields(2)=field_cdm
     fields(3)=field_gas
     fields(4)=field_star
     fields(5)=field_electron_pressure

     ! Allocate arrays for power
     ALLOCATE(pow_li(nk),pow_2h(nf,nf,nk),pow_1h(nf,nf,nk),pow_hm(nf,nf,nk))

     ! Assigns the cosmological model
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Assign the default halo model
     CALL assign_halomod(ihm,hmod,verbose)

     ! Set the redshift
     z=0.

     ! Loop forever
     DO

        ! Initiliasation for the halomodel calcualtion        
        CALL random_baryon_parameters(hmod)
        CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        CALL calculate_HMx_a(fields,nf,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

        ! Do the halo-model calculation for a range of halo types
        DO j1=1,nf
           DO j2=1,nf
              
              DO i=1,nk
                 IF(ISNAN(pow_hm(j1,j2,i))) THEN
                    CALL print_halomod(hmod,cosm,verbose)
                    WRITE(*,*) 'HMx_DRIVER: Halo types:', fields(j1), fields(j2)
                    WRITE(*,*) 'HMx_DRIVER: k [h/Mpc]:', k(i)
                    STOP 'HMX_DRIVER: Error, NaN found in pow_full array'
                 END IF                 
              END DO
              
           END DO
        END DO

        WRITE(*,*)

     END DO

  ELSE IF(imode==5) THEN

     ! Lensing diagnostics

     STOP 'HMx_DRIVER: Error, you need to actually code this up'

     CALL write_nz(lens,outfile)

     CALL write_lensing_efficiency(lens,cosm,outfile)

     DO i=1,2
        CALL write_projection_kernel(proj(i),cosm,outfile)
     END DO

  ElSE IF(imode==6) THEN

     ! n(z) normalisation check

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
        CALL read_nz(nz,lens)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (linear):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,1)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (quadratic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,1,lens%nnz,2)
        WRITE(*,*) 'HMx_DRIVER: n(z) integral (cubic):', integrate_table(lens%z_nz,lens%nz,lens%nnz,2,lens%nnz,3)
        WRITE(*,*)
     END DO

  ELSE IF(imode==7 .OR. imode==8 .OR. imode==9 .OR. imode==10 .OR. imode==11 .OR. imode==37 .OR. imode==42 .OR. imode==43 .OR. imode==47) THEN

     ! General stuff for various 2D projections
     !  7 - Do general angular cross correlation
     !  8 - Angular cross correlation as a function of cosmology
     !  9 - Breakdown angular correlations in halo mass
     ! 10 - Breakdown angular correlations in redshift
     ! 11 - Breakdown angular correlations in halo radius
     ! 37 - CFHTLenS angular correlations
     ! 42 - PAPER: breakdown lensing-lensing per ell
     ! 43 - PAPER: breakdown lensing-y per ell
     ! 47 - Make CMB-lensing data to compare with CAMB
     
     ! Set the fields
     ix=-1
     IF(imode==37) ix=tracer_CFHTLenS    ! CFHTLenS autospectrum
     IF(imode==42) ix=tracer_KiDS_450    ! KiDS-450 autospectrum
     IF(imode==47) ix=tracer_CMB_lensing ! CMB lensing autospectrum
     IF(imode==43) THEN
        ix(1)=tracer_KiDS_450  ! KiDS-450 z = 0.1->0.9
        ix(2)=tracer_Compton_y ! Compton y
     END IF
     CALL set_xcorr_type(ix,ip)
     IF(imode==37 .OR. imode==47) ip=field_dmonly ! Set DMONLY haloes

     ! Assign the cosmological model
     IF(imode==37) icosmo=4 ! WMAP9 fits CFHTLenS okay
     IF(imode==42 .OR. imode==43) icosmo=1 ! Boring
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

     ! Set the ell range
     lmin=1
     lmax=1e4 ! Problems if this is pushed up to 10^5
     IF(imode==11 .OR. imode==37) lmax=1e5 ! Need higher lmax for the correlation functions
     nl=128

     ! Allocate arrays for l and C(l)
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cl(nl))

     ! Set the angular arrays in degrees and allocate arrays for theta and xi(theta)
     thmin=0.01
     thmax=10.
     nth=128
     CALL fill_array(log(thmin),log(thmax),theta,nth)
     theta=exp(theta)
     ALLOCATE(xi(3,nth))

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

     IF(imode==7 .OR. imode==11 .OR. imode==37 .OR. imode==42 .OR. imode==43 .OR. imode==47) THEN

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        ! Write out diagnostics
        IF(imode==37 .OR. imode==47) ihm=1  ! HMcode (2016)
        IF(imode==42 .OR. imode==43) ihm=20 ! Standard halo-model in response
        CALL assign_halomod(ihm,hmod,verbose)
        CALL calculate_HMx(ip,2,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)

        base=TRIM(dir)//'power'
        CALL write_power_a_multiple(k,a,pows_li,pows_2h(1,2,:,:),pows_1h(1,2,:,:),pows_hm(1,2,:,:),nk,na,base,verbose)

        ! Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        !CALL write_projection_kernels(proj,cosm)

        ! Set the distance range for the Limber integral
        !r1=100.
        !r2=proj%rs
        r1=0.
        r2=maxdist(proj)

        ! Write information to screen
        WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
        WRITE(*,*) 'HMx_DRIVER: r min [Mpc/h]:', r1
        WRITE(*,*) 'HMx_DRIVER: r max [Mpc/h]:', r2
        WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
        WRITE(*,*) 'HMx_DRIVER: lower limit of Limber integral [Mpc/h]:', REAL(r1)
        WRITE(*,*) 'HMx_DRIVER: upper limit of Limber integral [Mpc/h]:', REAL(r2)
        WRITE(*,*)

        ALLOCATE(pow_ka(nk,na))

        ! Loop over all types of C(l) to create
        DO j=4,1,-1

           IF(ifull .AND. (j .NE. 4)) CYCLE

           ! Write information to screen
           IF(j==1) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing linear'
              outfile=TRIM(dir)//'cl_linear.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_linear.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_linear.dat'
              pow_ka=pows_li
           ELSE IF(j==2) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing 2-halo'
              outfile=TRIM(dir)//'cl_2h.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_2h.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_2h.dat'
              pow_ka=pows_2h(1,2,:,:)
           ELSE IF(j==3) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing 1-halo'
              outfile=TRIM(dir)//'cl_1h.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_1h.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_1h.dat'
              pow_ka=pows_1h(1,2,:,:)
           ELSE IF(j==4) THEN
              WRITE(*,*) 'HMx_DRIVER: Doing full'
              outfile=TRIM(dir)//'cl_hm.dat'
              IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_cl_hm.dat'
              IF(imode==47) outfile=TRIM(dir)//'CMBlensing_cl_hm.dat'
              pow_ka=pows_hm(1,2,:,:)
           ELSE
              STOP 'HMx_DRIVER: Something went wrong'
           END IF

           WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

           ! Actually calculate the C(ell)
           CALL calculate_Cl(r1,r2,ell,Cl,nl,k,a,pow_ka,nk,na,proj,cosm)
           CALL write_Cl(ell,Cl,nl,outfile)

           IF(j==4 .AND. (imode==7 .OR. imode==11 .OR. imode==42 .OR. imode==43)) THEN
              CALL Cl_contribution_ell(r1,r2,k,a,pow_ka,nk,na,proj,cosm)
              IF(imode==42 .OR. imode==43) EXIT
           END IF

           IF(imode==11 .OR. imode==37) THEN

              ! Set xi output files
              IF(j==1) THEN
                 outfile=TRIM(dir)//'xi_linear.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_linear.dat'
              ELSE IF(j==2) THEN
                 outfile=TRIM(dir)//'xi_2h.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_2h.dat'                 
              ELSE IF(j==3) THEN
                 outfile=TRIM(dir)//'xi_1h.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_1h.dat'
              ELSE IF(j==4) THEN
                 outfile=TRIM(dir)//'xi_hm.dat'
                 IF(imode==37) outfile=TRIM(dir)//'CFHTLenS_xi_hm.dat'
              ELSE
                 STOP 'HMx_DRIVER: Something went wrong'
              END IF
              
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)

              ! Actually calculate the xi(theta)
              CALL calculate_xi(theta,xi,nth,ell,Cl,nl,NINT(lmax))
              CALL write_xi(theta,xi,nth,outfile)

           END IF

        END DO
        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     ELSE IF(imode==8) THEN

        ! Assess cross-correlation as a function of cosmology

        ! Allocate array for power
        ALLOCATE(pow_ka(nk,na))

        ! Set range in sigma_8
        sig8min=0.7
        sig8max=0.9
        ncos=5 ! I may have changed this number inadvertantly

        ! Loop over cosmology
        DO i=1,ncos

           cosm%sig8=progression(sig8min,sig8max,i,ncos)
           CALL init_cosmology(cosm)
           CALL print_cosmology(cosm)

           CALL assign_halomod(ihm,hmod,verbose)
           CALL calculate_HMx(ip,2,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)

           ! Fill out the projection kernels
           CALL fill_projection_kernels(ix,proj,cosm)
           !CALL write_projection_kernels(proj,cosm)

           ! Write to screen
           WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
           WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
           WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
           WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
           WRITE(*,*)

           ! Loop over all types of C(l) to create
           dir='data/'
           base=TRIM(dir)//'cosmology_'
           DO j=1,4
              
              IF(j==1) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) linear'
                 ext='_cl_linear.dat'
                 pow_ka=pows_li
              ELSE IF(j==2) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) 2-halo'
                 ext='_cl_2h.dat'
                 pow_ka=pows_2h(1,2,:,:)
              ELSE IF(j==3) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) 1-halo'
                 ext='_cl_1h.dat'
                 pow_ka=pows_1h(1,2,:,:)
              ELSE IF(j==4) THEN
                 WRITE(*,*) 'HMx_DRIVER: Doing C(l) full'
                 ext='_cl_hm.dat'
                 pow_ka=pows_hm(1,2,:,:)
              END IF           
              outfile=number_file(base,i,ext)
              WRITE(*,*) 'HMx_DRIVER: Output: ', TRIM(outfile)
              
              ! Actually calculate the C(l)
              CALL calculate_Cl(0.,maxdist(proj),ell,Cl,nl,k,a,pow_ka,nk,na,proj,cosm)
              CALL write_Cl(ell,Cl,nl,outfile)
              
           END DO
           WRITE(*,*) 'HMx_DRIVER: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==9) THEN

        ! Breakdown cross-correlation in terms of mass

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        !CALL initialise_cosmology(verbose,cosm)
        CALL print_cosmology(cosm)

        ! Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        !CALL write_projection_kernels(proj,cosm)

        ! Allocate arrays for power
        ALLOCATE(pows_li(nk,na), pows_2h(2,2,nk,na), pows_1h(2,2,nk,na), pows_hm(2,2,nk,na))
        ALLOCATE(pow_ka(nk,na))

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
              ! Set the mass intervals
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

           ! Set the code to not 'correct' the two-halo power for missing
           ! mass when doing the calcultion binned in halo mass
           IF((icumulative .EQV. .TRUE.) .AND. i>0) hmod%ip2h=0
           
           WRITE(*,fmt='(A16)') 'HMx_DRIVER: Mass range'
           WRITE(*,fmt='(A16,I5)') 'HMx_DRIVER: Iteration:', i
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_min [Msun/h]:', m1
           WRITE(*,fmt='(A21,2ES15.7)') 'HMx_DRIVER: M_max [Msun/h]:', m2
           WRITE(*,*)

           !Loop over redshifts
           DO j=1,na

              !z=redshift_a(a(j))

              !Initiliasation for the halomodel calcualtion
              CALL assign_halomod(ihm,hmod,verbose)
              !CALL init_halomod(m1,m2,z,hmod,cosm,verbose)
              CALL init_halomod(m1,m2,a(j),hmod,cosm,verbose)
              CALL print_halomod(hmod,cosm,verbose)
              CALL calculate_HMx_a(ip,2,k,nk,pows_li(:,j),pows_2h(:,:,:,j),pows_1h(:,:,:,j),pows_hm(:,:,:,j),hmod,cosm,verbose,response=.FALSE.)

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
                 outfile=TRIM(dir)//'cl_linear.dat'
                 ext='_cl_linear.dat'
                 pow_ka=pows_li
              ELSE IF(j==2) THEN               
                 outfile=TRIM(dir)//'cl_2h.dat'
                 ext='_cl_2h.dat'
                 pow_ka=pows_2h(1,2,:,:)
              ELSE IF(j==3) THEN               
                 outfile=TRIM(dir)//'cl_1h.dat'
                 ext='_cl_1h.dat'
                 pow_ka=pows_1h(1,2,:,:)
              ELSE IF(j==4) THEN              
                 outfile=TRIM(dir)//'cl_hm.dat'
                 ext='_cl_hm.dat'
                 pow_ka=pows_hm(1,2,:,:)
              END IF

              IF(i>0) outfile=number_file2(base,NINT(log10(m1)),mid,NINT(log10(m2)),ext)

              WRITE(*,*) 'HMx_DRIVER: File: ', TRIM(outfile)

              CALL calculate_Cl(0.,maxdist(proj),ell,Cl,nl,k,a,pow_ka,nk,na,proj,cosm)
              CALL write_Cl(ell,Cl,nl,outfile)

           END DO
           WRITE(*,*) 'HMx_DRIVER: Done'
           WRITE(*,*)

        END DO

     ELSE IF(imode==10) THEN

        ! Break down cross-correlation in terms of redshift

        ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
        CALL print_cosmology(cosm)

        CALL assign_halomod(ihm,hmod,verbose)

        CALL calculate_HMx(ip,2,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)

        ! Allocate array for power
        ALLOCATE(pow_ka(nk,na))

        ! Fill out the projection kernels
        CALL fill_projection_kernels(ix,proj,cosm)
        !CALL write_projection_kernels(proj,cosm)

        ! Write to screen
        WRITE(*,*) 'HMx_DRIVER: Computing C(l)'
        WRITE(*,*) 'HMx_DRIVER: ell min:', REAL(ell(1))
        WRITE(*,*) 'HMx_DRIVER: ell max:', REAL(ell(nl))
        WRITE(*,*) 'HMx_DRIVER: number of ell:', nl
        WRITE(*,*)

        ! Set range of z
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
                 pow_ka=pows_li
              ELSE IF(j==2) THEN
                 ext='_cl_2h.dat'
                 outfile=TRIM(dir)//'cl_2h.dat'
                 pow_ka=pows_2h(1,2,:,:)
              ELSE IF(j==3) THEN
                 ext='_cl_1h.dat'
                 outfile=TRIM(dir)//'cl_1h.dat'
                 pow_ka=pows_1h(1,2,:,:)
              ELSE IF(j==4) THEN
                 ext='_cl_hm.dat'
                 outfile=TRIM(dir)//'cl_hm.dat'
                 pow_ka=pows_hm(1,2,:,:)
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
              CALL calculate_Cl(r1,r2,ell,Cl,nl,k,a,pow_ka,nk,na,proj,cosm)
              CALL write_Cl(ell,Cl,nl,outfile)

           END DO
           WRITE(*,*)

        END DO

        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)

     ELSE

        STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

     END IF

  ELSE IF(imode==12 .OR. imode==44) THEN

     ! Triad stuff
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
     ALLOCATE(Cl(nl))

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
              ix(1)=tracer_KiDS        ! KiDS (z = 0.1 -> 0.9)
              ix(2)=tracer_CMB_lensing ! CMB lensing
              outfile=TRIM(dir)//'/triad_Cl_gal-CMB.dat'
           ELSE IF(itriad==2 .OR. itriad==3) THEN
              ix(1)=tracer_KiDS_450    ! KiDS 450 (z = 0.1 -> 0.9)
              ix(2)=tracer_CMB_lensing ! CMB
              outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.9-CMB.dat'
           ELSE
              STOP 'HMX_DRIVER: Error, itriad specified incorrectly'
           END IF
        ELSE IF(i==2) THEN
           ix(1)=tracer_CMB_lensing ! CMB
           ix(2)=tracer_Compton_y   ! y
           outfile=TRIM(dir)//'/triad_Cl_CMB-y.dat'
        ELSE IF(i==3) THEN
           IF(itriad==1) THEN
              ix(1)=tracer_Compton_y ! y
              ix(2)=tracer_KiDS      ! KiDS (z = 0.1 -> 0.9)
              outfile=TRIM(dir)//'/triad_Cl_y-gal.dat'
           ELSE IF(itriad==2 .OR. itriad==3) THEN
              ix(1)=tracer_Compton_y ! y
              ix(2)=tracer_KiDS_450  ! KiDS 450 (z = 0.1 -> 0.9)
              outfile=TRIM(dir)//'/triad_Cl_y-gal_z0.1-0.9.dat'
           ELSE
              STOP 'HMX_DRIVER: Error, itriad specified incorrectly'
           END IF
        ELSE IF(i==4) THEN
           ix(1)=tracer_KiDS_450_bin1 ! KiDS 450 (z = 0.1 -> 0.5)
           ix(2)=tracer_CMB_lensing   ! CMB
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.5-CMB.dat'
        ELSE IF(i==5) THEN
           ix(1)=tracer_KiDS_450_bin2 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=tracer_CMB_lensing   ! CMB
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-CMB.dat'
        ELSE IF(i==6) THEN
           ix(1)=tracer_Compton_y     ! y
           ix(2)=tracer_KiDS_450_bin1 ! KiDS 450 (z = 0.1 -> 0.5)
           outfile=TRIM(dir)//'/triad_Cl_y-gal_z0.1-0.5.dat'
        ELSE IF(i==7) THEN
           ix(1)=tracer_Compton_y     ! y
           ix(2)=tracer_KiDS_450_bin2 ! KiDS 450 (z = 0.5 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_y-gal_z0.5-0.9.dat'
        ELSE IF(i==8) THEN
           ix(1)=tracer_KiDS_450_bin1 ! KiDS 450 (z = 0.1 -> 0.5)
           ix(2)=tracer_KiDS_450_bin1 ! KiDS 450 (z = 0.1 -> 0.5)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.5-gal_z0.1-0.5.dat'
        ELSE IF(i==9) THEN
           ix(1)=tracer_KiDS_450_bin1 ! KiDS 450 (z =  0.1 -> 0.5)
           ix(2)=tracer_KiDS_450      ! KiDS 450 (z =  0.1 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.5-gal_z0.1-0.9.dat'
        ELSE IF(i==10) THEN
           ix(1)=tracer_KiDS_450 ! KiDS 450 (z = 0.1 -> 0.9)
           ix(2)=tracer_KiDS_450 ! KiDS 450 (z = 0.1 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.1-0.9-gal_z0.1-0.9.dat'
        ELSE IF(i==11) THEN
           ix(1)=tracer_KiDS_450_bin2 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=tracer_KiDS_450_bin1 ! KiDS 450 (z = 0.1 -> 0.5)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-gal_z0.1-0.5.dat'
        ELSE IF(i==12) THEN
           ix(1)=tracer_KiDS_450_bin2 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=tracer_KiDS_450      ! KiDS 450 (z = 0.1 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-gal_z0.1-0.9.dat'
        ELSE IF(i==13) THEN
           ix(1)=tracer_KiDS_450_bin2 ! KiDS 450 (z = 0.5 -> 0.9)
           ix(2)=tracer_KiDS_450_bin2 ! KiDS 450 (z = 0.5 -> 0.9)
           outfile=TRIM(dir)//'/triad_Cl_gal_z0.5-0.9-gal_z0.5-0.9.dat'
        END IF

        ! Do the cross correlation
        CALL xcorr(ix,mmin,mmax,ell,Cl,nl,hmod,cosm,verbose)
        CALL write_Cl(ell,Cl,nl,outfile)

        WRITE(*,*) 'HMx_DRIVER: Done'
        WRITE(*,*)
        
     END DO

  ELSE IF(imode==13) THEN

     ! 13 - Calculate a cross-correlation coefficient

     ! Assign the cosmology
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)

     ! Assign the halo model
     CALL assign_halomod(ihm,hmod,verbose)

     ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     CALL print_cosmology(cosm)

     ! Set the ell range and allocate arrays for l and C(l)
     lmin=1e0
     lmax=1e4 ! Strange errors and crashes if this is increased to 10^5
     nl=64 
     CALL fill_array(log(lmin),log(lmax),ell,nl)
     ell=exp(ell)
     ALLOCATE(Cl(nl))

     ! Directory for outputs
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
           outfile=TRIM(dir)//'/cl_hm.dat'
        END IF

        ! Do the cross correlation
        CALL xcorr(ix,mmin,mmax,ell,Cl,nl,hmod,cosm,verbose)
        CALL write_Cl(ell,Cl,nl,outfile)
        
     END DO

  ELSE IF(imode==14 .OR. imode==33) THEN

     ! Make power spectra variations as a function of baryon parameter variations

     ! Number of values to try for each parameter
     m=9

     ! Set number of k points and k range (log spaced)
     kmin=1e-3
     kmax=1e1
     nk=128
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Set the fields
     nf=5
     ALLOCATE(fields(nf))
     fields(1)=field_matter
     fields(2)=field_cdm
     fields(3)=field_gas
     fields(4)=field_star
     fields(5)=field_electron_pressure

     ! Allocate arrays for the power spectra
     ALLOCATE(powd_li(nk),powd_2h(nk),powd_1h(nk),powd_hm(nk))
     ALLOCATE(pow_li(nk),pow_2h(nf,nf,nk),pow_1h(nf,nf,nk),pow_hm(nf,nf,nk))

     ! Set the redshift
     z=0.

     ! Assigns the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)

     ! Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
     !CALL initialise_cosmology(verbose,cosm)
     CALL print_cosmology(cosm)

     ! Initiliasation for the halo-model calcualtion
     IF(imode==33) ihm=3
     CALL assign_halomod(ihm,hmod,verbose)
     CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     ! DMONLY
     field=field_dmonly
     outfile='data/DMONLY.dat'
     CALL calculate_HMx_a(field,1,k,nk,powd_li,powd_2h,powd_1h,powd_hm,hmod,cosm,verbose,response=.FALSE.)
     CALL write_power(k,powd_li,powd_2h,powd_1h,powd_hm,nk,outfile,verbose)

     ! Prevents warning
     ilog=.FALSE.

     ! Number of parameters
     npa=10

     ! Loop over parameters     
     DO ipa=10,npa
        
        ! Set maximum and minimum parameter values and linear or log range
        IF(ipa==1) THEN
           ! alpha - virial temperature pre factor
           param_min=0.05
           param_max=0.65
           ilog=.FALSE.
        ELSE IF(ipa==2) THEN
           ! epsilon - concentration change due to gas
           param_min=0.5
           param_max=2.0
           ilog=.FALSE.
        ELSE IF(ipa==3) THEN
           ! Gamma - KS polytropic index
           param_min=1.12
           param_max=1.22
           ilog=.FALSE.
        ELSE IF(ipa==4) THEN
           ! M0 - bound gas transition in Msun/h
           param_min=1e13
           param_max=1e15
           ilog=.TRUE.
        ELSE IF(ipa==5) THEN
           ! A* - Stellar mass fraction
           param_min=0.02
           param_max=0.04
           ilog=.FALSE.
        ELSE IF(ipa==6) THEN
           ! WHIM temperature in K
           param_min=1e5
           param_max=1e7
           ilog=.TRUE.
        ELSE IF(ipa==7) THEN
           ! c* - stellar concentation
           param_min=10.
           param_max=100.
           ilog=.TRUE.
        ELSE IF(ipa==8) THEN
           ! fcold - fraction of bound gas that is cold
           param_min=0.0
           param_max=0.25
           ilog=.FALSE.
        ELSE IF(ipa==9) THEN
           ! alphap - alpha mass index
           param_min=-0.1
           param_max=0.1
           ilog=.FALSE.
        ELSE IF(ipa==10) THEN
           ! Gammap - Gamma mass index
           param_min=-0.01
           param_max=0.01
           ilog=.FALSE.
        END IF

        ! Loop over parameter values
        DO i=1,m

           ! Set the parameter value that is being varied
           IF(ilog) THEN
              param=progression_log(param_min,param_max,i,m)
           ELSE
              param=progression(param_min,param_max,i,m)
           END IF

           CALL assign_halomod(ihm,hmod,verbose)

           IF(hmod%fixed_HMx) THEN
              IF(ipa==1)  hmod%alpha=param
              IF(ipa==2)  hmod%eps=param
              IF(ipa==3)  hmod%Gamma=param
              IF(ipa==4)  hmod%M0=param
              IF(ipa==5)  hmod%Astar=param
              IF(ipa==6)  hmod%Twhim=param
              IF(ipa==7)  hmod%cstar=param
              IF(ipa==8)  hmod%fcold=param
              IF(ipa==9)  hmod%alphap=param
              IF(ipa==10) hmod%Gammap=param
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
              ELSE
                 STOP 'HMx_DRIVER: Error, unfixed not supported for this parameter'
              END IF
           END IF
           
           CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
           CALL print_halomod(hmod,cosm,verbose)

           ! DO NOT DELETE THIS
           ! It is only used to print values to the screen later
           ! For example, mass, which is inconvenient if written out in full
           IF(ilog) THEN
              param_neat=log10(param)
           ELSE
              param_neat=param
           END IF

           ! Write out halo matter and electron-pressure profile information
           ! All the string crap is in the loop for a reason
           DO j=10,16
              base='data/profile_mass_'
              ext='_param_'
              base=number_file(base,j,ext)
              mid='_value_'
              ext='.dat'
              outfile=number_file2(base,ipa,mid,i,ext)
              mass=10.**j
              CALL write_halo_profiles(mass,hmod,cosm,outfile)
           END DO

           ! Write out halo mass fraction information
           base='data/mass_fractions_param_'
           outfile=number_file2(base,ipa,mid,i,ext)
           CALL write_mass_fractions(hmod,cosm,outfile)           

           CALL calculate_HMx_a(fields,nf,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response=.FALSE.)
           
           ! Loop over field combinations and do the calculation
           DO j1=1,nf
              DO j2=j1,nf

                 ! Set output file
                 base='data/power_param_'
                 mid='_value_'
                 ext='_'
                 outfile=number_file2(base,ipa,mid,i,ext)

                 mid=''
                 ext='.dat'
                 outfile=number_file2(outfile,fields(j1),mid,fields(j2),ext)

                 ! Write progress to screen
                 WRITE(*,fmt='(4I5,F14.7,A50)') ipa, i, fields(j1), fields(j2), param_neat, TRIM(outfile)

                 ! Do the halo-model calculation and write to file
                 CALL write_power(k,pow_li,pow_2h(j1,j2,:),pow_1h(j1,j2,:),pow_hm(j1,j2,:),nk,outfile,verbose=.FALSE.)

              END DO
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
     nf=2
     DO i=1,nf
        WRITE(*,*) 'HMx_driver: Choose halo', i
        CALL set_halo_type(ip(i))
     END DO

     ! User chooses halo model   
     CALL calculate_HMx(ip,nf,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)

     base='data/power'
     CALL write_power_a_multiple(k,a,pows_li,pows_2h(ip(1),ip(2),:,:),pows_1h(ip(1),ip(2),:,:),pows_hm(ip(1),ip(2),:,:),nk,na,base,verbose)

  ELSE IF(imode==18 .OR. imode==48) THEN

     ! Create 3D bias function
     ! 18 - User choice
     ! 48 - HI bias

     ! Assign the cosmological model
     IF(imode==48) icosmo=27 ! Illustris TNG 75
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

     ! Allocate arrays for fields
     nf=2
     ALLOCATE(fields(nf))
     fields(1)=field_dmonly

     ! Select field type for bias study
     IF(imode==18) CALL set_halo_type(fields(2))
     IF(imode==48) fields(2)=12 ! HI

     CALL calculate_HMx(fields,nf,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)

     DO j1=1,nf
        DO j2=j1,nf

           IF(j1==1 .AND. j2==1) THEN
              ! DMONLY-DMONLY
              base='data/power_mm'
           ELSE IF(j1==1 .AND. j2==2) THEN
              ! DMONLY-field
              base='data/power_mf'
           ELSE IF(j1==2 .AND. j2==2) THEN
              ! field-field
              base='data/power_ff'
           END IF
           
           CALL write_power_a_multiple(k,a,pows_li,pows_2h(j1,j2,:,:),pows_1h(j1,j2,:,:),pows_hm(j1,j2,:,:),nk,na,base,verbose)

        END DO
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
     ALLOCATE(pow_li(nk),pow_2h(1,1,nk),pow_1h(1,1,nk),pow_hm(1,1,nk))

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
           CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
           CALL print_halomod(hmod,cosm,verbose)

           ! Do the halo-model calculation
           field=field_dmonly
           CALL calculate_HMx_a(field,1,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response=.FALSE.)

           ! Write out the results
           CALL write_power(k,pow_li,pow_2h(1,1,:),pow_1h(1,1,:),pow_hm(1,1,:),nk,outfile,verbose)

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
     CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
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

     ! Mead (2017) dark energy results

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e1
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)
     ALLOCATE(pow_li(nk),pow_2h(1,1,nk),pow_1h(1,1,nk),pow_hm(1,1,nk))

     ! Set the redshift
     z=0.

     ! Directory for output
     dir='data'

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
        CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Do the halo-model calculation
        field=field_dmonly
        CALL calculate_HMx_a(field,1,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response=.FALSE.)

        ! Write out the results
        outfile=TRIM(dir)//'/power_'//TRIM(outfile)//'.dat'
        CALL write_power(k,pow_li,pow_2h(1,1,:),pow_1h(1,1,:),pow_hm(1,1,:),nk,outfile,verbose)

     END DO

  ELSE IF(imode==24) THEN

     STOP 'HMx_DRIVER: Error, imode=24 no longer supported'

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
        ALLOCATE(k(nk),a(na),pow_ka(nk,na))

        ! Read-in test data
        OPEN(7,file=infile,status='old')
        DO i=0,nk
           IF(i==0) THEN
              READ(7,*) crap, (a(j), j=1,na)
           ELSE
              READ(7,*) k(i), (pow_ka(i,j), j=1,na)
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

        field=field_dmonly
        CALL calculate_HMx(field,1,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose_tests,response=.FALSE.)
        
        ! Loop over k and a
        error_max=0.
        DO j=1,na
           DO i=1,nk
              error=ABS(-1.+pows_hm(1,1,i,j)/pow_ka(i,j))
              IF(error>error_max) error_max=error
              IF(error>tolerance) THEN
                 WRITE(*,*) 'HMx_DRIVER: Test:', itest
                 WRITE(*,*) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                 WRITE(*,*) 'HMx_DRIVER: Scale-factor:', a(j)
                 WRITE(*,*) 'HMx_DRIVER: Expected power:', pow_ka(i,j)
                 WRITE(*,*) 'HMx_DRIVER: Model power:', pows_hm(1,1,i,j)
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
        !CALL write_power_a_multiple(k,a,powa_lin,powa_2h,powa_1h,powa_full,nk,na,base,verbose_tests)
        CALL write_power_a_multiple(k,a,pows_li,pows_2h(1,1,:,:),pows_1h(1,1,:,:),pows_hm(1,1,:,:),nk,na,base,verbose_tests)

        DEALLOCATE(k,a,pow_ka)
        
     END DO

     IF(ifail) THEN
        STOP 'HMx_DRIVER: Error, tests failed'
     ELSE
        WRITE(*,*) 'HMx_DRIVER: Tests should take around 0.55 seconds to run'
        WRITE(*,*) 'HMx_DRIVER: Tests passed'
        WRITE(*,*)
     END IF

  ELSE IF(imode==27 .OR. imode==28 .OR. imode==29 .OR. imode==30)  THEN

     ! Comparison with FrankenEmu or Mira Titan

     ! Number of cosmological models (+1)
     IF(imode==28 .OR. imode==30) n=37
     IF(imode==27 .OR. imode==29) n=10

     ! Set redshifts/scale factors
     nz=4
     na=nz
     ALLOCATE(zs(nz),a(na))
     zs(1)=0.0
     zs(2)=0.5
     zs(3)=1.0
     zs(4)=2.0
     DO i=1,na
        a(i)=scale_factor_z(zs(i))
     END DO

     base='data/cosmo'
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

           IF(imode==27 .OR. imode==29) CALL read_Mira_Titan_power(k_sim,pow_sim,nk,zs(j),cosm,rebin=.FALSE.)
           IF(imode==28 .OR. imode==30) CALL read_FrankenEmu_power(k_sim,pow_sim,nk,zs(j),cosm,rebin=.FALSE.)

           ALLOCATE(pow_li(nk),pow_2h(1,1,nk),pow_1h(1,1,nk),pow_hm(1,1,nk))

           CALL init_halomod(mmin,mmax,a(j),hmod,cosm,verbose=.FALSE.)
           CALL print_halomod(hmod,cosm,verbose=.FALSE.)
           field=field_dmonly
           CALL calculate_HMx_a(field,1,k_sim,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

           ! Write data to disk
           outfile=number_file2(base,i,mid,j,ext)
           OPEN(7,file=outfile)
           DO ii=1,nk
              WRITE(7,*) k_sim(ii), pow_hm(1,1,ii), pow_sim(ii)
           END DO
           CLOSE(7)

           DEALLOCATE(pow_li,pow_2h,pow_1h,pow_hm)
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

     ! Set the fields
     nf=2
     ALLOCATE(fields(2))
     fields(1)=field_matter
     fields(2)=field_electron_pressure

     ! Allocate arrays for power
     ALLOCATE(pow_li(nk),pow_2h(nf,nf,nk),pow_1h(nf,nf,nk),pow_hm(nf,nf,nk))

     ! Assign the cosmological model
     icosmo=4 ! BAHAMAS cosmology
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the halo model
     ihm=3
     CALL assign_halomod(ihm,hmod,verbose)

     ! Set the redshift
     z=0.

     ! Loop over upper limit of mass integral
     DO i=10,16

        ! Set the upper limit for the mass integration
        ! Needs to be 10. to enforce that m2 is real
        m2=10.**i

        ! Initiliasation for the halomodel calcualtion
        CALL init_halomod(mmin,m2,scale_factor_z(z),hmod,cosm,verbose)
        CALL print_halomod(hmod,cosm,verbose)

        ! Do the halo-model calculation
        CALL calculate_HMx_a(fields,nf,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response=.FALSE.)

        ! Loop over fields
        DO j1=1,nf
           DO j2=j1,nf

              ! Write out the results
              base='data/power_'
              mid=''
              ext='_m'
              base=number_file2(base,fields(j1),mid,fields(j2),ext)
              ext='.dat'
              outfile=number_file(base,i,ext)
              CALL write_power(k,pow_li,pow_2h(j1,j2,:),pow_1h(j1,j2,:),pow_hm(j1,j2,:),nk,outfile,verbose)

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
     CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
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
           WRITE(7,*) hmod%z, HMx_alpha(1e14,hmod), HMx_eps(hmod), HMx_Gamma(1e14,hmod), HMx_M0(hmod), HMx_Astar(hmod), HMx_Twhim(hmod)
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

     ! Allocate arrays for the power calculation
     ALLOCATE(pow_li(nk),pow_2h(1,1,nk),pow_1h(1,1,nk),pow_hm(1,1,nk))

     ! Assign the cosmological model
     icosmo=1 ! Boring
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the redshift
     z=0.

     ! Initiliasation for the halomodel calcualtion
     ihm=21
     CALL assign_halomod(ihm,hmod,verbose)

     ! Range of cores to explore
     rcore_min=0.
     rcore_max=0.1
     ncore=16
     
     DO i=1,ncore

        ! Set the core radius
        hmod%rcore=progression(rcore_min,rcore_max,i,ncore)

        ! Initialise the halo-model calculation
        IF(i==1) THEN
           verbose2=verbose
        ELSE
           verbose2=.FALSE.
        END IF
        CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose2)
        CALL print_halomod(hmod,cosm,verbose2)
        
        ! Do the halo-model calculation
        field=field_dmonly
        CALL calculate_HMx_a(field,1,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose=.FALSE.,response=.FALSE.)

        ! Write out the results
        base='data/power_cored_'
        ext='.dat'
        outfile=number_file(base,i,ext)
        CALL write_power(k,pow_li,pow_2h(1,1,:),pow_1h(1,1,:),pow_hm(1,1,:),nk,outfile,verbose)

     END DO

  ELSE IF(imode==36) THEN

     ! Automated testing of hydro models

     ! Assigns the cosmological model
     icosmo=4
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the field combinations
     nf=5
     ALLOCATE(fields(nf))
     fields(1)=field_matter
     fields(2)=field_cdm
     fields(3)=field_gas
     fields(4)=field_star
     fields(5)=field_electron_pressure

     ! Set the redshifts
     na=4
     ALLOCATE(a(na))
     DO i=1,na
        IF(i==1) z=0.0
        IF(i==2) z=0.5
        IF(i==3) z=1.0
        IF(i==4) z=2.0
        a(i)=scale_factor_z(z)
     END DO

     ! File naming things
     mid=''
     inext='.txt'

     !! Read benchmark data

     ! Loop over fields
     DO itest=1,nf
        DO jtest=itest,nf

           ! Loop over redshifts
           DO j=1,na

              ! Set the redshift and input and output file bases
              IF(j==1) THEN
                 inbase='benchmarks/power_z0.0_'
              ELSE IF(j==2) THEN
                 inbase='benchmarks/power_z0.5_'
              ELSE IF(j==3) THEN
                 inbase='benchmarks/power_z1.0_'
              ELSE IF(j==4) THEN
                 inbase='benchmarks/power_z2.0_'
              ELSE
                 STOP 'HMX_DRIVER: Error, iz specified incorrectly'
              END IF

              ! Input file name
              infile=number_file2(inbase,fields(itest),mid,fields(jtest),inext)

              ! Allocate arrays            
              IF(itest==1 .AND. jtest==1 .AND. j==1) THEN
                 nk=file_length(infile,verbose=.FALSE.)
                 ALLOCATE(k(nk),powb_hm(nf,nf,nk,na))
              END IF

              ! Loop over k and read in benchmark data
              ! TODO: Add tests for one- and two-halo terms individually
              OPEN(7,file=infile,status='old')
              DO i=1,nk
                 READ(7,*) k(i), spam, spam, spam, powb_hm(itest,jtest,i,j)
              END DO
              CLOSE(7)
        
           END DO
        END DO
     END DO

     !!
     
     ! Assign the halo model
     ihm=3
     CALL assign_halomod(ihm,hmod,verbose)
     CALL calculate_HMx(fields,nf,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)
     
     ! Initially assume that all the tests will pass
     ifail=.FALSE.

     ! File naming things
     outext='.dat'

     ! Loop over fields
     DO itest=1,nf
        DO jtest=itest,nf

           ! Loop over redshift
           DO j=1,na

              ! Set the redshift and input and output file bases
              IF(j==1) THEN
                 outbase='data/power_z0.0_'
              ELSE IF(j==2) THEN
                 outbase='data/power_z0.5_'
              ELSE IF(j==3) THEN
                 outbase='data/power_z1.0_'
              ELSE IF(j==4) THEN
                 outbase='data/power_z2.0_'
              ELSE
                 STOP 'HMX_DRIVER: Error, iz specified incorrectly'
              END IF

              ! Write data to file
              outfile=number_file2(outbase,fields(itest),mid,fields(jtest),outext)
              CALL write_power(k,pows_li(:,j),pows_2h(itest,jtest,:,j),pows_1h(itest,jtest,:,j),pows_hm(itest,jtest,:,j),nk,outfile,verbose)

              ! Set the error max to be zero before each test
              error_max=0.
              verbose_tests=.TRUE.
              
              ! Loop over k and check for accuracy and write out data
              DO i=1,nk

                 ! This is what I use as the error
                 error=ABS(-1.+pows_hm(itest,jtest,i,j)/powb_hm(itest,jtest,i,j))

                 ! Update the maximium error if it is exceeded
                 IF(error>error_max) error_max=error

                 ! If the test has failed then write out this diagnostic stuff
                 IF(verbose_tests .AND. error>tolerance) THEN
                    WRITE(*,*) 'HMx_DRIVER: Test failing'
                    WRITE(*,*) 'HMx_DRIVER: Test fields:', fields(itest), fields(jtest)
                    WRITE(*,*) 'HMx_DRIVER: Wavenumber [h/Mpc]:', k(i)
                    WRITE(*,*) 'HMx_DRIVER: Redshift:', redshift_a(a(j))
                    WRITE(*,*) 'HMx_DRIVER: Benchmark power:', powb_hm(itest,jtest,i,j)
                    WRITE(*,*) 'HMx_DRIVER: HMx power:', pows_hm(itest,jtest,i,j)
                    WRITE(*,*) 'HMx_DRIVER: Tolerance:', tolerance
                    WRITE(*,*) 'HMx_DRIVER: Error:', error
                    WRITE(*,*)
                    verbose_tests=.FALSE.
                    ifail=.TRUE.
                 END IF
                 
              END DO

              ! Write test information to the screen
              WRITE(*,*) 'HMx_DRIVER: Test fields:', fields(itest), fields(jtest)
              WRITE(*,*) 'HMx_DRIVER: Redshift:', redshift_a(a(j))
              WRITE(*,*) 'HMx_DRIVER: Tolerance:', tolerance
              WRITE(*,*) 'HMx_DRIVER: Max error:', error_max
              WRITE(*,*)

           END DO

        END DO
     END DO

     ! Write pass/fail information to the screen
     IF(ifail) THEN
        WRITE(*,*) 'HMx_DRIVER: Hydro tests failed'
     ELSE
        WRITE(*,*) 'HMx_DRIVER: Hydro tests should take around 3.00 seconds to run'
        WRITE(*,*) 'HMx_DRIVER: Hydro tests passed'        
     END IF
     WRITE(*,*)

  ELSE IF(imode==45) THEN

     ! Comparison of power spectra from Sheth & Tormen (1999) vs. Tinker et al. (2010) mass function

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Allocate arrays for power
     ALLOCATE(pow_li(nk),pow_2h(1,1,nk),pow_1h(1,1,nk),pow_hm(1,1,nk))

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Sets the redshift
     z=0.

     DO j=1,3

        IF(j==1) THEN
           ! Sheth & Tormen (1999) mass function and bias
           ihm=3
           outfile='data/power_ShethTormen.dat'
        ELSE IF(j==2) THEN
           ihm=23 ! Tinker et al. (2010) mass function and bias
           outfile='data/power_Tinker.dat'
        ELSE IF(j==3) THEN
           ihm=27 ! Press & Schecter (1974) mass function and bias
           outfile='data/power_PressSchecter.dat'
        END IF

        IF(j==1) THEN
           verbose2=verbose
        ELSE
           verbose2=.FALSE.
        END IF
        
        !Initiliasation for the halomodel calcualtion
        CALL assign_halomod(ihm,hmod,verbose2)
        CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose2)
        CALL print_halomod(hmod,cosm,verbose2)

        ! Do the halo-model calculation
        field=field_dmonly
        CALL calculate_HMx_a(field,1,k,nk,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose2,response=.FALSE.)

        ! Write out data to disk
        CALL write_power(k,pow_li,pow_2h(1,1,:),pow_1h(1,1,:),pow_hm(1,1,:),nk,outfile,verbose)

     END DO

  ELSE IF(imode==46) THEN

     ! Mass function plots for different mass functions

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
     CALL init_halomod(mmin,mmax,scale_factor_z(z),hmod,cosm,verbose)
     CALL print_halomod(hmod,cosm,verbose)

     ! Range in nu
     nu_min=0.1
     nu_max=6.
     n=256

     ! Loop over nu and write out mass function
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

     ! HI mass fractions
     
     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set the range of redshifts
     z1=0.
     z2=5.
     nz=6
     ALLOCATE(hmods(nz),HI_frac(nz),z_tab(nz),a(nz))
     
     ! Initiliasation for the halomodel calcualtion at all zs
     ihm=25
     DO j=1,nz
        z_tab(j)=progression(z1,z2,j,nz)
        a(j)=scale_factor_z(z_tab(j))
        CALL assign_halomod(ihm,hmods(j),verbose)
        CALL init_halomod(mmin,mmax,a(j),hmods(j),cosm,verbose)
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

  ELSE IF(imode==50) THEN

     ! Mass function plots as Lbox is varied

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose=.FALSE.)

     ! Minimum and maximum box sizes [Mpc/h] and number of cosmologies
     Lmin=32.
     Lmax=2048.
     ncos=8

     ! Set cosmological models
     ALLOCATE(cosms(ncos))
     cosms=cosm
     DO i=1,ncos
        IF(i .NE. 1) THEN
           cosms(i)%box=.TRUE.
        END IF
        cosms(i)%Lbox=2**(12-i)
        CALL init_cosmology(cosms(i))
     END DO
     CALL print_cosmology(cosms(1))

     ! Set the redshift
     z=0.

     !Initilisation for the halomodel calcualtion
     ihm=3
     ALLOCATE(hmods(ncos))
     DO i=1,ncos
        IF(i==1) THEN
           verbose2=verbose
        ELSE
           verbose2=.FALSE.
        END IF
        CALL assign_halomod(ihm,hmods(i),verbose2)
        CALL init_halomod(mmin,mmax,scale_factor_z(z),hmods(i),cosms(i),verbose2)        
     END DO
     CALL print_halomod(hmods(1),cosms(1),verbose)

     ! Range in nu
     nu_min=0.1
     nu_max=6.
     n=256

     ALLOCATE(masses(ncos))
     OPEN(7,file='data/nu_mass_Lbox.dat')
     DO i=1,n
        DO j=1,ncos
           nu=progression(nu_min,nu_max,i,n)
           masses(j)=M_nu(nu,hmods(j))
           !mf, mf*mass/comoving_matter_density(cosm), mf*mass**2/comoving_matter_density(cosm)
        END DO
        WRITE(7,*) nu, (masses(j), j=1,ncos)
     END DO
     CLOSE(7)

  ELSE IF(imode==51) THEN

     ! Assigns the cosmological model
     icosmo=1
     CALL assign_cosmology(icosmo,cosm,verbose)
     CALL init_cosmology(cosm)
     CALL print_cosmology(cosm)

     ! Set number of k points and k range (log spaced)
     nk=128
     kmin=1e-3
     kmax=1e2
     CALL fill_array(log(kmin),log(kmax),k,nk)
     k=exp(k)

     ! Set the number of redshifts and range (linearly spaced) and convert z -> a
     na=16
     amin=0.1
     amax=1.
     CALL fill_array(amin,amax,a,na)

     ! Set the field
     field=field_dmonly

     DO i=1,2

        ! Assign halo model
        IF(i==1) THEN
           ihm=3
           base='data/power'
        ELSE IF(i==2) THEN
           ihm=8
           base='data/power_scatter'
        ELSE
           STOP 'HMX_DRIVER: Error, something went wrong'
        END IF
           
        CALL assign_halomod(ihm,hmod,verbose)

        CALL calculate_HMx(field,1,mmin,mmax,k,nk,a,na,pows_li,pows_2h,pows_1h,pows_hm,hmod,cosm,verbose,response=.FALSE.)
        
        CALL write_power_a_multiple(k,a,pows_li,pows_2h(1,1,:,:),pows_1h(1,1,:,:),pows_hm(1,1,:,:),nk,na,base,verbose)

     END DO
     
  ELSE
        
     STOP 'HMx_DRIVER: Error, you have specified the mode incorrectly'

  END IF

CONTAINS

  SUBROUTINE read_k_values(infile,k,nk)

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
    
  END SUBROUTINE read_k_values

  SUBROUTINE YinZhe_Fig1(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmological model
    INTEGER :: i
    REAL :: r, rs, rv, c, Mh, rh, r500c, m500c

    REAL, PARAMETER :: M=1e15    ! Halo virial? mass [Msun]
    REAL, PARAMETER :: rmin=1e-3 ! Minimum radius [Mpc]
    REAL, PARAMETER :: rmax=8    ! Maximum radius [Mpc] 
    INTEGER, PARAMETER :: nr=512 ! Number of points in radius

    LOGICAL, PARAMETER :: real_space=.TRUE.
    INTEGER, PARAMETER :: itype=field_electron_pressure ! electron pressure
    !INTEGER, PARAMETER :: ipnh=1 ! Type of term to consider

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod,cosm)

    Mh=M*cosm%h ! This virial mass is now [Msun/h]

    rv=exp(find(log(Mh),hmod%log_m,log(hmod%rv),hmod%n,3,3,2)) ! [Mpc/h]
    c=find(log(Mh),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c ! [Mpc/h]

    m500c=exp(find(log(Mh),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ! [Mpc/h]
    r500c=exp(find(log(Mh),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ! [Mpc/h]

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
       r=progression(rmin,rmax,i,nr) ! Radius [Mpc]
       rh=r*cosm%h ! Convert [Mpc/h]
       WRITE(7,*) r, UPP(real_space,rh,Mh,rv,rs,hmod,cosm)*r**2, win_type(real_space,itype,rh,Mh,rv,rs,hmod,cosm)*r**2
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
          output=TRIM(base)//'_2h.dat'
          pow=pow_2h
       ELSE IF(i==3) THEN
          output=TRIM(base)//'_1h.dat'
          pow=pow_1h
       ELSE IF(i==4) THEN
          output=TRIM(base)//'_hm.dat'
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

  SUBROUTINE xcorr(ix,mmin,mmax,ell,Cl,nl,hmod,cosm,verbose)

    ! Calculates the C(l) for the cross correlation of fields ix(1) and ix(2)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    REAL, INTENT(IN) :: mmin, mmax
    REAL, INTENT(IN) :: ell(nl)
    REAL, INTENT(OUT) :: Cl(nl)
    INTEGER, INTENT(IN) :: nl
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    REAL, ALLOCATABLE :: a(:), k(:), pow_li(:,:), pow_2h(:,:,:,:), pow_1h(:,:,:,:), pow_hm(:,:,:,:)
    TYPE(projection) :: proj(2)
    REAL :: lmin, lmax
    INTEGER :: ip(2), nk, na
    REAL :: r1, r2

    ! Set the k range
    nk=nk_xcorr
    CALL fill_array(log(kmin_xcorr),log(kmax_xcorr),k,nk)
    k=exp(k)   

    ! Set the a range
    na=na_xcorr
    CALL fill_array(amin_xcorr,amax_xcorr,a,na)

    ! Set the ell range
    lmin=ell(1)
    lmax=ell(nl)

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Cross-correlation information'
       WRITE(*,*) 'XCORR: P(k) minimum k [h/Mpc]:', REAL(kmin_xcorr)
       WRITE(*,*) 'XCORR: P(k) maximum k [h/Mpc]:', REAL(kmax_xcorr)
       WRITE(*,*) 'XCORR: Number of k:', nk
       WRITE(*,*) 'XCORR: minimum a:', REAL(amin_xcorr)
       WRITE(*,*) 'XCORR: maximum a:', REAL(amax_xcorr)
       WRITE(*,*) 'XCORR: number of a:', na
       WRITE(*,*) 'XCORR: minimum ell:', REAL(lmin)
       WRITE(*,*) 'XCORR: maximum ell:', REAL(lmax)
       WRITE(*,*) 'XCORR: number of ell:', nl
       WRITE(*,*)
    END IF

    ! Use the xcorrelation type to set the necessary halo profiles
    CALL set_xcorr_type(ix,ip)

    ! Do the halo model power spectrum calculation
    CALL calculate_HMx(ip,2,mmin,mmax,k,nk,a,na,pow_li,pow_2h,pow_1h,pow_hm,hmod,cosm,verbose,response=.FALSE.)

    ! Fill out the projection kernels
    CALL fill_projection_kernels(ix,proj,cosm)
    !IF(verbose) CALL write_projection_kernels(proj,cosm)

    ! Set the range in comoving distance for the Limber integral
    r1=0.
    r2=maxdist(proj)

    ! Actually calculate the C(ell), but only for the full halo model part
    ! TODO: Array temporary
    CALL calculate_Cl(r1,r2,ell,Cl,nl,k,a,pow_hm(1,2,:,:),nk,na,proj,cosm)

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE xcorr

  SUBROUTINE set_xcorr_type(ix,ip)

    ! Set the cross-correlation type
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    INTEGER, INTENT(OUT) :: ip(2)
    INTEGER :: i, j

    ! Loop over two-components of xcorr
    DO i=1,2

       IF(ix(i)==-1) THEN
          WRITE(*,fmt='(A30,I3)') 'SET_XCORR_TYPE: Choose field: ', i
          WRITE(*,*) '========================='
          DO j=1,n_tracers
             WRITE(*,fmt='(I3,A3,A30)') j, '- ', TRIM(xcorr_type(j))
          END DO
          READ(*,*) ix(i)
          WRITE(*,*) '========================='
          WRITE(*,*)
       END IF

       IF(ix(i)==tracer_Compton_y) THEN
          ! Compton y
          ip(i)=field_electron_pressure
       ELSE IF(ix(i)==tracer_gravity_wave) THEN
          ! Gravitational waves
          ip(i)=field_matter
       ELSE
          ! Gravitational lensing
          ip(i)=field_matter
       END IF
       
    END DO

  END SUBROUTINE set_xcorr_type
  
END PROGRAM HMx_driver
