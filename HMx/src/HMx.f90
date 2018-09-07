MODULE HMx

  !Module usage statements
  USE constants
  USE solve_equations
  USE special_functions
  USE string_operations
  USE calculus_table
  USE cosmology_functions
  
  IMPLICIT NONE

  !Halo-model stuff that needs to be recalculated for each new z
  TYPE halomod
     INTEGER :: ip2h, ibias, imf, iconc, iDolag, iAs, ip2h_corr, ikb
     INTEGER :: idc, iDv, ieta, ikstar, i2hdamp, i1hdamp, itrans, iscatter
     LOGICAL :: voids, use_UPP, smooth_freegas
     REAL :: z
     REAL :: alpha, eps, Gamma, M0, Astar, Twhim, rstar, sstar, mstar, Theat ! HMx baryon parameters
     REAL :: A_alpha, B_alpha, C_alpha, D_alpha
     REAL :: A_eps, B_eps, C_eps, D_eps
     REAL :: A_Gamma, B_Gamma, C_Gamma, D_Gamma
     REAL :: A_M0, B_M0, C_M0, D_M0
     REAL :: A_Astar, B_Astar, C_Astar, D_Astar
     REAL :: A_Twhim, B_Twhim, C_Twhim, D_Twhim
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:), log_m(:)
     REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
     REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
     REAL :: sigv, sigv100, c3, knl, rnl, mnl, neff, sig8z
     REAL :: gmin, gmax, gbmin, gbmax
     REAL :: n_c, n_s, n_g, rho_HI
     REAL :: Dv0, Dv1, dc0, dc1, eta0, eta1, f0, f1, ks, A, alp0, alp1 ! HMcode parameters
     REAL :: mgal, HImin, HImax ! HOD parameters
     INTEGER :: n
     LOGICAL :: has_HI, has_galaxies, has_mass_conversions, safe_negative, has_dewiggle
     LOGICAL :: fixed_HMx
     !LOGICAL :: verbose
     REAL :: acc_HMx, large_nu
     CHARACTER(len=256) :: name
     REAL, ALLOCATABLE :: log_k_pdamp(:), log_pdamp(:)
     INTEGER :: n_pdamp
  END TYPE halomod

  !Global parameters
  REAL, PARAMETER :: acc_win=1e-3 ! Halo-model accuracy parameter
  INTEGER, PARAMETER :: imeth_win=3 ! Winint method

CONTAINS

  SUBROUTINE dewiggle_init(z,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: kv(4), pv(4), sigv, a
    !REAL, ALLOCATABLE :: log_Pk(:), log_Pkraw(:), log_k(:)
    REAL, ALLOCATABLE :: Pk(:), Pkraw(:), k(:)
    INTEGER :: i, nk

    IF(.NOT. ALLOCATED(cosm%log_k_plin)) STOP 'DEWIGGLE_INIT: Error, P(k) needs to be tabulated for this to work'

    nk=cosm%n_plin
    hmod%n_pdamp=nk
    sigv=hmod%sigv
    a=scale_factor_z(z)

    ! Allocate a raw array so as to plot the processed and unprocessed results
    !ALLOCATE(log_Pkraw(nk), log_Pk(nk), log_k(nk))
    ALLOCATE(Pkraw(nk), Pk(nk), k(nk))

    ! Allocate the internal arrays from the cosmology arrays
    k=exp(cosm%log_k_plin)
    Pk=exp(cosm%log_plin)
    !log_k=cosm%log_k_plin
    !log_Pk=cosm%log_plin

    ! Fixed k values - CAMB!
    kv(1)=0.008
    kv(2)=0.01
    kv(3)=0.8
    kv(4)=1.0

    ! Fix p values
    DO i=1,4
       pv(i)=exp(find(log(kv(i)),log(k),log(Pk),nk,3,3,2))
       !pv(i)=exp(find(log(kv(i)),log_k,log_Pk,nk,3,3,2))       
    END DO

    ! Create new 'raw' spectrum which has no wiggles
    DO i=1,nk
       IF(k(i)<=kv(1) .OR. k(i)>=kv(4)) THEN
       !IF(log_k(i)<=log(kv(1)) .OR. log_k(i)>=log(kv(4))) THEN
          Pkraw(i)=Pk(i)
          !log_Pkraw(i)=log_Pk(i)
       ELSE
          Pkraw(i)=exp(Lagrange_polynomial(log(k(i)),3,log(kv),log(pv)))
          !log_Pkraw(i)=Lagrange_polynomial(log_k(i),3,log(kv),log(pv))
       END IF
    END DO

    ! Isolate just the wiggles
    Pk=Pk-Pkraw
    !log_Pk=log(exp(log_Pk)-exp(log_Pkraw))

    ! Damp the wiggles
    Pk=Pk*exp(-(sigv*k)**2)
    !log_Pk=log(exp(log_Pk)*exp(-(sigv*k)**2))
    !log_Pk=log_Pk-(sigv*exp(log_k))**2

    ! Add the damped wiggles back in
    Pk=Pk+Pkraw
    !log_Pk=log(exp(log_Pk)+exp(log_Pkraw))

    ! Create the damped power array
    IF(ALLOCATED(hmod%log_k_pdamp)) DEALLOCATE(hmod%log_k_pdamp)
    IF(ALLOCATED(hmod%log_pdamp))   DEALLOCATE(hmod%log_pdamp)
    ALLOCATE(hmod%log_k_pdamp(nk),hmod%log_pdamp(nk))

    ! Fill the k array
    hmod%log_k_pdamp=log(k)
    !hmod%log_k_pdamp=cosm%log_k_plin
    
    ! Grow damped power to the correct redshift and fill array
    hmod%log_pdamp=log(Pk*grow(a,cosm)**2)
    !hmod%log_pdamp=log(exp(log_Pk)*grow(a,cosm)**2)
    !hmod%log_pdamp=log_Pk+2.*log(grow(a,cosm))

    ! Set the flag
    hmod%has_dewiggle=.TRUE.

  END SUBROUTINE dewiggle_init

  FUNCTION p_dewiggle(k,z,hmod,cosm)

    IMPLICIT NONE
    REAL :: p_dewiggle
    REAL, INTENT(IN) :: k, z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%has_dewiggle .EQV. .FALSE.) CALL dewiggle_init(z,hmod,cosm)
    p_dewiggle=exp(find(log(k),hmod%log_k_pdamp,hmod%log_pdamp,hmod%n_pdamp,3,3,2))

  END FUNCTION p_dewiggle

  FUNCTION halo_type(i)

    !Name halo types
    IMPLICIT NONE
    CHARACTER(len=256) :: halo_type
    INTEGER :: i
    
    halo_type=''
    IF(i==-1) halo_type='DMONLY'
    IF(i==0)  halo_type='Matter'
    IF(i==1)  halo_type='CDM'
    IF(i==2)  halo_type='Gas'
    IF(i==3)  halo_type='Star'
    IF(i==4)  halo_type='Bound gas'
    IF(i==5)  halo_type='Free gas'
    IF(i==6)  halo_type='Electron pressure'
    IF(i==7)  halo_type='Void'
    IF(i==8)  halo_type='Compensated void'
    IF(i==9)  halo_type='Central galaxies'
    IF(i==10) halo_type='Satellite galaxies'
    IF(i==11) halo_type='Galaxies'
    IF(i==12) halo_type='HI'
    IF(halo_type=='') STOP 'HALO_TYPE: Error, i not specified correctly'
    
  END FUNCTION halo_type

  SUBROUTINE set_halo_type(ip)

    !Set the halo types
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

  SUBROUTINE calculate_HMx(itype,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa_full,hmod,cosm,verbose,response)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nk, na, itype(2)
    REAL, INTENT(IN) :: k(:), a(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: powa_2h(:,:), powa_1h(:,:), powa_full(:,:), powa_lin(:,:)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose, response
    REAL, INTENT(IN) :: mmin, mmax
    INTEGER :: i
    REAL :: z
    LOGICAL :: verbose2

    !To avoid splurge of stuff printed to screen
    verbose2=verbose

    IF(ALLOCATED(powa_lin))  DEALLOCATE(powa_lin)
    IF(ALLOCATED(powa_2h))   DEALLOCATE(powa_2h)
    IF(ALLOCATED(powa_1h))   DEALLOCATE(powa_1h)
    IF(ALLOCATED(powa_full)) DEALLOCATE(powa_full)

    !Allocate power arrays
    ALLOCATE(powa_lin(nk,na),powa_2h(nk,na),powa_1h(nk,na),powa_full(nk,na))
    !IF(response) THEN
    !   ALLOCATE(powb_lin(nk,na),powb_2h(nk,na),powb_1h(nk,na),powb_full(nk,na))
    !END IF

    !Do the halo-model calculation
    DO i=na,1,-1
       z=redshift_a(a(i))
       CALL init_halomod(mmin,mmax,z,hmod,cosm,verbose2)
       CALL print_halomod(hmod,cosm,verbose2)
       CALL calculate_halomod(itype(1),itype(2),k,nk,z,powa_lin(:,i),powa_2h(:,i),powa_1h(:,i),powa_full(:,i),hmod,cosm,verbose2,response)
       !IF(response) THEN
       !   CALL calculate_halomod(-1,-1,k,nk,z,powb_lin(:,i),powb_2h(:,i),powb_1h(:,i),powb_full(:,i),hmod,cosm,verbose=.FALSE.)
       !END IF
       IF(i==na .and. verbose) WRITE(*,*) 'CALCULATE_HMx: Doing calculation'       
       IF(verbose) WRITE(*,fmt='(A15,I5,F10.2)') 'CALCULATE_HMx:', i, REAL(z)
       verbose2=.FALSE.
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HMx: Done'
       WRITE(*,*)
    END IF

    !IF(response) THEN
    !   powa_lin=powa_lin/powb_lin
    !   powa_2h=powa_2h/powb_2h
    !   powa_1h=powa_1h/powb_1h
    !   powa_full=powa_full/powb_full
    !END IF

  END SUBROUTINE calculate_HMx

  SUBROUTINE calculate_halomod(itype1,itype2,k,nk,z,pow_lin,pow_2h,pow_1h,pow,hmod,cosm,verbose,response)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: itype1, itype2
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: k(nk), z
    REAL, INTENT(OUT) :: pow_lin(nk), pow_2h(nk), pow_1h(nk), pow(nk)
    REAL :: powg_2h(nk), powg_1h(nk), powg(nk)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod
    LOGICAL, INTENT(IN) :: verbose, response
    INTEGER :: i
    REAL :: plin, a

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'CALCUALTE_HALOMOD: Halo type 1:', itype1
       WRITE(*,*) 'CALCUALTE_HALOMOD: Halo type 2:', itype2
       WRITE(*,*) 'CALCULATE_HALOMOD: k min:', REAL(k(1))
       WRITE(*,*) 'CALCULATE_HALOMOD: k max:', REAL(k(nk))
       WRITE(*,*) 'CALCULATE_HALOMOD: number of k:', nk
       WRITE(*,*) 'CALCULATE_HALOMOD: z:', REAL(z)
       WRITE(*,*) 'CALCULATE_HALOMOD: Calculating halo-model power spectrum'
    END IF

    a=scale_factor_z(z)
    
    ! Loop over k values
    !TODO: add OMP support properly. What is private and what is shared? CHECK THIS!
!!$OMP PARALLEL DO DEFAULT(SHARED)!, private(k,plin,pow_2h,pow_1h,pow,pow_lin)
!!$OMP PARALLEL DO DEFAULT(PRIVATE)
!!$OMP PARALLEL DO FIRSTPRIVATE(nk,cosm,compute_p_lin,k,a,pow_lin,plin,itype1,itype2,z,pow_2h,pow_1h,pow,hmod)
!!$OMP PARALLEL DO
    DO i=1,nk
       
       plin=p_lin(k(i),a,cosm)
       pow_lin(i)=plin

       ! Do the halo model calculation
       CALL calculate_halomod_k(itype1,itype2,k(i),z,pow_2h(i),pow_1h(i),pow(i),plin,hmod,cosm)

       ! If doing a response
       IF(response) THEN
          CALL calculate_halomod_k(-1,-1,k(i),z,powg_2h(i),powg_1h(i),powg(i),plin,hmod,cosm)
          pow_lin(i)=1.
          pow_2h(i)=pow_2h(i)/powg_2h(i)
          pow_1h(i)=pow_1h(i)/powg_1h(i)
          pow(i)=pow(i)/powg(i)
       END IF

    END DO
!!$OMP END PARALLEL DO
    
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_HALOMOD: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_halomod

  SUBROUTINE calculate_halomod_k(ih1,ih2,k,z,p2h,p1h,pfull,plin,hmod,cosm)

    !Gets the one- and two-halo terms and combines them
    IMPLICIT NONE
    REAL, INTENT(OUT) :: p1h, p2h, pfull
    REAL, INTENT(IN) :: plin, k, z
    INTEGER, INTENT(IN) :: ih1, ih2
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: alp, et, nu, a
    REAL :: wk(hmod%n,2), wk2(hmod%n), m, rv, rs
    INTEGER :: i, j, ih(2)
    REAL :: c, dc
    
    !Get the scale factor
    a=scale_factor_z(z)
    
    !Initially fill this small array 
    ih(1)=ih1
    ih(2)=ih2

    !For the i's
    !-1 - DMonly
    ! 0 - All matter
    ! 1 - CDM
    ! 2 - Gas
    ! 3 - Stars
    ! 4 - Bound gas
    ! 5 - Free gas
    ! 6 - Electron pressure
    ! 7 - Voids
    ! 8 - Compensated voids
    ! 9 - Central galaxies
    !10 - Satellite galaxies

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN

       !This should really never be called for k=0
       p1h=0.
       p2h=0.

    ELSE

       !Get eta
       et=eta(a,hmod,cosm)
       
       !No scatter in halo properties
       IF(hmod%iscatter==1) THEN

          !Calculate the halo window functions
          DO j=1,2
             DO i=1,hmod%n
                m=hmod%m(i)
                rv=hmod%rv(i)
                c=hmod%c(i)
                rs=rv/c
                nu=hmod%nu(i)
                wk(i,j)=win_type(.FALSE.,ih(j),1,k*nu**et,z,m,rv,rs,hmod,cosm)
             END DO
             IF(ih(2)==ih(1)) THEN
                !Avoid having to call win_type twice if doing auto spectrum
                wk(:,2)=wk(:,1)
                EXIT
             END IF
          END DO

          !wk(1)*wk(2) in the case of no scatter
          wk2=wk(:,1)*wk(:,2)
          
       ELSE IF(hmod%iscatter==2) THEN

          !Scatter in log concentration: sigma_ln(c)
          dc=0.25

          !Scatter in halo properties
          !TODO: include scatter in two-halo term
          DO i=1,hmod%n
             m=hmod%m(i)
             rv=hmod%rv(i)
             c=hmod%c(i)             
             wk2(i)=integrate_scatter(c,dc,ih,k,z,m,rv,hmod,cosm,hmod%acc_HMx,3)
          END DO
          
       END IF

       !Get the one-halo term
       p1h=p_1h(wk2,k,z,hmod,cosm)

       !If linear theory if used for two-halo term we need to recalculate the window
       !functions for the two-halo term with k=0 fixed
       IF(hmod%ip2h==1 .OR. hmod%smooth_freegas) THEN
          DO j=1,2
             DO i=1,hmod%n
                m=hmod%m(i)
                rv=hmod%rv(i)
                rs=rv/hmod%c(i)
                nu=hmod%nu(i)
                IF(hmod%ip2h==1) wk(i,j)=win_type(.FALSE.,ih(j),2,0.,z,m,rv,rs,hmod,cosm)
                IF(hmod%smooth_freegas) wk(i,j)=win_type(.FALSE.,ih(j),2,k*nu**et,z,m,rv,rs,hmod,cosm)
             END DO
             IF(ih(2)==ih(1)) THEN
                !Avoid having to call win_type twice if doing auto spectrum
                wk(:,2)=wk(:,1)
                EXIT
             END IF
          END DO
       END IF

       !Get the two-halo term
       p2h=p_2h(ih,wk,k,z,plin,hmod,cosm)

    END IF

    ! alpha is set to one sometimes, which is just the standard halo-model sum of terms
    ! No need to have an IF statement around this
    IF(hmod%itrans==2 .OR. hmod%itrans==4) THEN
       IF(p2h<0. .OR. p1h<0.) THEN
          IF(hmod%safe_negative) THEN
             pfull=p2h+p1h
          ELSE
             WRITE(*,*) 'P_2h:', p2h
             WRITE(*,*) 'P_1h:', p1h
             STOP 'CALCULATE_HALOMOD_K: Error, either p2h or p1h is less than zero, this is a problem for the smoothed transition'
          END IF
       ELSE
          alp=alpha_transition(hmod,cosm)    
          pfull=(p2h**alp+p1h**alp)**(1./alp)
       END IF     
    ELSE IF(hmod%itrans==5) THEN
       pfull=p2h+sigmoid_log(1.*k/hmod%knl,1.)*(p1h-p2h)
    ELSE
       pfull=p2h+p1h
    END IF

    !If we are worrying about voids...
    IF(hmod%voids) THEN
       pfull=pfull+p_1void(k,hmod)
    END IF

  END SUBROUTINE calculate_halomod_k
    
  SUBROUTINE halo_diagnostics(z,hmod,cosm,dir)

    !Writes out to file a whole set of halo diagnostics
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL, INTENT(IN) :: z
    CHARACTER(len=*), INTENT(IN) :: dir
    REAL :: mass    
    CHARACTER(len=64) :: ext
    CHARACTER(len=512) :: base
    CHARACTER(len=1024) :: outfile
    INTEGER :: m

    !Integer 10^m to produce haloes between
    INTEGER, PARAMETER :: m1=10
    INTEGER, PARAMETER :: m2=16

    WRITE(*,*) 'HALO_DIAGNOSTICS: Outputting diagnostics'

    outfile=TRIM(dir)//'/mass_fractions.dat'
    WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
    CALL write_mass_fractions(hmod,cosm,outfile)

    IF(z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_DIAGNOSTICS: Error, need to make this better with z'
    END IF

    DO m=m1,m2

       mass=10.**m

       base=TRIM(dir)//'/halo_profile_m'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_profiles(mass,z,hmod,cosm,outfile)

       base=TRIM(dir)//'/halo_window_m'
       outfile=number_file(base,m,ext)
       WRITE(*,*) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
       CALL write_halo_transforms(mass,z,hmod,cosm,outfile)

    END DO

    WRITE(*,*) 'HALO_DIAGNOSTICS: Done'
    WRITE(*,*)

  END SUBROUTINE halo_diagnostics

  SUBROUTINE halo_definitions(z,hmod,cosm,dir)

    !Writes out to files the different halo definitions
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: dir
    CHARACTER(len=256) :: fradius, fmass, fconc!, ext
    CHARACTER(len=64) :: ext
    INTEGER :: i

    WRITE(*,*) 'HALO_DEFINITIONS: Outputting definitions'

    IF(z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(z==2.0) THEN
       ext='_z2.0.dat'
    ELSE
       STOP 'HALO_DIAGNOSTICS: Error, need to make this better with z'
    END IF

    fradius=TRIM(dir)//'/radius'//TRIM(ext)
    fmass=TRIM(dir)//'/mass'//TRIM(ext)
    fconc=TRIM(dir)//'/concentration'//TRIM(ext)

    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fradius)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fmass)
    WRITE(*,*) 'HALO_DEFINITIONS: ', TRIM(fconc)
    
    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(z,hmod,cosm)

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

  SUBROUTINE halo_properties(z,hmod,dir)

    !Writes out to files the different halo definitions
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    CHARACTER(len=*), INTENT(IN) :: dir
    CHARACTER(len=256) :: output
    CHARACTER(len=64) :: ext
    INTEGER :: i

    WRITE(*,*) 'HALO_PROPERTIES: Outputting definitions'

    IF(z==0.0) THEN
       ext='_z0.0.dat'
    ELSE IF(z==0.5) THEN
       ext='_z0.5.dat'
    ELSE IF(z==1.0) THEN
       ext='_z1.0.dat'
    ELSE IF(z==2.0) THEN
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

    !Writes out the halo mass fractions
    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile   
    REAL :: m, mmin, mmax
    INTEGER :: i, j, n

    mmin=1e10
    mmax=1e16
    n=101

    OPEN(7,file=outfile)
    DO i=1,n
       m=exp(progression(log(mmin),log(mmax),i,n))
       WRITE(7,*) m, (halo_fraction(j,m,hmod,cosm), j=1,5)
    END DO
    CLOSE(7)

  END SUBROUTINE write_mass_fractions

  SUBROUTINE write_halo_profiles(m,z,hmod,cosm,outfile)

    !Writes out the halo density profiles
    IMPLICIT NONE
    REAL, INTENT(IN) :: m, z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=*), INTENT(IN) :: outfile    
    REAL :: r, rv, rs, c
    INTEGER :: i, j   

    REAL, PARAMETER :: rmin=1e-3 !Mininum r/rv
    REAL, PARAMETER :: rmax=1.1e0 !Maximum r/rv
    INTEGER, PARAMETER :: n=512 !Number of points
    LOGICAL, PARAMETER :: rsp=.TRUE. !Real profiles

    !Calculate halo attributes
    rv=exp(find(log(m),hmod%log_m,log(hmod%rv),hmod%n,3,3,2))
    c=find(log(m),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c

    !Write file
    OPEN(7,file=outfile)
    DO i=1,n
       r=exp(progression(log(rmin),log(rmax),i,n))
       r=r*rv
       WRITE(7,*) r/rv, (win_type(rsp,j,1,r,z,m,rv,rs,hmod,cosm)*rv**3, j=1,6) !rv**3 here is from r^2 dr in integral
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_profiles

  SUBROUTINE write_halo_transforms(m,z,hmod,cosm,outfile)

    !Writes out to file the Fourier transform of the halo density profiles
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: outfile
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, INTENT(IN) :: m, z
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: x, rv, c, rs, k, rhobar
    INTEGER :: i, j

    REAL, PARAMETER :: xmin=1e-1 !Minimum r/rv
    REAL, PARAMETER :: xmax=1e2 !Maximum r/rv
    INTEGER, PARAMETER :: n=512 !Number of points
    LOGICAL, PARAMETER :: rsp=.FALSE. !Fourier profiles

    !Calculate halo attributes
    rv=exp(find(log(m),hmod%log_m,log(hmod%rv),hmod%n,3,3,2))
    c=find(log(m),hmod%log_m,hmod%c,hmod%n,3,3,2)
    rs=rv/c

    !Need mean density
    rhobar=comoving_matter_density(cosm)

    !Write file
    OPEN(7,file=outfile)
    DO i=1,n
       x=exp(progression(log(xmin),log(xmax),i,n))
       k=x/rv
       WRITE(7,*) x, (win_type(rsp,j,1,k,z,m,rv,rs,hmod,cosm)*rhobar/m, j=1,6)
    END DO
    CLOSE(7)

  END SUBROUTINE write_halo_transforms

  FUNCTION delta_c(a,hmod,cosm)

    !Linear collapse density
    IMPLICIT NONE
    REAL :: delta_c
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%idc==1) THEN
       !Fixed value
       delta_c=1.686
    ELSE IF(hmod%idc==2) THEN
       !From Nakamura & Suto (1997) LCDM fitting function
       delta_c=dc_NakamuraSuto(a,cosm)
    ELSE IF(hmod%idc==3) THEN
       !From Mead et al. (2015, 2016)
       !delta_c=1.59+0.0314*log(sigma(8.,a,cosm))
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

  FUNCTION Delta_v(a,hmod,cosm)

    !Virialised overdensity
    IMPLICIT NONE
    REAL :: Delta_v
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%iDv==1) THEN
       !Fixed value
       Delta_v=200.
    ELSE IF(hmod%iDv==2) THEN
       !From Bryan & Norman (1998; arXiv:astro-ph/9710107) fitting functions
       Delta_v=Dv_BryanNorman(a,cosm)    
    ELSE IF(hmod%iDv==3) THEN
       !From Mead et al. (2015, 2016)
       !Delta_v=418.*Omega_m(a,cosm)**(-0.352)
       Delta_v=hmod%Dv0*Omega_m(a,cosm)**hmod%Dv1
    ELSE IF(hmod%iDv==4) THEN
       !From Mead (2017) fitting function
       Delta_v=Dv_Mead(a,cosm)
    ELSE IF(hmod%iDv==5) THEN
       !From spheircal-collapse calculation
       Delta_v=Dv_spherical(a,cosm)
    ELSE
       STOP 'DELTA_V: Error, iDv defined incorrectly'
    END IF

  END FUNCTION Delta_v

  FUNCTION eta(a,hmod,cosm)

    !Calculates the eta that comes into the bastardised one-halo term
    IMPLICIT NONE
    REAL :: eta
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(hmod%ieta==1) THEN
       eta=0.
    ELSE IF(hmod%ieta==2) THEN
       !From Mead et al. (2015; arXiv 1505.07833, 2016)
       !eta=0.603-0.3*(sigma(8.,a,cosm))
       eta=hmod%eta0-hmod%eta1*(sigma(8.,a,cosm))
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION eta

  FUNCTION kstar(hmod,cosm)

    !Calculates the one-halo damping wave number
    IMPLICIT NONE
    REAL :: kstar
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: crap

    !To prevent compile-time warnings
    crap=cosm%A

    IF(hmod%ikstar==1) THEN
       !Set to zero for the standard Poisson one-halo term
       kstar=0.
    ELSE IF(hmod%ikstar==2) THEN
       !One-halo cut-off wavenumber from Mead et al. (2015, 2016)
       kstar=hmod%ks/hmod%sigv
    ELSE
       STOP 'KSTAR: Error, ihm defined incorrectly'
    END IF

  END FUNCTION kstar

  FUNCTION As(hmod,cosm)

    !Halo concentration pre-factor
    IMPLICIT NONE
    REAL :: As
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    !To prevent compile-time warnings
    crap=cosm%A
    
    IF(hmod%iAs==1) THEN
       !Set to 4 for the standard Bullock value
       As=4.
    ELSE IF(hmod%iAs==2) THEN
       !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833, 2016)
       !As=3.13
       As=hmod%A
    ELSE
       STOP 'AS: Error, iconc defined incorrectly'
    END IF

    ! Now this is divided by 4 so as to be relative to the Bullock base result
    As=As/4.

  END FUNCTION As

  FUNCTION fdamp(a,hmod,cosm)

    !Calculates the linear-theory damping factor
    IMPLICIT NONE
    REAL ::fdamp
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: crap

    REAL, PARAMETER :: fdamp_min=1e-3
    REAL, PARAMETER :: fdamp_max=0.99

    !To prevent compile-time warnings
    crap=cosm%A
    crap=a

    IF(hmod%i2hdamp==1) THEN
       !Set to 0 for the standard linear theory two halo term
       fdamp=0.
    ELSE IF(hmod%i2hdamp==2) THEN
       !Mead et al. (2015)
       !fdamp=0.188*hmod%sig8z**4.29
       fdamp=hmod%f0*hmod%sig8z**hmod%f1
    ELSE IF(hmod%i2hdamp==3) THEN
       !Mead et al. (2016)
       !fdamp=0.0095*hmod%sigv100**1.37
       fdamp=hmod%f0*hmod%sigv100**hmod%f1
    ELSE
       STOP 'FDAMP: Error, i2hdamp defined incorrectly'
    END IF

    !Catches extreme values of fdamp that occur for ridiculous cosmologies
    IF(fdamp<fdamp_min) fdamp=0.
    IF(fdamp>fdamp_max) fdamp=fdamp_max

  END FUNCTION fdamp

  FUNCTION alpha_transition(hmod,cosm)

    !Calculates the alpha to smooth the two- to one-halo transition
    IMPLICIT NONE
    REAL :: alpha_transition
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    REAL, PARAMETER :: alpha_min=0.5
    REAL, PARAMETER :: alpha_max=2.0

    !To prevent compile-time warnings
    crap=cosm%A

!!$    IF(hmod%itrans==2) THEN
!!$       !From Mead et al. (2015)   
!!$       !This uses the top-hat defined neff in contrast to the neff in HALOFIT
!!$       !alpha_transition=2.93*1.77**hmod%neff
!!$       alpha_transition=hmod%alp0*hmod%alp1**hmod%neff 
!!$    ELSE IF(hmod%itrans==3) THEN
!!$       !From Mead et al. (2016)
!!$       !This uses the top-hat defined neff in contrast to the neff in HALOFIT
!!$       !alpha_transition=3.24*1.85**hmod%neff
!!$       alpha_transition=hmod%alp0*hmod%alp1**hmod%neff
    IF(hmod%itrans==2) THEN
       !From Mead et al. (2015, 2016)   
       alpha_transition=hmod%alp0*hmod%alp1**hmod%neff          
    ELSE IF(hmod%itrans==4) THEN
       !Specially for HMx, exponentiated Mead et al. (2016) result
       !alpha_transition=(3.24*1.85**hmod%neff)**2.5
       alpha_transition=(hmod%alp0*hmod%alp1**hmod%neff)**2.5
    ELSE
       alpha_transition=1.
    END IF

    !Catches values of alpha that are crazy
    IF(alpha_transition<alpha_min) alpha_transition=alpha_min
    IF(alpha_transition>alpha_max) alpha_transition=alpha_max 

  END FUNCTION alpha_transition

  SUBROUTINE print_halomod(hmod,cosm,verbose)

    !This subroutine writes out the physical halo-model parameters at some redshift 
    !(e.g., Delta_v) rather than the model parameters
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod
    LOGICAL, INTENT(IN) :: verbose
    REAL :: a

    IF(verbose) THEN

       a=scale_factor_z(hmod%z)

       WRITE(*,*) 'HALOMODEL: Writing out halo-model parameters'
       WRITE(*,*) '============================================'
       WRITE(*,*) TRIM(hmod%name)

       WRITE(*,*) '============================================'
       WRITE(*,*) 'HALOMODEL: redshift:', hmod%z
       WRITE(*,*) '============================================'
       WRITE(*,*) 'HALOMODEL: Number of points in look-up tables:', hmod%n
       WRITE(*,*) 'HALOMODEL: Halo model accuracy parameter:', hmod%acc_HMx
       WRITE(*,*) 'HALOMODEL: Large value of nu:', hmod%large_nu
       WRITE(*,*) '============================================'

       ! Form of the two-halo term
       IF(hmod%ip2h==1) WRITE(*,*) 'HALOMODEL: Linear two-halo term'
       IF(hmod%ip2h==2) WRITE(*,*) 'HALOMODEL: Standard two-halo term (Seljak 2000)'
       IF(hmod%ip2h==3) WRITE(*,*) 'HALOMODEL: Linear two-halo term with damped wiggles'

       ! Order to go to in halo bias
       IF(hmod%ip2h .NE. 1) THEN
          IF(hmod%ibias==1) WRITE(*,*) 'HALOMODEL: Linear halo bias'
          IF(hmod%ibias==2) WRITE(*,*) 'HALOMODEL: Second-order halo bias'
       END IF

       ! Correction for missing low-mass haloes
       IF(hmod%ip2h .NE. 1) THEN
          IF(hmod%ip2h_corr==1) WRITE(*,*) 'HALOMODEL: No two-halo correction applied for missing low-mass haloes'
          IF(hmod%ip2h_corr==2) WRITE(*,*) 'HALOMODEL: Two-halo term corrected by adding missing g(nu)b(nu)' 
          IF(hmod%ip2h_corr==3) WRITE(*,*) 'HALOMODEL: Two-halo term corrected via delta function at low mass end'      
       END IF

       ! Scale-dependent halo bias
       IF(hmod%ikb==1) WRITE(*,*) 'HALOMODEL: Scale-independent halo bias'
       IF(hmod%ikb==2) WRITE(*,*) 'HALOMODEL: Scale-dependent halo bias fudge from Fedeli (2014b)'
       IF(hmod%ikb==3) WRITE(*,*) 'HALOMODEL: Scale-dependent halo bias fudge'

       !Halo mass function
       IF(hmod%imf==1) WRITE(*,*) 'HALOMODEL: Press & Schecter (1974) mass function'
       IF(hmod%imf==2) WRITE(*,*) 'HALOMODEL: Sheth & Tormen (1999) mass function'
       IF(hmod%imf==3) WRITE(*,*) 'HALOMODEL: Tinker et al. (2010) mass function'

       !Concentration-mass relation
       IF(hmod%iconc==1) WRITE(*,*) 'HALOMODEL: Full Bullock et al. (2001) concentration-mass relation'
       IF(hmod%iconc==2) WRITE(*,*) 'HALOMODEL: Simple Bullock et al. (2001) concentration-mass relation'
       IF(hmod%iconc==3) WRITE(*,*) 'HALOMODEL: Full-sample mean density Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==4) WRITE(*,*) 'HALOMODEL: Full-sample virial denity Duffy et al. (2008) concentration-mass relation'
       IF(hmod%iconc==5) WRITE(*,*) 'HALOMODEL: Relaxed-sample x200 critical density Duffy et al. (2008) concentration-mass relation'

       !Concentration-mass relation correction
       IF(hmod%iDolag==1) WRITE(*,*) 'HALOMODEL: No concentration-mass correction for dark energy'
       IF(hmod%iDolag==2) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction'
       IF(hmod%iDolag==3) WRITE(*,*) 'HALOMODEL: Dolag (2004) dark energy halo concentration correction with 1.5 exponent'

       !delta_c
       IF(hmod%idc==1) WRITE(*,*) 'HALOMODEL: Fixed delta_c = 1.686'
       IF(hmod%idc==2) WRITE(*,*) 'HALOMODEL: delta_c from Nakamura & Suto (1997) fitting function'
       IF(hmod%idc==3) WRITE(*,*) 'HALOMODEL: delta_c from Mead et al. (2015, 2016) power spectrum fit'
       IF(hmod%idc==4) WRITE(*,*) 'HALOMODEL: delta_c from Mead (2017) fitting function'
       IF(hmod%idc==5) WRITE(*,*) 'HALOMODEL: delta_c from spherical-collapse calculation'

       !Delta_v
       IF(hmod%iDv==1) WRITE(*,*) 'HALOMODEL: Fixed Delta_v = 200'
       IF(hmod%iDv==2) WRITE(*,*) 'HALOMODEL: Delta_v from Bryan & Norman (1998) fitting function'
       IF(hmod%iDv==3) WRITE(*,*) 'HALOMODEL: Delta_v from Mead et al. (2015, 2016) power spectrum fit'
       IF(hmod%iDv==4) WRITE(*,*) 'HALOMODEL: Delta_v from Mead (2017) fitting function'
       IF(hmod%iDv==5) WRITE(*,*) 'HALOMODEL: Delta_v from spherical-collapse calculation'

       !Eta for halo window function
       IF(hmod%ieta==1) WRITE(*,*) 'HALOMODEL: eta = 0 fixed'
       IF(hmod%ieta==2) WRITE(*,*) 'HALOMODEL: eta from Mead et al. (2015, 2016) power spectrum fit'

       !Small-scale two-halo term damping coefficient
       IF(hmod%i2hdamp==1) WRITE(*,*) 'HALOMODEL: No two-halo term damping at small scales'
       IF(hmod%i2hdamp==2) WRITE(*,*) 'HALOMODEL: Two-halo term damping from Mead et al. (2015)'
       IF(hmod%i2hdamp==3) WRITE(*,*) 'HALOMODEL: Two-halo term damping from Mead et al. (2016)'

       !Large-scale one-halo term damping function
       IF(hmod%i1hdamp==1) WRITE(*,*) 'HALOMODEL: No damping in one-halo term at large scales'
       IF(hmod%i1hdamp==2) WRITE(*,*) 'HALOMODEL: One-halo term large-scale damping via an exponential'
       IF(hmod%i1hdamp==3) WRITE(*,*) 'HALOMODEL: One-halo term large-scale damping like Delta^2 ~ k^7'

       !Large-scale one-halo term damping coefficient
       IF(hmod%i1hdamp .NE. 1) THEN
          IF(hmod%ikstar==1) WRITE(*,*) 'HALOMODEL: No damping in one-halo term at large scales'
          IF(hmod%ikstar==2) WRITE(*,*) 'HALOMODEL: One-halo term damping function from Mead et al. (2015, 2016)'
       END IF

       !Concentration-mass scaling
       IF(hmod%iAs==1) WRITE(*,*) 'HALOMODEL: No rescaling of concentration-mass relation'
       IF(hmod%iAs==2) WRITE(*,*) 'HALOMODEL: Concentration-mass relation rescaled mass independetly (Mead et al. 2015, 2016)'

       !Scatter in halo properties
       IF(hmod%iscatter==1) WRITE(*,*) 'HALOMODEL: No scatter in halo properties at fixed mass'
       IF(hmod%iscatter==2) WRITE(*,*) 'HALOMODEL: Scatter in halo concentration at fixed mass'

       !Two- to one-halo transition region
       IF(hmod%itrans==1) WRITE(*,*) 'HALOMODEL: Standard sum of two- and one-halo terms'
       IF(hmod%itrans==2) WRITE(*,*) 'HALOMODEL: Smoothed transition using alpha'
       IF(hmod%itrans==4) WRITE(*,*) 'HALOMODEL: Experimental smoothed transition for HMx'
       IF(hmod%itrans==5) WRITE(*,*) 'HALOMODEL: Tanh transition with k_nl'

       IF(hmod%use_UPP) WRITE(*,*) 'HALOMODEL: Using UPP for all electron pressure calculations'

       !Numerical parameters
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMcode parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'Dv:', Delta_v(a,hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'dc:', delta_c(a,hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'eta:', eta(a,hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'k*:', kstar(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'A:', As(hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'fdamp:', fdamp(a,hmod,cosm)
       WRITE(*,fmt='(A30,F10.5)') 'alpha:', alpha_transition(hmod,cosm)
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HMx parameters'
       WRITE(*,*) '======================================='
       IF(.NOT. hmod%fixed_HMx) WRITE(*,fmt='(A30,F10.5)') 'log10(T_heat) [K]:', log10(hmod%Theat)
       WRITE(*,fmt='(A30,F10.5)') 'alpha:', HMx_alpha(hmod)
       WRITE(*,fmt='(A30,F10.5)') 'epsilon:', HMx_eps(hmod)
       WRITE(*,fmt='(A30,F10.5)') 'Gamma:', HMx_Gamma(hmod)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M0) [Msun/h]:', log10(HMx_M0(hmod))
       WRITE(*,fmt='(A30,F10.5)') 'A*:', HMx_Astar(hmod)
       WRITE(*,fmt='(A30,F10.5)') 'log10(T_WHIM) [K]:', log10(HMx_Twhim(hmod))
       WRITE(*,fmt='(A30,F10.5)') 'R*:', hmod%rstar
       WRITE(*,fmt='(A30,F10.5)') 'sigma*:', hmod%sstar
       WRITE(*,fmt='(A30,F10.5)') 'log10(M*) [Msun/h]:', log10(hmod%Mstar)
       WRITE(*,*) '======================================='
       WRITE(*,*) 'HALOMODEL: HOD/HI parameters'
       WRITE(*,*) '======================================='
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_gal) [Msun/h]:', log10(hmod%mgal)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_HI_min) [Msun/h]:', log10(hmod%HImin)
       WRITE(*,fmt='(A30,F10.5)') 'log10(M_HI_max) [Msun/h]:', log10(hmod%HImax)
       WRITE(*,*) '======================================='
       WRITE(*,*)
       
    END IF

  END SUBROUTINE print_halomod

  REAL FUNCTION HMx_alpha(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%fixed_HMx) THEN

       HMx_alpha=hmod%alpha

    ELSE

       A=hmod%A_alpha
       B=hmod%B_alpha
       C=hmod%C_alpha
       D=hmod%D_alpha

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_alpha=A*(1.+z)*T+B*(1.+z)+C*T+D

    END IF
       
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

  REAL FUNCTION HMx_Gamma(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%fixed_HMx) THEN

       HMx_Gamma=hmod%Gamma

    ELSE

       A=hmod%A_Gamma
       B=hmod%B_Gamma
       C=hmod%C_Gamma
       D=hmod%D_Gamma

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_Gamma=A*(1.+z)*T+B*(1.+z)+C*T+D

    END IF
       
  END FUNCTION HMx_Gamma

  REAL FUNCTION HMx_M0(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%fixed_HMx) THEN

       HMx_M0=hmod%M0

    ELSE

       A=hmod%A_M0
       B=hmod%B_M0
       C=hmod%C_M0
       D=hmod%D_M0

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_M0=A*(1.+z)*T+B*(1.+z)+C*T+D
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
       
  END FUNCTION HMx_Astar

  REAL FUNCTION HMx_Twhim(hmod)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: z, T, A, B, C, D

    IF(hmod%fixed_HMx) THEN

       HMx_Twhim=hmod%Twhim

    ELSE

       A=hmod%A_Twhim
       B=hmod%B_Twhim
       C=hmod%C_Twhim
       D=hmod%D_Twhim

       z=hmod%z
       T=log10(hmod%Theat)
       HMx_Twhim=A*(1.+z)*T+B*(1.+z)+C*T+D
       HMx_Twhim=10**HMx_Twhim

    END IF
       
  END FUNCTION HMx_Twhim

  FUNCTION r_nl(hmod)

    !Calculates R_nl where nu(R_nl)=1.
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL :: r_nl  

    IF(hmod%nu(1)>1.) THEN
       !This catches some very strange values
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

    ! Experimental window look-up table
    !hmod%nk=nk
    !ALLOCATE(hmod%log_m(n),hmod%log_k(nk),hmod%log_win(n,nk))
    !hmod%log_k=0.
    !hmod%log_win=0.
    !hmod%iwin=.FALSE.

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

  SUBROUTINE assign_halomod(ihm,hmod,verbose)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ihm
    TYPE(halomod), INTENT(OUT) :: hmod
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    !Names of pre-defined halo models
    INTEGER, PARAMETER :: nhalomod=17 !Number of pre-defined halo-model types
    CHARACTER(len=256):: names(1:nhalomod)    
    names(1)='Accurate halo-model calculation (Mead et al. 2016)'
    names(2)='Basic halo-model calculation (Two-halo term is linear)'
    names(3)='Standard halo-model calculation (Seljak 2000)'
    names(4)='Standard halo-model calculation but with Mead et al. (2015) transition'
    names(5)='Standard halo-model calculation but with Delta_v=200 and delta_c=1.686 and Bullock c(M)'
    names(6)='Half-accurate halo-model calculation (Mead et al. 2015, 2016)'
    names(7)='Accurate halo-model calculation (Mead et al. 2015)'
    names(8)='Including scatter in halo properties at fixed mass'
    names(9)='Parameters for CCL tests'
    names(10)='Comparison of mass conversions with Wayne Hu code'
    names(11)='Standard halo-model calculation (Seljak 2000) but with UPP for electron pressure'
    names(12)='Spherical collapse used for Mead (2017) results'
    names(13)='Experimental log-tanh transition'
    names(14)='Experimental scale-dependent halo bias'
    names(15)='Accurate halo-model calculation (Mead et al. 2018) ...'
    names(16)='Halo-void model'
    names(17)='HMx'

    IF(verbose) WRITE(*,*) 'ASSIGN_HALOMOD: Assigning halo model'
    
    ! Default options

    ! Number of points in integration (128 is okay, 1024 is better)
    hmod%n=128

    ! Accuracy for continuous integrals (1e-3 is okay, 1e-4 is better
    hmod%acc_HMx=1e-3

    ! A large value for nu (6 is okay, corrections are suppressed by exp(-large_nu^2)
    hmod%large_nu=6.

    ! Two-halo term
    ! 1 - Linear theory
    ! 2 - Standard from Seljak (2000)
    ! 3 - Linear theory with damped wiggles
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

    ! Scatter in halo properties at fixed mass
    ! 1 - No
    ! 2 - Scatter in halo concentration 
    hmod%iscatter=1

    ! Do voids?
    hmod%voids=.FALSE.

    ! Use UPP for electron pressure?
    hmod%use_UPP=.FALSE.

    ! Smoothly distribute free gas?
    hmod%smooth_freegas=.TRUE.

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
    hmod%A=3.13
    hmod%alp0=3.24
    hmod%alp1=1.85

    ! HMx parameters
    hmod%fixed_HMx=.TRUE.
    hmod%alpha=1.
    hmod%eps=1.
    hmod%Gamma=1.17
    hmod%M0=1e14 ! Halo mass that has lost half gas
    hmod%Astar=0.02 ! Maximum star-formation efficiency
    hmod%Twhim=1e6 ! WHIM temperature [K]
    hmod%rstar=0.1
    hmod%sstar=1.2
    hmod%Mstar=5e12

    !$\alpha$ z and Theat variation
    hmod%A_alpha=0.862
    hmod%B_alpha=-6.127
    hmod%C_alpha=-0.234
    hmod%D_alpha=1.670

    !$\log_{10} \epsilon_c$ z and Theat variation
    hmod%A_eps=0.074
    hmod%B_eps=-0.343
    hmod%C_eps=0.827
    hmod%D_eps=-5.929

    !$\Gamma$ z and Theat variation
    hmod%A_Gamma=0.600
    hmod%B_Gamma=-4.391
    hmod%C_Gamma=-0.302
    hmod%D_Gamma=3.380

    !$\log_{10}M_0$ z and Theat variation
    hmod%A_M0=0.473
    hmod%B_M0=-4.172
    hmod%C_M0=1.352
    hmod%D_M0=3.738

    !$A_*$ z and Theat variation
    hmod%A_Astar=0.006
    hmod%B_Astar=-0.057
    hmod%C_Astar=-0.020
    hmod%D_Astar=0.197

    !$\log_{10}T_{WHIM}$ z and Theat variation
    hmod%A_Twhim=0.564
    hmod%B_Twhim=-4.676
    hmod%C_Twhim=-0.583
    hmod%D_Twhim=10.871

    ! Default values of the HOD parameters
    hmod%mgal=1e13

    ! Default values for the HI parameters
    hmod%HImin=1e9
    hmod%HImax=1e12
    
    IF(ihm==-1) THEN
       WRITE(*,*) 'ASSIGN_HALOMOD: Choose your halo model'
       DO i=1,nhalomod
          WRITE(*,*) i, TRIM(names(i))
       END DO
       READ(*,*) ihm
       WRITE(*,*)
    END IF
       
    IF(ihm==1 .OR. ihm==7 .OR. ihm==15) THEN
       ! 1 - Accurate HMcode halo-model calculation (Mead et al. 2016)
       ! 7 - Accurate HMcode halo-model calculation (Mead et al. 2015)
       !15 - Accurate HMcode halo-model calculation (Mead et al. 2018)
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
          hmod%ip2h=3 ! Linear theory with damped wiggles
       END IF       
    ELSE IF(ihm==2) THEN
       ! 2 - Basic halo model with linear two halo term (Delta_v = 200, delta_c = 1.686))
       hmod%ip2h=1
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
    ELSE IF(ihm==3) THEN
       !3 - Standard halo-model calculation (Seljak 2000)
       !This is the default, so do nothing here
    ELSE IF(ihm==4) THEN
       !4 - Standard halo-model calculation but with Mead et al. (2015) smoothed two- to one-halo transition and one-halo damping       
       hmod%itrans=4
       hmod%ikstar=2
       hmod%i1hdamp=3
       hmod%safe_negative=.TRUE.
    ELSE IF(ihm==5) THEN
       !5 - Standard halo-model calculation but with Delta_v = 200 and delta_c = 1.686 fixed and Bullock c(M)
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
    ELSE IF(ihm==6) THEN
       !6 - Half-accurate halo-model calculation (Mead et al. 2015, 2016)
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
       !hmod%ip2h=1
       !hmod%i1hdamp=2
       !hmod%iconc=1
       !hmod%idc=3
       !hmod%iDv=3
       !hmod%ieta=2
       !hmod%ikstar=2
       !hmod%iAs=2
       !hmod%i2hdamp=3
       !hmod%itrans=2
       !hmod%iDolag=3
    ELSE IF(ihm==8) THEN
       !8 - Include scatter in halo properties
       hmod%idc=1
       hmod%iDv=1
       hmod%iconc=1
       hmod%iscatter=2
    ELSE IF(ihm==9) THEN
       !9 - For CCL comparison
       hmod%n=2048 ! Increase accuracy for the CCL benchmarks
       hmod%acc_HMx=1e-5 ! Increase accuracy for the CCL benchmarks
       hmod%ip2h=2
       hmod%ip2h_corr=3
       hmod%ibias=1
       hmod%i1hdamp=1
       hmod%imf=2
       hmod%iconc=4 ! Virial Duffy relation
       hmod%idc=2 ! Virial dc
       hmod%iDv=2 ! Virial Dv
       hmod%ieta=1
       hmod%ikstar=1
       hmod%iAs=1
       hmod%i2hdamp=1
       hmod%itrans=1
       hmod%iDolag=1
       hmod%iscatter=1
       hmod%use_UPP=.FALSE.
       hmod%smooth_freegas=.FALSE.
    ELSE IF(ihm==10) THEN
       !10 - For mass conversions comparison with Wayne Hu's code
       hmod%iconc=2
       hmod%idc=1
       hmod%iDv=1
    ELSE IF(ihm==11) THEN
       !11 - UPP
       hmod%use_UPP=.TRUE.
    ELSE IF(ihm==12) THEN
       !12 - Spherical-collapse model to produce Mead (2017) results
       hmod%iconc=1
       hmod%idc=5
       hmod%iDv=5
       hmod%iDolag=2 ! This seems not to be important for these results
    ELSE IF(ihm==13) THEN
       !13 - Experimental log-tanh transition
       hmod%itrans=5
    ELSE IF(ihm==14) THEN
       !14 - Experimental scale-dependent halo bias
       hmod%ikb=3
       hmod%ikstar=2
       hmod%i1hdamp=3
    ELSE IF(ihm==16) THEN
       !16 - Halo-void model
       hmod%voids=.TRUE.
    ELSE IF(ihm==17) THEN
       !17 - HMx
       hmod%fixed_HMx=.FALSE.
       hmod%Theat=10**7.8
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

  SUBROUTINE init_halomod(mmin,mmax,z,hmod,cosm,verbose)

    ! Halo-model initialisation routine
    ! The computes other tables necessary for the one-halo integral
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    REAL, INTENT(IN) :: mmin, mmax
    LOGICAL, INTENT(IN) :: verbose
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: Dv, dc, m, nu, R, sig, A0, frac, a

    LOGICAL, PARAMETER :: slow=.FALSE.

    ! Set the redshift (this routine needs to be called anew for each z)
    hmod%z=z

    IF(ALLOCATED(hmod%log_m)) CALL deallocate_HMOD(hmod)
    CALL allocate_HMOD(hmod)

    ! Set flags to false
    hmod%has_galaxies=.FALSE.
    hmod%has_HI=.FALSE.
    hmod%has_mass_conversions=.FALSE.
    hmod%has_dewiggle=.FALSE.

    ! Get the scale factor
    a=scale_factor_z(z)

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

    DO i=1,hmod%n

       m=exp(progression(log(mmin),log(mmax),i,hmod%n))
       R=radius_m(m,cosm)
       sig=sigma(R,a,cosm)
       nu=nu_R(R,a,hmod,cosm)

       hmod%m(i)=m
       hmod%rr(i)=R
       hmod%sig(i)=sig
       hmod%nu(i)=nu

    END DO

    hmod%log_m=log(hmod%m)

    IF(verbose) WRITE(*,*) 'INIT_HALOMOD: M, R, nu, sigma tables filled'

    ! Get delta_c
    dc=delta_c(a,hmod,cosm)

    ! Fill virial radius table using real radius table
    Dv=Delta_v(a,hmod,cosm)
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
       WRITE(*,*) 'INIT_HALOMOD: Minimum log10(M/[Msun/h]):', REAL(log10(hmod%m(1)))
       WRITE(*,*) 'INIT_HALOMOD: Maximum log10(M/[Msun/h]):', REAL(log10(hmod%m(hmod%n)))
    END IF

    ! Calculate missing mass things if necessary
    IF(hmod%ip2h_corr==2 .OR. hmod%ip2h_corr==3) THEN 
       IF(slow) hmod%gmin=1.-integrate_hmod(hmod%nu(1),hmod%large_nu,g_nu,hmod,hmod%acc_HMx,3)
       IF(slow) hmod%gmax=integrate_hmod(hmod%nu(hmod%n),hmod%large_nu,g_nu,hmod,hmod%acc_HMx,3)
       hmod%gbmin=1.-integrate_hmod(hmod%nu(1),hmod%large_nu,gb_nu,hmod,hmod%acc_HMx,3)
       IF(slow) hmod%gbmax=integrate_hmod(hmod%nu(hmod%n),hmod%large_nu,gb_nu,hmod,hmod%acc_HMx,3)
       IF(verbose) THEN
          IF(slow) WRITE(*,*) 'INIT_HALOMOD: Missing g(nu) at low end:', REAL(hmod%gmin)
          IF(slow) WRITE(*,*) 'INIT_HALOMOD: Missing g(nu) at high end:', REAL(hmod%gmax)
          WRITE(*,*) 'INIT_HALOMOD: Missing g(nu)b(nu) at low end:', REAL(hmod%gbmin)
          IF(slow) WRITE(*,*) 'INIT_HALOMOD: Missing g(nu)b(nu) at high end:', REAL(hmod%gbmax)
       END IF
    END IF
       
    ! Calculate the total stellar mass fraction
    IF(slow) THEN
       !ALLOCATE(integrand(hmod%n))
       !DO i=1,hmod%n
       !   integrand(i)=halo_fraction(3,hmod%m(i),cosm)*g_nu(hmod%nu(i),hmod)
       !END DO
       !frac=integrate_table(hmod%nu,integrand,hmod%n,1,hmod%n,3)
       !frac=total_stellar_mass_fraction(hmod,cosm)
       !WRITE(*,*) 'INIT_HALOMOD: Total stellar mass fraction:', frac
       frac=total_stellar_mass_fraction(hmod,cosm)
       IF(verbose) WRITE(*,*) 'INIT_HALOMOD: Total stellar mass fraction:', frac
    END IF

    ! Find non-linear radius and scale
    ! This is defined as nu(M_star)=1 *not* sigma(M_star)=1, so depends on delta_c
    hmod%rnl=r_nl(hmod)
    hmod%mnl=mass_r(hmod%rnl,cosm)
    hmod%knl=1./hmod%rnl

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Non-linear mass [log10(M*/[Msun/h])]:', REAL(log10(hmod%mnl))
       WRITE(*,*) 'INIT_HALOMOD: Non-linear halo virial radius [Mpc/h]:', REAL(virial_radius(hmod%mnl,a,hmod,cosm))
       WRITE(*,*) 'INIT_HALOMOD: Non-linear Lagrangian radius [Mpc/h]:', REAL(hmod%rnl)
       WRITE(*,*) 'INIT_HALOMOD: Non-linear wavenumber [h/Mpc]:', REAL(hmod%knl)
    END IF

    hmod%neff=effective_index(hmod,cosm)

    IF(verbose) WRITE(*,*) 'INIT_HALOMOD: Collapse n_eff:', REAL(hmod%neff)

    CALL fill_halo_concentration(z,hmod,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Halo concentration tables filled'
       WRITE(*,*) 'INIT_HALOMOD: Minimum concentration:', REAL(hmod%c(hmod%n))
       WRITE(*,*) 'INIT_HALOMOD: Maximum concentration:', REAL(hmod%c(1))
    END IF

    IF(slow) THEN
       A0=one_halo_amplitude(hmod,cosm)
       IF(verbose) THEN
          WRITE(*,*) 'INIT_HALOMOD: One-halo amplitude [Mpc/h]^3:', REAL(A0)
          WRITE(*,*) 'INIT_HALOMOD: One-halo amplitude [log10(M/[Msun/h])]:', REAL(log10(A0*comoving_matter_density(cosm)))
       END IF
    END IF

    IF(verbose) THEN
       WRITE(*,*) 'INIT_HALOMOD: Done'
       WRITE(*,*)
    END IF

    !IF(verbose) CALL print_halomodel_parameters(a,hmod,cosm)

  END SUBROUTINE init_halomod

  SUBROUTINE convert_mass_definitions(z,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rhom, rhoc, Dv, a

    !Scale factor
    a=scale_factor_z(z)

    !Get the densities
    rhom=comoving_matter_density(cosm)
    rhoc=comoving_critical_density(a,cosm)
    Dv=Delta_v(a,hmod,cosm)

    !Calculate Delta = 200, 500 and Delta_c = 200, 500 quantities
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,1.,hmod%r500,hmod%c500,hmod%m500,500.,1.,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,1.,hmod%r200,hmod%c200,hmod%m200,200.,1.,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,rhom,hmod%r500c,hmod%c500c,hmod%m500c,500.,rhoc,hmod%n)
    CALL convert_mass_definition(hmod%rv,hmod%c,hmod%m,Dv,rhom,hmod%r200c,hmod%c200c,hmod%m200c,200.,rhoc,hmod%n)

    hmod%has_mass_conversions=.TRUE.
    
  END SUBROUTINE convert_mass_definitions

  REAL FUNCTION total_stellar_mass_fraction(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    total_stellar_mass_fraction=rhobar(hmod%nu(1),hmod%large_nu,rhobar_star_integrand,hmod,cosm)
    total_stellar_mass_fraction=total_stellar_mass_fraction/comoving_matter_density(cosm)

  END FUNCTION total_stellar_mass_fraction

  SUBROUTINE init_galaxies(a,hmod,cosm)

    !Calcuate the number densities of galaxies
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu_min

    LOGICAL, PARAMETER :: verbose=.FALSE.
    
    nu_min=nu_M(hmod%mgal,a,hmod,cosm)
    hmod%n_c=rhobar(nu_min,hmod%large_nu,rhobar_central_integrand,hmod,cosm)
    hmod%n_s=rhobar(nu_min,hmod%large_nu,rhobar_satellite_integrand,hmod,cosm)
    hmod%n_g=hmod%n_c+hmod%n_s
    IF(verbose) THEN
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of central galaxies [(Mpc/h)^-3]:', REAL(hmod%n_c)
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of satellite galaxies [(Mpc/h)^-3]:', REAL(hmod%n_s)
       WRITE(*,*) 'INIT_GALAXIES: Comoving density of all galaxies [(Mpc/h)^-3]:', REAL(hmod%n_g)
       WRITE(*,*)
    END IF

    hmod%has_galaxies=.TRUE.
  
  END SUBROUTINE init_galaxies

  SUBROUTINE init_HI(a,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: nu_min, nu_max

    LOGICAL, PARAMETER :: verbose=.FALSE.

    !Calculate the normalisation constant for HI
    nu_min=nu_M(hmod%HImin,a,hmod,cosm)
    nu_max=nu_M(hmod%HImax,a,hmod,cosm)
    hmod%rho_HI=comoving_matter_density(cosm)*integrate_hmod(nu_min,nu_max,g_nu,hmod,hmod%acc_HMx,3)
    IF(verbose) THEN
       WRITE(*,*) 'INIT_HI: HI normalisation factor [log10(rho/(Msun/h)/(Mpc/h)^3)]:', REAL(log10(hmod%rho_HI))
       WRITE(*,*)
    END IF

    hmod%has_HI=.TRUE.

  END SUBROUTINE init_HI
  
  FUNCTION nu_R(R,a,hmod,cosm)

    !Calculates nu(R) where R is the Lagrangian halo radius
    IMPLICIT NONE
    REAL :: nu_R
    REAL, INTENT(IN) :: R, a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    nu_R=delta_c(a,hmod,cosm)/sigma(R,a,cosm)

  END FUNCTION nu_R

  FUNCTION nu_M(M,a,hmod,cosm)

    !Calculates nu(M) where M is the halo mass
    IMPLICIT NONE
    REAL :: nu_M
    REAL, INTENT(IN) :: M, a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: R

    R=radius_m(M,cosm)
    nu_M=nu_R(R,a,hmod,cosm)
    
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

    !Integrand for the number density of satellite galaxies
    IMPLICIT NONE
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: M

    M=M_nu(nu,hmod)    
    !rhobar_star_integrand=M*halo_star_fraction(M,cosm)*g_nu(nu,hmod)
    rhobar_star_integrand=halo_star_fraction(M,hmod,cosm)*g_nu(nu,hmod)
    
  END FUNCTION rhobar_star_integrand

!!$  REAL FUNCTION rhobar_HI_integrand(nu,hmod,cosm)
!!$
!!$    !Integrand for the mean HI density
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: nu
!!$    TYPE(halomod), INTENT(INOUT) :: hmod
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    REAL :: M
!!$
!!$    M=M_nu(nu,hmod)    
!!$    rhobar_HI_integrand=HI_fraction(M,cosm)*g_nu(nu,hmod)
!!$    
!!$  END FUNCTION rhobar_HI_integrand

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

  END FUNCTIOn one_halo_amplitude

  SUBROUTINE convert_mass_definition(r1,c1,m1,D1,rho1,r2,c2,m2,D2,rho2,n)

    !Converts mass definition from Delta_1 rho_1 overdense to Delta_2 rho_2 overdense
    !r1(n) - Array of initial virial radii [Mpc/h]
    !c1(n) - Array of initial halo concentration
    !M1(n) - Array of initial halo mass [Msun/h]
    !D1 - Initial halo overdensity definition (e.g., 200, 500)
    !rho1 - Initial halo overdensity defintion (critical or mass)
    !r1(n) - Output array of virial radii with new definition [Mpc/h]
    !c1(n) - Output array of halo concentration with new definition
    !M1(n) - Output array of halo mass with new definition [Msun/h]
    !D1 - Final halo overdensity definition (e.g., 200, 500)
    !rho1 - Final halo overdensity defintion (critical or mass)
    !n - Number of entries in tables
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: r1(n), c1(n), m1(n)
    REAL, INTENT(OUT) :: r2(n), c2(n), m2(n)
    REAL, INTENT(IN) :: D1, D2, rho1, rho2
    REAL :: f(n)!, r(n)
    REAL :: rmin, rmax, rs
    REAL, ALLOCATABLE :: r(:)
    INTEGER :: i, j

    LOGICAL, PARAMETER :: verbose = .FALSE.

    IF(verbose) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Converting mass definitions:'
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Initial overdensity:', D1
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Final overdensity:', D2
    END IF

    !Make an array of general 'r' values for the solution later
    rmin=r1(1)/10. !Should be sufficient for reasonable cosmologies
    rmax=r1(n)*10. !Should be sufficient for reasonable cosmologies
    CALL fill_array(log(rmin),log(rmax),r,n) !Necessary to log space
    r=exp(r)
    
    IF(verbose) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: rmin:', rmin
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: rmax:', rmax
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: nr:', n
    END IF
    
    !Now use the find algorithm to invert L(r_i)=R(r_j) so that r_j=R^{-1}[L(r_i)]
    IF(verbose) THEN
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

       IF(verbose) THEN
          WRITE(*,fmt='(ES15.7,6F15.7)') m1(i), r1(i), c1(i), m2(i)/m1(i), r2(i)/r1(i), c2(i)/c1(i), NFW_factor(c2(i))/NFW_factor(c1(i))
       END IF

    END DO

    IF(verbose) WRITE(*,*) '========================================================================================================'

    IF(verbose) THEN
       WRITE(*,*) 'CONVERT_MASS_DEFINITION: Done'
       WRITE(*,*)
       STOP
    END IF

  END SUBROUTINE convert_mass_definition

  REAL FUNCTION NFW_factor(x)

    !The NFW 'mass' factor that crops up all the time
    !This is X(c) in M(r) = M X(r/rs) / X(c)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    NFW_factor=log(1.+x)-x/(1.+x)

  END FUNCTION NFW_factor

  FUNCTION radius_m(m,cosm)

    !The comoving radius corresponding to mass M in a homogeneous universe
    IMPLICIT NONE
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(INOUT) :: cosm

    radius_m=(3.*m/(4.*pi*comoving_matter_density(cosm)))**(1./3.)

  END FUNCTION radius_m

  FUNCTION virial_radius(m,a,hmod,cosm)

    !The comoving halo virial radius
    IMPLICIT NONE
    REAL :: virial_radius
    REAL, INTENT(IN) :: m, a
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    virial_radius=(3.*m/(4.*pi*comoving_matter_density(cosm)*Delta_v(a,hmod,cosm)))**(1./3.)

  END FUNCTION virial_radius

  REAL FUNCTION effective_index(hmod,cosm)

    !Power spectrum slope a the non-linear scale
    IMPLICIT NONE
    !REAL :: effective_index
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod

    !Numerical differentiation to find effective index at collapse
    effective_index=-3.-derivative_table(log(hmod%rnl),log(hmod%rr),log(hmod%sig**2),hmod%n,3,3)

    !For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
    !In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    IF(effective_index<cosm%n-4.) effective_index=cosm%n-4.
    IF(effective_index>cosm%n)      effective_index=cosm%n

  END FUNCTION effective_index

  SUBROUTINE fill_halo_concentration(z,hmod,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: mstar, m, zc
    INTEGER :: i

    !iconc = 1: Full Bullock et al. (2001)
    !iconc = 2: Simple Bullock et al. (2001)
    !iconc = 3: Duffy et al. (2008): mean
    !iconc = 4: Duffy et al. (2008): virial

    !Any initialisation for the c(M) relation goes here
    IF(hmod%iconc==1) THEN
       !Fill the collapse z look-up table
       CALL zcoll_Bullock(z,hmod,cosm)
    ELSE IF(hmod%iconc==2) THEN
       mstar=hmod%mnl
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
          hmod%c(i)=conc_Bullock_simple(m,mstar)
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

  FUNCTION gas_correction(m,hmod,cosm)

    IMPLICIT NONE
    REAL :: gas_correction
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    gas_correction=(1.+(hmod%eps-1.)*halo_boundgas_fraction(m,hmod,cosm)/(cosm%Om_b/cosm%Om_m))
    
  END FUNCTION gas_correction

  SUBROUTINE Dolag_correction(hmod,cosm)

    IMPLICIT NONE
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: g_LCDM, g_wCDM, ainf, f
    TYPE(cosmology) :: cosm_LCDM

    REAL, PARAMETER :: zinf=100. ! I changed this from 10 -> 100 to make it more infinite, I used zinf=10 for Mead (2017)!

    ! The 'infinite' scale factor
    ainf=scale_factor_z(zinf)

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

  FUNCTION conc_Bullock(z,zc)

    IMPLICIT NONE
    REAL :: conc_Bullock
    REAL, INTENT(IN) :: z, zc

    !PARAMETERS
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
       dc=delta_c(a,hmod,cosm)

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

  FUNCTION p_2h(ih,wk,k,z,plin,hmod,cosm)

    !Produces the 'two-halo' power
    IMPLICIT NONE
    REAL :: p_2h    
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL, INTENT(IN) :: k, z, plin, wk(hmod%n,2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(IN) :: ih(2)
    REAL :: sigv, frac, rhom, a
    REAL :: nu, m, m0, wki(2)
    REAL :: integrand1(hmod%n,2), integrand2(hmod%n,2)
    REAL :: sum1(2), sum2(2)
    INTEGER :: i, j

    !Get the scale factor
    a=scale_factor_z(z)

    rhom=comoving_matter_density(cosm)

    IF(hmod%ip2h==1) THEN

       ! Simply linear theory
       p_2h=plin

    ELSE IF(hmod%ip2h==3) THEN

       ! Damped BAO linear theory
       p_2h=p_dewiggle(k,z,hmod,cosm)
       
    ELSE

       DO i=1,hmod%n

          !Some variables to make equations cleaner below
          m=hmod%m(i)
          nu=hmod%nu(i)
          wki=wk(i,:)

          DO j=1,2

             !Linear bias term, standard two-halo term integral
             integrand1(i,j)=g_nu(nu,hmod)*b_nu(nu,hmod)*wki(j)/m

             IF(hmod%ibias==2) THEN
                !Second-order bias term
                integrand2(i,j)=g_nu(nu,hmod)*b2_nu(nu,hmod)*wki(j)/m
             END IF

          END DO

       END DO

       !Evaluate these integrals from the tabled values
       DO j=1,2
          sum1(j)=integrate_table(hmod%nu,integrand1(:,j),hmod%n,1,hmod%n,3)
       END DO

       IF(hmod%ip2h_corr==1) THEN
          !Do nothing in this case
          !There will be large errors if any signal is from low-mass haloes
          !e.g., for the matter power spectrum
       ELSE IF(hmod%ip2h_corr==2) THEN
          !Add on the value of integral b(nu)*g(nu) assuming W(k)=1
          !Advised by Yoo et al. (????) and Cacciato et al. (2012)
          !THIS WILL NOT WORK FOR FIEDS THAT DO NOT HAVE MASS FUNCTIONS DEFINED
          STOP 'P_2H: This will not work for fields that do not have mass fractions defined'
          DO j=1,2
             sum1(j)=sum1(j)+hmod%gbmin*halo_fraction(ih(j),m,hmod,cosm)/rhom
          END DO
       ELSE IF(hmod%ip2h_corr==3) THEN
          !Put the missing part of the integrand as a delta function at the low-mass limit of the integral
          !I think this is the best thing to do
          m0=hmod%m(1)
          wki=wk(1,:)
          DO j=1,2             
             sum1(j)=sum1(j)+hmod%gbmin*wki(j)/m0
          END DO
       ELSE
          STOP 'P_2h: Error, ip2h_corr not specified correctly'
       END IF

       p_2h=plin*sum1(1)*sum1(2)*(rhom**2)

       IF(hmod%ibias==2) THEN
          !Second-order bias correction
          !This needs to have the property that \int f(nu)b2(nu) du = 0
          !This means it is hard to check that the normalisation is correct
          !e.g., how much do low mass haloes matter
          !Varying mmin *does* make a difference to the values of the integrals
          !sum21=integrate_table(hmod%nu,integrand21,hmod%n,1,hmod%n,3)
          !sum22=integrate_table(hmod%nu,integrand22,hmod%n,1,hmod%n,3)
          DO j=1,2
             sum2(j)=integrate_table(hmod%nu,integrand2(:,j),hmod%n,1,hmod%n,3)
          END DO
          p_2h=p_2h+(plin**2)*sum2(1)*sum2(2)*rhom**2
       END IF

    END IF

    !Apply the damping to the two-halo term
    IF(hmod%i2hdamp .NE. 1) THEN
       !Two-halo damping parameters
       sigv=hmod%sigv
       frac=fdamp(a,hmod,cosm)
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

    !For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0.
    IF(p_2h<0.) STOP 'P_2h: Caution! P_2h < 0., this used to be fixed by setting p_2h=0. explicitly'

  END FUNCTION p_2h

  FUNCTION p_1h(wk2,k,z,hmod,cosm)

    !Calculates the one-halo term
    IMPLICIT NONE
    REAL :: p_1h    
    TYPE(halomod), INTENT(INOUT) :: hmod
    REAL, INTENT(IN) :: k, z, wk2(hmod%n)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m, g, fac, ks, a
    REAL, ALLOCATABLE :: integrand(:)
    INTEGER :: i

    !REAL, PARAMETER :: zmax_1halo=5.

    !IF(z>zmax_1halo) THEN
    !
    !   p_1h=0.
    !
    !ELSE

    !Get the scale factor
    a=scale_factor_z(z)

    ALLOCATE(integrand(hmod%n))
    integrand=0.

    !Calculates the value of the integrand at all nu values!
    DO i=1,hmod%n
       g=g_nu(hmod%nu(i),hmod)
       m=hmod%m(i)
       integrand(i)=g*wk2(i)/m
    END DO

    !Carries out the integration
    !Important to use basic trapezium rule because the integrand is messy due to rapid oscillations in W(k)
    p_1h=comoving_matter_density(cosm)*integrate_table(hmod%nu,integrand,hmod%n,1,hmod%n,1)*(4.*pi)*(k/(2.*pi))**3

    DEALLOCATE(integrand)

    !Damping of the 1-halo term at very large scales
    ks=kstar(hmod,cosm)       

    IF(ks>0.) THEN

       IF(hmod%i1hdamp==1) THEN
          !Do nothing in this case
       ELSE IF(hmod%i1hdamp==2) THEN
          IF((k/ks)**2>7.) THEN
             !Prevents problems if k/ks is very large
             fac=0.
          ELSE
             fac=exp(-((k/ks)**2))
          END IF
          p_1h=p_1h*(1.-fac)
       ELSE IF(hmod%i1hdamp==3) THEN
          !Note that the power here should be 4 because it multiplies Delta^2(k) ~ k^3 at low k (NOT 7)
          !Want f(k<<ks) ~ k^4; f(k>>ks) = 1
          fac=1./(1.+(ks/k)**4)
          p_1h=p_1h*fac
       ELSE
          STOP 'P_1H: Error, i1hdamp not specified correctly'          
       END IF

    END IF

    !END IF

  END FUNCTION p_1h

  FUNCTION p_1void(k,hmod)!,cosm)

    IMPLICIT NONE
    REAL :: p_1void
    REAL, INTENT(IN) :: k
    TYPE(halomod), INTENT(INOUT) :: hmod
    !TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: dc, wk, V, rvoid, rcomp, nu
    REAL :: integrand(hmod%n)
    INTEGER :: i, n

    !Parameters
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

    p_1void=p_1void*(4.*pi)*(k/(2.*pi))**3

  END FUNCTION p_1void

  FUNCTION win_type(real_space,itype,ipnh,k,z,m,rv,rs,hmod,cosm)

    !Selects the halo profile type
    IMPLICIT NONE
    REAL :: win_type
    REAL, INTENT(IN) :: k, z, m, rv, rs
    INTEGER, INTENT(IN) :: itype, ipnh
    LOGICAL, INTENT(IN) :: real_space
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    

    IF(itype==-1) THEN
       !Overdensity if all the matter were CDM
       win_type=win_DMONLY(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==0) THEN
       !Matter overdensity (sum of CDM, gas, stars)
       win_type=win_total(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==1) THEN
       !CDM overdensity
       win_type=win_CDM(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==2) THEN
       !All gas, both bound and free overdensity
       win_type=win_gas(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==3) THEN
       !Stellar overdensity
       win_type=win_star(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==4) THEN
       !Bound gas overdensity
       win_type=win_boundgas(real_space,1,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==5) THEN
       !Free gas overdensity
       win_type=win_freegas(real_space,1,ipnh,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==6) THEN
       !Electron pressure
       win_type=win_electron_pressure(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==7) THEN
       !Void
       win_type=win_void(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==8) THEN
       !Compensated void
       win_type=win_compensated_void(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==9) THEN
       !Central galaxies
       win_type=win_centrals(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==10) THEN
       !Satellite galaxies
       win_type=win_satellites(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==11) THEN
       !All galaxies
       win_type=win_galaxies(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE IF(itype==12) THEN
       !Neutral hydrogen - HI
       win_type=win_HI(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE
       STOP 'WIN_TYPE: Error, itype not specified correclty' 
    END IF

  END FUNCTION win_type

  FUNCTION win_total(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)

    !The halo profile of all the matter
    IMPLICIT NONE
    REAL :: win_total
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: ipnh
    REAL, INTENT(IN) :: k, z, rv, rs, m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    win_total=win_CDM(real_space,k,z,m,rv,rs,hmod,cosm)+win_gas(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)+win_star(real_space,k,z,m,rv,rs,hmod,cosm)

  END FUNCTION win_total

  FUNCTION win_DMONLY(real_space,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for all matter under the assumption that it is all CDM
    IMPLICIT NONE
    REAL :: win_DMONLY
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax, c, rss
    REAL :: crap

    !Set the DMONLY halo model
    !1 - Analyical NFW
    !2 - Non-analytical NFW (good for testing W(k) functions)
    !3 - Tophat
    !4 - Delta function
    INTEGER, PARAMETER :: imod=1

    !Prevent compile-time warnings
    crap=z
    crap=hmod%A

    IF(imod==1) THEN
       irho=5 !Analytical NFW
    ELSE IF(imod==2) THEN
       irho=4 !Non-analyical NFW
    ELSE IF(imod==3) THEN
       irho=2 !Tophat
    ELSE IF(imod==4) THEN
       irho=0 !Delta function
    ELSE
       STOP 'WIN_DMONLY: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=rv

    ! Force it to use the gravity-only concentration relation (unapply gas correction)
    c=rv/rs
    c=c/gas_correction(m,hmod,cosm)
    rss=rv/c
    !rss=rs

    IF(real_space) THEN
       r=k
       win_DMONLY=rho(r,rmin,rmax,rv,rss,zero,zero,irho)
       win_DMONLY=win_DMONLY/normalisation(rmin,rmax,rv,rss,zero,zero,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_DMONLY=m*win_norm(k,rmin,rmax,rv,rss,zero,zero,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_DMONLY

  FUNCTION win_CDM(real_space,k,z,m,rv,rs,hmod,cosm)

    !The halo profile for CDM
    IMPLICIT NONE
    REAL :: win_CDM
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax
    REAL :: crap

    !Set the model
    !1 - NFW
    INTEGER, PARAMETER :: imod=1

    IF(imod==1) THEN
       !Analytical NFW
       irho=5
    ELSE
       STOP 'WIN_CDM: Error, imod specified incorrectly'
    END IF

    !Prevent compile-time warnings
    crap=z
    crap=hmod%A

    rmin=0.
    rmax=rv

    IF(real_space) THEN
       r=k
       win_CDM=rho(r,rmin,rmax,rv,rs,zero,zero,irho)
       win_CDM=win_CDM/normalisation(rmin,rmax,rv,rs,zero,zero,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_CDM=m*win_norm(k,rmin,rmax,rv,rs,zero,zero,irho)/comoving_matter_density(cosm)
    END IF

    win_CDM=halo_CDM_fraction(m,hmod,cosm)*win_CDM

  END FUNCTION win_CDM

  FUNCTION win_gas(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for gas
    IMPLICIT NONE
    REAL :: win_gas
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: ipnh
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    win_gas=win_boundgas(real_space,1,k,z,m,rv,rs,hmod,cosm)+win_freegas(real_space,1,ipnh,k,z,m,rv,rs,hmod,cosm)

  END FUNCTION win_gas

  FUNCTION win_boundgas(real_space,itype,k,z,m,rv,rs,hmod,cosm)

    ! Halo profile for the electron pressure of the bound component
    IMPLICIT NONE
    REAL :: win_boundgas
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rho0, T0, r
    REAL :: rmin, rmax, p1, p2
    INTEGER :: irho_density, irho_electron_pressure
    REAL :: crap

    ! Select model
    ! 1 - Simplified Komatsu & Seljak (2001) gas model
    ! 2 - Isothermal beta model
    ! 3 - Full Komatsu & Seljak (2001) gas model
    INTEGER, PARAMETER :: imod=1

    ! Stop compile-time warnings
    crap=z
    crap=hmod%A

    ! Initially set the halo parameters to zero
    p1=0.
    p2=0.

    IF(imod==1 .OR. imod==3) THEN
       ! Set KS profile
       IF(imod==1) THEN
          irho_density=11
          irho_electron_pressure=13
       ELSE IF(imod==3) THEN
          irho_density=21
          irho_electron_pressure=23
       END IF
       rmin=0.
       rmax=rv
       p1=hmod%Gamma
    ELSE IF(imod==2) THEN
       irho_density=6 ! Set cored isothermal profile with beta=2/3 
       irho_electron_pressure=irho_density ! okay to use density for electron pressure because temperature is constant
       rmin=0.
       rmax=rv
    ELSE        
       STOP 'WIN_BOUNDGAS: Error, imod not specified correctly'
    END IF

    IF(itype==1) THEN

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

    ELSE IF(itype==2) THEN

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
       rho0=rho0*msun/mpc/mpc/mpc !Overflow with REAL*4 if you use mpc**3
       rho0=rho0*cosm%h**2 !Absorb factors of h, so now [kg/m^3]

       ! Calculate the value of the temperature prefactor [K]
       T0=HMx_alpha(hmod)*virial_temperature(m,rv,cosm)

       ! Convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
       win_boundgas=win_boundgas*(rho0/(mp*cosm%mue))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
       win_boundgas=win_boundgas/(eV*cm**(-3)) ! Change units to pressure in [eV/cm^3]
       win_boundgas=win_boundgas*cosm%mue/cosm%mup ! Convert from total thermal pressure to electron pressure

    ELSE

       STOP 'WIN_BOUNDGAS: Error, itype not specified correctly'

    END IF

  END FUNCTION win_boundgas

  FUNCTION win_freegas(real_space,itype,ipnh,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for the free gas component
    IMPLICIT NONE
    REAL :: win_freegas
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: itype, ipnh
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: re, rmin, rmax, r, A, rho0, rhov, T0, p1, p2, beta, c, thing, m0
    INTEGER :: irho_density, irho_electron_pressure
    LOGICAL :: match_electron_pressure

    !Set the model
    !1 - Isothermal model (out to 2rv)
    !2 - Ejected gas model from Schneider (2015)
    !3 - Isothermal shell that connects electron pressure and density to boundgas at rv
    !4 - Komatsu-Seljak continuation
    !5 - Power-law continuation
    !6 - Cubic profile
    !7 - Smoothly distributed (physically dubious)
    !8 - Delta function (physically dubious)
    INTEGER, PARAMETER :: imod=8

    IF(hmod%smooth_freegas .AND. ipnh==1) THEN

       !This is only used for the one-halo term if one is considering a smooth free gas fraction
       win_freegas=0.

    ELSE

       !Enable to force the electron pressure to be matched at the virial radius
       !This is enabled by default for some halo gas/pressure models
       match_electron_pressure=.FALSE.

       !Set the halo 'parameter' variables to zero initially
       p1=0.
       p2=0.

       IF(halo_freegas_fraction(m,hmod,cosm)==0.) THEN

          !Sometimes the freegas fraction will be zero, in which case this avoids problems
          win_freegas=0.

       ELSE

          IF(imod==1) THEN

             !Simple isothermal model, motivated by constant velocity and rate expulsion
             irho_density=1
             irho_electron_pressure=irho_density !Okay because T is constant
             rmin=0.
             rmax=2.*rv

          ELSE IF(imod==2) THEN

             !Ejected gas model from Schneider (2015)
             irho_density=10
             irho_electron_pressure=irho_density !Okay because T is constant
             rmin=rv
             re=rv
             p1=re
             rmax=15.*re !Needs to be such that integral converges (15rf seems okay)

          ELSE IF(imod==3) THEN

             !Now do isothermal shell connected to the KS profile continuously
             irho_density=16
             irho_electron_pressure=irho_density !Okay because T is constant

             !Isothermal model with continuous link to KS
             rhov=win_boundgas(.TRUE.,1,rv,z,m,rv,rs,hmod,cosm) !This is the value of rho at the halo boundary for the bound gas           
             A=rhov/rho(rv,0.,rv,rv,rs,p1,p2,irho_density) !This is A, as in A/r^2

             rmin=rv
             rmax=rv+halo_freegas_fraction(m,hmod,cosm)/(4.*pi*A) !This ensures density continuity and mass conservation

             c=10. !How many times larger than the virial radius can the gas cloud go?          
             IF(rmax>c*rv) rmax=c*rv !This needs to be set otherwise get huge decrement in gas power at large scales
             match_electron_pressure=.TRUE. !Match the electron pressure at the boundary

          ELSE IF(imod==4) THEN

             !Ejected gas is a continuation of the KS profile
             irho_density=11 !KS
             irho_electron_pressure=13 !KS
             rmin=rv
             rmax=2.*rv
             p1=hmod%Gamma

          ELSE IF(imod==5) THEN

             m0=1e14

             IF(m<m0) THEN

                irho_density=0
                irho_electron_pressure=irho_density
                rmin=0.
                rmax=rv

             ELSE

                !Set the density profile to be the power-law profile
                irho_density=17
                irho_electron_pressure=irho_density !Not okay

                !Calculate the KS index at the virial radius
                c=rv/rs
                beta=(c-(1.+c)*log(1.+c))/((1.+c)*log(1.+c))
                beta=beta/(hmod%Gamma-1.) !This is the power-law index at the virial radius for the KS gas profile
                p1=beta
                !WRITE(*,*) 'Beta:', beta, log10(m)
                IF(beta<=-3.) beta=-2.9 !If beta<-3 then there is only a finite amount of gas allowed in the free component

                !Calculate the density at the boundary of the KS profile
                rhov=win_boundgas(.TRUE.,1,rv,z,m,rv,rs,hmod,cosm)
                !WRITE(*,*) 'rho_v:', rhov

                !Calculate A as in rho(r)=A*r**beta
                A=rhov/rho(rv,0.,rv,rv,rs,p1,p2,irho_density)
                !WRITE(*,*) 'A:', A

                !Set the minimum radius for the power-law to be the virial radius
                rmin=rv
                !WRITE(*,*) 'rmin:', rmin

                !Set the maximum radius so that it joins to KS profile seamlessly
                thing=(beta+3.)*halo_freegas_fraction(m,hmod,cosm)/(4.*pi*A)+(rhov*rv**3)/A
                !WRITE(*,*) 'thing:', thing
                IF(thing>0.) THEN
                   !This then fixes the condition of contiunity in amplitude and gradient
                   rmax=thing**(1./(beta+3.))
                ELSE
                   !If there are no sohmodions then fix to 10rv and accept discontinuity
                   !There may be no sohmodion if there is a lot of free gas and if beta<-3
                   rmax=10.*rv
                END IF
                !WRITE(*,*) 'rmax 2:', rmax

             END IF

          ELSE IF(imod==6) THEN

             !Cubic profile
             rmin=rv
             rmax=3.*rv
             irho_density=18
             irho_electron_pressure=irho_density

          ELSE IF(imod==7) THEN

             !Smooth profile
             rmin=0.
             rmax=rv
             irho_density=19
             irho_electron_pressure=irho_density

          ELSE IF(imod==8) THEN

             !Delta function
             rmin=0.
             rmax=rv
             irho_density=0
             irho_electron_pressure=irho_density

          ELSE
             STOP 'WIN_FREEGAS: Error, imod_freegas specified incorrectly'
          END IF

          !Density profile
          IF(itype==1) THEN

             !Density profile of free gas
             IF(real_space) THEN
                r=k
                win_freegas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_density)
                win_freegas=win_freegas/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density)
             ELSE
                !Properly normalise and convert to overdensity
                win_freegas=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_density)/comoving_matter_density(cosm)
             END IF

             win_freegas=halo_freegas_fraction(m,hmod,cosm)*win_freegas

          !Electron pressure profile
          ELSE IF(itype==2) THEN

             !If we are applying a pressure-matching condition
             IF(match_electron_pressure) THEN

                STOP 'WIN_FREEGAS: Check the units and stuff here *very* carefully'

                r=k
                IF(r>rmin .AND. r<rmax) THEN
                   !Only works for isothermal profile
                   win_freegas=win_boundgas(.TRUE.,2,rv,z,m,rv,rs,hmod,cosm)*(r/rv)**(-2)
                ELSE
                   win_freegas=0.
                END IF

             ELSE

                !Electron pressure profile of free gas
                IF(real_space) THEN
                   r=k
                   win_freegas=rho(r,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
                ELSE  
                   win_freegas=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)*normalisation(rmin,rmax,rv,rs,p1,p2,irho_electron_pressure)
                END IF

                !Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
                rho0=m*halo_freegas_fraction(m,hmod,cosm)/normalisation(rmin,rmax,rv,rs,p1,p2,irho_density) !rho0 in [(Msun/h)/(Mpc/h)^3]
                rho0=rho0*msun/Mpc/Mpc/Mpc !Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
                rho0=rho0*cosm%h**2 !Absorb factors of h, so now [kg/m^3]

                !This is the total thermal pressure of the WHIM
                T0=hmod%Twhim !Units are [K]

                !Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density) 
                win_freegas=win_freegas*(rho0/(mp*cosm%mup))*(kb*T0) !Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
                win_freegas=win_freegas/(eV*cm**(-3)) !Change units to pressure in [eV/cm^3]
                win_freegas=win_freegas*cosm%mue/cosm%mup !Convert from total thermal pressure to electron pressure

             END IF

          ELSE

             STOP 'WIN_FREEGAS: Error, itype not specified correctly'

          END IF

       END IF

    END IF

  END FUNCTION win_freegas

  FUNCTION win_star(real_space,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for stars
    IMPLICIT NONE
    REAL :: win_star
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: rstar, r, rmin, rmax, p1, p2
    REAL :: crap

    !Set the model
    !1 - Fedeli (2014) stellar distribution
    !2 - Schneider (2015) stellar distribution
    !3 - Delta function
    INTEGER, PARAMETER :: imod=1 !Set the model

    !To prevent compile-time warnings
    crap=rs
    crap=z
    crap=hmod%A

    !Initially set p1, p2
    p1=0.
    p2=0.

    IF(imod==1) THEN
       !Fedeli (2014)
       irho=7
       !rstar=0.1*rv
       rstar=hmod%rstar*rv
       p1=rstar
       rmax=rv !Set so that not too much bigger than rstar, otherwise bumps integration goes tits
    ELSE IF(imod==2) THEN
       !Schneider (2015), following Mohammed (2014)
       irho=9
       rstar=0.01*rv
       p1=rstar
       rmax=10.*rstar !Set so that not too much bigger than rstar, otherwise bumps integration goes tits
    ELSE IF(imod==3) THEN
       !Delta function
       irho=0
       !rmax=rv !Set this although it does not matter
    ELSE
       STOP 'WIN_STAR: Error, imod_star specified incorrectly'
    END IF

    rmin=0.

    IF(real_space) THEN
       r=k
       win_star=rho(r,rmin,rmax,rv,rs,p1,p2,irho)
       win_star=win_star/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
    ELSE
       !Properly normalise and convert to overdensity
       win_star=m*win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)/comoving_matter_density(cosm)
    END IF

    win_star=halo_star_fraction(m,hmod,cosm)*win_star

  END FUNCTION win_star

  FUNCTION win_electron_pressure(real_space,ipnh,k,z,m,rv,rs,hmod,cosm)

    !Halo electron pressure profile function for the sum of bound + unbound electron gas
    IMPLICIT NONE
    REAL :: win_electron_pressure
    LOGICAL, INTENT(IN) :: real_space
    INTEGER, INTENT(IN) :: ipnh
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%use_UPP) THEN
       !This overrides everything and just uses the UPP
       win_electron_pressure=UPP(real_space,k,z,m,rv,rs,hmod,cosm)
    ELSE
       win_electron_pressure=win_boundgas(real_space,2,k,z,m,rv,rs,hmod,cosm)+win_freegas(real_space,2,ipnh,k,z,m,rv,rs,hmod,cosm)
       !win_pressure=win_pressure*(1.+z)**3
    END IF

  END FUNCTION win_electron_pressure

  FUNCTION win_void(real_space,k,z,m,rv,rs,hmod,cosm)

    !Void profile
    IMPLICIT NONE
    REAL :: win_void
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax
    REAL :: crap

    !Set the void model
    !1 - Top-hat void
    INTEGER, PARAMETER :: imod=1

    !Stop compile-time warnings
    crap=z
    crap=hmod%A

    IF(imod==1) THEN
       !Top-hat
       irho=2
    ELSE
       STOP 'WIN_VOID: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=10.*rv

    IF(real_space) THEN
       r=k
       win_void=rho(r,rmin,rmax,rv,rs,zero,zero,irho)
       win_void=win_void/normalisation(rmin,rmax,rv,rs,zero,zero,irho)
    ELSE       
       win_void=m*win_norm(k,rmin,rmax,rv,rs,zero,zero,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_void

  FUNCTION win_compensated_void(real_space,k,z,m,rv,rs,hmod,cosm)

    !Profile for compensated voids
    IMPLICIT NONE
    REAL:: win_compensated_void
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax
    REAL :: crap

    !Set the void model
    !1 - Top-hat void
    INTEGER, PARAMETER :: imod=1

    !Stop compile-time warnings
    crap=z
    crap=hmod%A

    IF(imod==1) THEN
       !Top-hat
       irho=2
    ELSE
       STOP 'WIN_COMPENSATED_VOID: Error, imod specified incorrectly'
    END IF

    rmin=0.
    rmax=10.*rv

    IF(real_space) THEN
       r=k
       win_compensated_void=rho(r,rmin,rmax,rv,rs,zero,zero,irho)
       win_compensated_void=win_compensated_void/normalisation(rmin,rmax,rv,rs,zero,zero,irho)
    ELSE       
       win_compensated_void=m*win_norm(k,rmin,rmax,rv,rs,zero,zero,irho)/comoving_matter_density(cosm)
    END IF

  END FUNCTION win_compensated_void

  FUNCTION win_centrals(real_space,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for central galaxies
    IMPLICIT NONE
    REAL :: win_centrals
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax
    REAL :: crap
    REAL :: a

    a=scale_factor_z(z)    
    IF(hmod%has_galaxies .EQV. .FALSE.) CALL init_galaxies(a,hmod,cosm)

    !Stop compile-time warnings
    crap=z
    crap=hmod%A
    crap=cosm%A

    !Delta functions
    irho=0

    rmin=0.
    rmax=rv

    IF(real_space) THEN
       r=k
       win_centrals=rho(r,rmin,rmax,rv,rs,zero,zero,irho)
       win_centrals=win_centrals/normalisation(rmin,rmax,rv,rs,zero,zero,irho)
    ELSE      
       win_centrals=win_norm(k,rmin,rmax,rv,rs,zero,zero,irho)/hmod%n_c
    END IF

    win_centrals=N_centrals(m,hmod)*win_centrals

  END FUNCTION win_centrals

  FUNCTION win_satellites(real_space,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for satellite galaxies
    IMPLICIT NONE
    REAL :: win_satellites
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax
    REAL :: crap
    REAL :: a

    a=scale_factor_z(z)
    IF(hmod%has_galaxies .EQV. .FALSE.) CALL init_galaxies(a,hmod,cosm)
    
    !Stop compile-time warnings
    crap=z
    crap=hmod%A
    crap=cosm%A

    !NFW profile
    irho=5

    rmin=0.
    rmax=rv

    IF(real_space) THEN
       r=k
       win_satellites=rho(r,rmin,rmax,rv,rs,zero,zero,irho)
       win_satellites=win_satellites/normalisation(rmin,rmax,rv,rs,zero,zero,irho)
    ELSE
       win_satellites=win_norm(k,rmin,rmax,rv,rs,zero,zero,irho)/hmod%n_s
    END IF

    win_satellites=N_satellites(m,hmod)*win_satellites
    
  END FUNCTION win_satellites

  FUNCTION win_galaxies(real_space,k,z,m,rv,rs,hmod,cosm)

    !Halo profile for all galaxies
    IMPLICIT NONE
    REAL :: win_galaxies
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    win_galaxies=win_centrals(real_space,k,z,m,rv,rs,hmod,cosm)+win_satellites(real_space,k,z,m,rv,rs,hmod,cosm)
    
  END FUNCTION win_galaxies

  FUNCTION N_centrals(m,hmod)

    !The number of central galaxies as a function of halo mass
    IMPLICIT NONE
    INTEGER :: N_centrals
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m<hmod%mgal) THEN
       N_centrals=0
    ELSE
       N_centrals=1
    END IF
    
  END FUNCTION N_centrals

  FUNCTION N_satellites(m,hmod)

    !The number of satellite galxies as a function of halo mass
    IMPLICIT NONE
    INTEGER :: N_satellites
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m<hmod%mgal) THEN
       N_satellites=0
    ELSE
       N_satellites=CEILING(m/hmod%mgal)-1
    END IF
    
  END FUNCTION N_satellites

  FUNCTION N_galaxies(m,hmod)

    !The number of central galaxies as a function of halo mass
    IMPLICIT NONE
    INTEGER :: N_galaxies
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    N_galaxies=N_centrals(m,hmod)+N_satellites(m,hmod)
    
  END FUNCTION N_galaxies

  FUNCTION win_HI(real_space,k,z,m,rv,rs,hmod,cosm)

    IMPLICIT NONE
    REAL :: win_HI
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: irho
    REAL :: r, rmin, rmax
    REAL :: crap
    REAL :: a

    a=scale_factor_z(z)
    IF(hmod%has_HI .EQV. .FALSE.) CALL init_HI(a,hmod,cosm)

    !Stop compile-time warnings
    crap=z
    crap=hmod%A
    crap=cosm%A

    !NFW profile
    irho=5

    rmin=0.
    rmax=rv

    IF(real_space) THEN
       r=k
       win_HI=rho(r,rmin,rmax,rv,rs,zero,zero,irho)
       win_HI=win_HI/normalisation(rmin,rmax,rv,rs,zero,zero,irho)
    ELSE
       !win_HI=win_norm(k,rmin,rmax,rv,rs,zero,zero,irho) !Wrong, but is what I first sent Richard and Kiyo
       win_HI=m*win_norm(k,rmin,rmax,rv,rs,zero,zero,irho)/hmod%rho_HI
    END IF

    win_HI=HI_fraction(m,hmod)*win_HI

  END FUNCTION win_HI

  REAL FUNCTION HI_fraction(m,hmod)

    IMPLICIT NONE
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(m>=hmod%HImin .AND. m<=hmod%HImax) THEN
       HI_fraction=1.
    ELSE
       HI_fraction=0.
    END IF
    
  END FUNCTION HI_fraction

  FUNCTION virial_temperature(M,rv,cosm)

    !Halo virial temperature in Kelvin
    !Calculates the temperature as if pristine gas falls into the halo
    !Energy is equally distributed between the particles
    IMPLICIT NONE
    REAL :: virial_temperature
    REAL :: M, rv !Virial mass and radius
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: modes=3. !0.5 k_BT per mode, 3 modes for 3 dimensions 

    virial_temperature=bigG*((M*msun)*mp*cosm%mup)/(rv*mpc)
    virial_temperature=virial_temperature/(kb*modes/2.) !Convert to temperature from energy

  END FUNCTION virial_temperature

  FUNCTION UPP(real_space,k,z,m,rv,rs,hmod,cosm)

    !Universal electron pressure profile (Arnaud et al. 2010; arxiv:0910.1234)
    !Note *very* well that this is for *electron* pressure (see Arnaud 2010 and my notes on this)
    !Note that is is also the physical pressure, and relates to the comoving pressure via (1+z)^3
    IMPLICIT NONE
    REAL :: UPP
    LOGICAL, INTENT(IN) :: real_space
    REAL, INTENT(IN) :: k, z, m, rv, rs
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(halomod), INTENT(INOUT) :: hmod  
    REAL :: r500c, rmin, rmax, a, r, m500c, E

    REAL, PARAMETER :: alphap=0.12 !Exponent correction
    REAL, PARAMETER :: b=0. !Hydrostatic mass bias
    INTEGER, PARAMETER :: irho=14 !Set UPP profile

    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(z,hmod,cosm)

    !Get r500 for UPP
    r500c=exp(find(log(m),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ![Mpc/h]

    !Set the radius range for the profile
    rmin=0.
    rmax=rv

    a=scale_factor_z(z)
    IF(real_space) THEN
       r=k
       UPP=rho(r,rmin,rmax,rv,rs,r500c,zero,irho)
    ELSE
       UPP=winint(k,rmin,rmax,rv,rs,r500c,zero,irho,imeth_win)
    END IF

    !Upp, P(x), equation 4.1 in Ma et al. (2015)
    m500c=exp(find(log(m),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ![Msun/h]
    m500c=m500c*(1.-b) ![Msun/h]

    !Dimensionless Hubble parameter
    E=sqrt(Hubble2(a,cosm))

    !Pre-factors from equation 4.1 in Ma et al. (2015) [eV cm^-3 - no h factors]
    UPP=UPP*((m500c/2.1e14)**(alphap+2./3.))*(E**(8./3.))*1.65*(cosm%h/0.7)**2

    !The standard UPP is written is in physical units
    !scale to comoving using (1+z)^3 because pressure ~energy density
    UPP=UPP/(1.+z)**3

  END FUNCTION UPP

!!$  REAL FUNCTION m500c_m(m,z,hmod,cosm)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: m, z
!!$    TYPE(halomod), INTENT(INOUT) :: hmod
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(z,hmod,cosm)
!!$
!!$    m500c_m=exp(find(log(m),hmod%log_m,log(hmod%m500c),hmod%n,3,3,2)) ![Msun/h]
!!$
!!$  END FUNCTION m500c_m

!!$  REAL FUNCTION r500c_m(m,z,hmod,cosm)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: m, z
!!$    TYPE(halomod), INTENT(INOUT) :: hmod
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    IF(hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(z,hmod,cosm)
!!$
!!$    r500c_m=exp(find(log(m),hmod%log_m,log(hmod%r500c),hmod%n,3,3,2)) ![Msun/h]
!!$
!!$  END FUNCTION r500c_m

  FUNCTION win_norm(k,rmin,rmax,rv,rs,p1,p2,irho)

    !Calculates the normalised spherical Fourier Transform of the density profile
    !Note that this means win_norm(k->0)=1
    !and that win must be between 0 and 1
    IMPLICIT NONE
    REAL :: win_norm
    REAL, INTENT(IN) :: rmin, rmax, k, rv, rs, p1, p2
    REAL :: re, f1, f2, rstar, kstar
    INTEGER, INTENT(IN) :: irho

    IF(k==0.) THEN

       !If called for the zero mode (e.g. for the normalisation)
       win_norm=1.

    ELSE

       IF(irho==0) THEN
          !Delta function profile is not localised in Fourier Space
          win_norm=1.
       ELSE IF(irho==1) THEN
          win_norm=wk_isothermal(k*rmax)
       ELSE IF(irho==2) THEN
          !Analytic for top hat
          win_norm=wk_tophat(k*rmax)
       ELSE IF(irho==5) THEN
          !Analytic for NFW
          win_norm=win_NFW(k,rmax,rs)
       ELSE IF(irho==7) THEN
          !Analytic for Fedeli (2014) stellar profile
          rstar=p1
          kstar=k*rstar
          f1=kstar-exp(-rmax/rstar)*(sin(k*rmax)+kstar*cos(k*rmax))
          f2=kstar*(1.+kstar**2)
          win_norm=f1/f2
          !win_norm=1./(1.+kstar**2) !bigRstar -> infinity limit (rmax >> rstar)
       ELSE IF(irho==9) THEN
          !STOP 'WIN_NORM: Double check this result for irho=9'
          !Only valid if rmin=0 and rmax=inf
          rstar=p1
          win_norm=(sqrt(pi)/2.)*erf(k*rstar)/(k*rstar)
       ELSE IF(irho==10) THEN
          !For ejected gas profile
          re=p1
          win_norm=exp(-1.5*(k*re)**2.)
       ELSE IF(irho==16) THEN
          !Isothermal shells
          win_norm=wk_isothermal_2(k*rmax,k*rmin)
       ELSE IF(irho==19) THEN
          !Smooth profile (not sure this is physical)
          win_norm=0.
       ELSE IF(irho==20) THEN
          !Exponential profile
          re=p1
          win_norm=1./(1.+(k*re)**2)**2
       ELSE
          !Numerical integral over the density profile (slower)
          win_norm=winint(k,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)/normalisation(rmin,rmax,rv,rs,p1,p2,irho)
       END IF

    END IF

  END FUNCTION win_norm

  FUNCTION rhor2at0(irho)

    !This is the value of r^2 * rho(r) at r=0.
    !For most profiles this is zero, BUT not if rho(r->0) -> r^-2
    !Note if rho(r->0) -> r^n with n<-2 then the profile mass would diverge!

    IMPLICIT NONE
    REAL :: rhor2at0
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

  FUNCTION rho(r,rmin,rmax,rv,rs,p1,p2,irho)

    !This is an UNNORMALISED halo profile of any sort

    !Types of profile
    !================
    ! 0 - Delta function at r=0
    ! 1 - Isothermal: r^-2
    ! 2 - Top hat: r^0
    ! 3 - Moore (1999)
    ! 4 - NFW (1997)
    ! 5 - Analytic NFW
    ! 6 - Beta model with beta=2/3
    ! 7 - Star profile
    ! 8 - Komatsu & Seljak (2001) according to Schneider (2015)
    ! 9 - Stellar profile from Schneider (2015)
    !10 - Ejected gas profile (Schneider 2015)
    !11 - Simplified Komatsu & Seljak (2001) density
    !12 - Simplified Komatsu & Seljak (2001) temperature
    !13 - Simplified Komatsu & Seljak (2001) pressure
    !14 - Universal pressure profile
    !15 - Isothermal beta model, beta=0.86 (Ma et al. 2015)
    !16 - Isothermal shell
    !17 - Power-law profile
    !18 - Cubic profile: r^-3
    !19 - Smooth profile (rho = 0, not really physical)
    !20 - Exponential profile
    !21 - Full Komatsu & Seljak (2001) density
    !22 - Full Komatsu & Seljak (2001) temperature
    !23 - Full Komatsu & Seljak (2001) pressure

    IMPLICIT NONE
    REAL :: rho
    REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2 !Standard profile parameters
    INTEGER, INTENT(IN) :: irho
    REAL :: y, ct, t, c, beta, Gamma, r500c, rt, A, re, rstar, B!, eta0 !Derived parameters
    REAL :: f1, f2
    REAL :: crap

    !UPP parameters
    REAL, PARAMETER :: P0=6.41
    REAL, PARAMETER :: c500=1.81
    REAL, PARAMETER :: alpha_UPP=1.33
    REAL, PARAMETER :: beta_UPP=4.13
    REAL, PARAMETER :: gamma_UPP=0.31

    !To stop compile-time warnings
    crap=p2

    IF(r<rmin .OR. r>rmax) THEN
       !The profile is considered to be zero outside this region
       rho=0.
    ELSE
       IF(irho==0) THEN
          !Delta function
          !Do not assign any value to rho as this gets handled properly elsewhere
          !STOP 'RHO: You should not be here for a delta-function profile'
          rho=0. 
       ELSE IF(irho==1 .OR. irho==16) THEN
          !Isothermal
          rho=1./r**2
       ELSE IF(irho==2) THEN
          !Top hat
          rho=1.
       ELSE IF(irho==3) THEN
          !Moore (1999)
          y=r/rs
          rho=1./((y**1.5)*(1.+y**1.5))
       ELSE IF(irho==4 .OR. irho==5) THEN
          !NFW (1997)
          y=r/rs
          rho=1./(y*(1.+y)**2)
       ELSE IF(irho==6) THEN
          !Isothermal beta model (X-ray gas; SZ profiles; beta=2/3 fixed)
          !AKA 'cored isothermal profile'
          y=r/rs
          beta=2./3.
          rho=1./((1.+y**2)**(3.*beta/2.))
       ELSE IF(irho==7) THEN
          !Stellar profile from Fedeli (2014a)
          rstar=p1
          y=r/rstar
          rho=(1./y)*exp(-y)
       ELSE IF(irho==8) THEN
          !Komatsu & Seljak (2001) profile with NFW transition radius
          !VERY slow to calculate the W(k) for some reason
          !Also creates a weird upturn in P(k) that I do not think can be correct
          STOP 'RHO: This is fucked'
          t=sqrt(5.)
          rt=rv/t
          y=r/rs
          c=rs/rv
          ct=c/t
          Gamma=(1.+3.*ct)*log(1.+ct)/((1.+ct)*log(1.+ct)-ct)
          IF(r<=rt) THEN
             !Komatsu Seljak in the interior
             rho=(log(1.+y)/y)**Gamma
          ELSE
             !NFW in the outskirts
             A=((rt/rs)*(1.+rt/rs)**2)*(log(1.+rt/rs)/(rt/rs))**Gamma
             rho=A/(y*(1.+y)**2)
          END IF
       ELSE IF(irho==9) THEN
          !Stellar profile from Schneider (2015) via Mohammed (2014)
          rstar=p1
          rho=exp(-(r/(2.*rstar))**2)/r**2
          !Converting to y caused the integration to crash for some reason !?!
          !y=r/rs
          !rho=exp(-(y/2.)**2.)/y**2.
       ELSE IF(irho==10) THEN
          !Ejected gas profile from Schneider (2015)
          re=p1
          rho=exp(-0.5*(r/re)**2)
       ELSE IF(irho==11 .OR. irho==12 .OR. irho==13) THEN
          !Komatsu & Seljak (2001) profile
          !Gamma=1.18 !Recommended by Rabold (2017)
          Gamma=p1
          y=r/rs
          rho=log(1.+y)/y
          IF(irho==11) THEN
             !KS density profile
             rho=rho**(1./(Gamma-1.))
          ELSE IF(irho==12) THEN
             !KS temperature profile
             rho=rho
          ELSE IF(irho==13) THEN
             !KS pressure profile
             rho=rho**(Gamma/(Gamma-1.))
          END IF
       ELSE IF(irho==14) THEN
          !UPP is in terms of r500c, not rv
          r500c=p1
          !UPP funny-P(x), equation 4.2 in Ma et al. (2015)
          f1=(c500*r/r500c)**gamma_UPP
          f2=(1.+(c500*r/r500c)**alpha_UPP)**((beta_UPP-gamma_UPP)/alpha_UPP)
          rho=P0/(f1*f2)
       ELSE IF(irho==15) THEN
          !Isothermal beta model
          !beta=0.86 !from Ma et al. (2015)
          beta=p1
          rho=(1.+(r/rs)**2)**(-3.*beta/2.)
       ELSE IF(irho==16) THEN
          !Isothermal (shell)
          rho=1./r**2
       ELSE IF(irho==17) THEN
          !Power-law profile
          beta=p1
          rho=r**beta
       ELSE IF(irho==18) THEN
          !Cubic profile
          rho=r**(-3)
       ELSE IF(irho==19) THEN
          !Smooth profile
          rho=0.
       ELSE IF(irho==20) THEN
          !Exponential profile (HI from Padmanabhan et al. 2017)
          re=p1
          rho=exp(-r/re)
       ELSE IF(irho==21 .OR. irho==22 .OR. irho==23) THEN
          !Komatsu & Seljak (2001) profile
          Gamma=p1
          !c=rv/rs
          !eta0=2.235+0.202*(c-5.)-1.16e-3*(c-5.)**2
          !f1=(3./eta0)*((Gamma-1.)/Gamma)
          !f2=log(1.+c)/c-1./(1.+c)
          !B=f1/f2
          B=1.
          y=r/rs
          rho=1.-B*(1.-log(1.+y)/y)
          IF(irho==21) THEN
             !KS density profile
             rho=rho**(1./(Gamma-1.))
          ELSE IF(irho==22) THEN
             !KS temperature profile
             rho=rho
          ELSE IF(irho==23) THEN
             !KS pressure profile
             rho=rho**(Gamma/(Gamma-1.))
          END IF
       ELSE
          STOP 'RHO: Error, irho not specified correctly'
       END IF

    END IF

  END FUNCTION rho

  FUNCTION winint(k,rmin,rmax,rv,rs,p1,p2,irho,imeth)

    !Calculates W(k,M)
    IMPLICIT NONE
    REAL :: winint
    REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho, imeth

    !Integration order
    INTEGER, PARAMETER :: iorder=3

    !Integration method
    !imeth = 1 - normal integration
    !imeth = 2 - bumps with normal integration
    !imeth = 3 - storage integration
    !imeth = 4 - bumps with storage integration
    !imeth = 5 - linear bumps
    !imeth = 6 - cubic bumps
    !imeth = 7 - Hybrid with storage and cubic bumps

    !Bump methods go crazy with some star profiles (those that drop too fast)
    !You need to make sure that the rmax for the integration does not extend too far out

    !The hybrid method seems not to be faster for practical calculations here

    IF(imeth==1) THEN
       winint=winint_normal(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc_win)
    ELSE IF(imeth==2 .OR. imeth==4 .OR. imeth==5 .OR. imeth==6 .OR. imeth==7) THEN
       winint=winint_bumps(k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc_win,imeth)
    ELSE IF(imeth==3) THEN
       winint=winint_store(rmin,rmax,k,rmin,rmax,rv,rs,p1,p2,irho,iorder,acc_win)
    ELSE
       STOP 'WININT: Error, imeth not specified correctly'
    END IF

  END FUNCTION winint

  SUBROUTINE winint_diagnostics(rmin,rmax,rv,rs,p1,p2,irho,outfile)

    !Write out the winint integrand as a function of k
    IMPLICIT NONE
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j
    REAL :: r, k
    REAL, ALLOCATABLE :: integrand(:)

    REAL, PARAMETER :: kmin=1d-1
    REAL, PARAMETER :: kmax=1d2
    INTEGER, PARAMETER :: nr=256 !Number of points in r
    INTEGER, PARAMETER :: nk=16 !Number of points in k

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

    !Integration routine using 'normal' method to calculate the normalised halo FT
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

    INTEGER, PARAMETER :: nlim=3 !Do the bumps approximation after this number of bumps

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

  FUNCTION win_NFW(k,rv,rs)

    !The analytic normalised (W(k=0)=1) Fourier Transform of the NFW profile
    IMPLICIT NONE
    REAL :: win_NFW
    REAL, INTENT(IN) :: k, rv, rs
    REAL :: c, ks
    REAL :: si1, si2, ci1, ci2
    REAL :: p1, p2, p3
    REAL :: rmin, rmax

    c=rv/rs
    ks=k*rv/c

    si1=si(ks)
    si2=si((1.+c)*ks)
    ci1=ci(ks)
    ci2=ci((1.+c)*ks)

    p1=cos(ks)*(ci2-ci1)
    p2=sin(ks)*(si2-si1)
    p3=sin(ks*c)/(ks*(1.+c))

    win_NFW=p1+p2-p3
    rmin=0.
    rmax=rv
    win_NFW=4.*pi*win_NFW*(rs**3.)/normalisation(rmin,rmax,rv,rs,zero,zero,4)

  END FUNCTION win_NFW

  FUNCTION normalisation(rmin,rmax,rv,rs,p1,p2,irho)

    ! This calculates the normalisation of a halo
    ! This is the integral of 4pir^2*rho(r)*dr between rmin and rmax

    ! Profile results
    !  0 - Delta function (M = 1)
    !  1 - Isothermal (M = 4pi*rv)
    !  2 - Top hat (M = (4pi/3)*rv^3)
    !  3 - Moore (M = (8pi/3)*rv^3*ln(1+c^1.5)/c^3)
    !  4,5 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
    !  6 - Beta model with beta=2/3 (M = 4*pi*rs^3*(rv/rs-atan(rv/rs)))
    !  7 - Fedeli stellar model (M = 4*pi*rstar^2 * [1-exp(-rmax/rstar)*(1.+rmax/rstar)]
    !  9 - Stellar profile (Schneider (2015)
    ! 10 - Ejected gas profile (Schneider 2015)
    ! 16 - Isothermal shell (M = 4pi*(rmax-rmin))
    ! 18 - Cubic profile
    ! 19 - Smooth profile (physically dubious)

    IMPLICIT NONE
    REAL :: normalisation
    REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
    INTEGER, INTENT(IN) :: irho
    REAL :: cmax, re, rstar, beta

    IF(irho==0) THEN
       !Delta function
       normalisation=1.
    ELSE IF(irho==1 .OR. irho==16) THEN
       !Isothermal
       normalisation=4.*pi*(rmax-rmin)
    ELSE IF(irho==2) THEN
       !Top hat
       normalisation=4.*pi*(rmax**3-rmin**3)/3.
    ELSE IF(irho==3) THEN
       !Moore et al. (1999)
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=(2./3.)*4.*pi*(rs**3)*log(1.+cmax**1.5)
       END IF
    ELSE IF(irho==4 .OR. irho==5) THEN
       !NFW (1997)
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=4.*pi*(rs**3)*NFW_factor(cmax)!(log(1.+cmax)-cmax/(1.+cmax))
       END IF
    ELSE IF(irho==6) THEN
       !Beta model with beta=2/3
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          cmax=rmax/rs
          normalisation=4.*pi*(rs**3)*(cmax-atan(cmax))
       END IF
    ELSE IF(irho==7) THEN
       !Fedeli (2014) stellar model
       IF(rmin .NE. 0) THEN
          !I could actually derive an analytical expression here if this was ever necessary
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          !This would be even easier if rmax -> infinity (just 4*pi*rstar^2)
          rstar=p1
          normalisation=4.*pi*(rstar**3)*(1.-exp(-rmax/rstar)*(1.+rmax/rstar))
          !normalisation=4.*pi*rstar**3 !rmax/rstar -> infinity limit (rmax >> rstar)
       END IF
    ELSE IF(irho==9) THEN
       !Stellar profile from Schneider & Teyssier (2015)       
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          !Assumed to go on to r -> infinity
          rstar=p1
          normalisation=4.*(pi**(3./2.))*rstar
       END IF
    ELSE IF(irho==10) THEN
       !Ejected gas profile from Schneider (2015)
       !Assumed to go on to r -> infinity
       !IF(rmin .NE. 0.) STOP 'NORMALISATION: Error, normalisation of Schneider (2015) gas profile assumed to be from 0 -> inf
       IF(rmin .NE. 0.) THEN
          normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
       ELSE
          !Assumed to go on to r -> infinity
          re=p1
          normalisation=4.*pi*sqrt(pi/2.)*re**3
       END IF
    ELSE IF(irho==17) THEN
       !Power-law profile
       beta=p1
       normalisation=(4.*pi/(beta+3.))*(rmax**(beta+3.)-rmin**(beta+3.))
    ELSE IF(irho==18) THEN
       !Cubic profile
       normalisation=4.*pi*log(rmax/rmin)
    ELSE IF(irho==19) THEN
       normalisation=1.
    ELSE
       !Otherwise need to do the integral numerically
       !k=0 gives normalisation
       normalisation=winint(0.,rmin,rmax,rv,rs,p1,p2,irho,imeth_win)
    END IF

  END FUNCTION normalisation

  FUNCTION b_nu(nu,hmod)

    !Bias function selection!
    IMPLICIT NONE
    REAL :: b_nu
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       b_nu=b_ps(nu)
    ELSE IF(hmod%imf==2) THEN
       b_nu=b_st(nu)
    ELSE IF(hmod%imf==3) THEN
       b_nu=b_Tinker(nu)
    ELSE
       STOP 'B_NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b_nu

  FUNCTION b_ps(nu)

    !Press & Scheter (1974) halo bias
    IMPLICIT NONE
    REAL :: b_ps
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: dc=1.686 ! TODO: This is probably not consistent

    b_ps=1.+(nu**2-1.)/dc

  END FUNCTION b_ps

  FUNCTION b_st(nu)

    !Sheth & Tormen (1999) halo bias (equation 12 in 9901122)
    !Comes from peak-background split
    !Haloes defined with SO relative to mean matter density with SC Delta_v relation
    !A redshift dependent delta_c is used for barrier height, again from SC
    IMPLICIT NONE
    REAL :: b_st
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: dc=1.686 ! TODO: This is probably not consistent

    b_st=1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/dc

  END FUNCTION b_st

  FUNCTION b_Tinker(nu)

    !Tinker et al. (2010; 1001.3162) halo bias
    IMPLICIT NONE
    REAL :: b_Tinker
    REAL, INTENT(IN) :: nu

    !Delta_v=200,m and delta_c=1.686 are hard-coded
    REAL, PARAMETER :: Delta_v=337.2
    REAL, PARAMETER :: delta_c=1.686
    REAL, PARAMETER :: y=log10(Delta_v) !This is the Delta_v dependence
    REAL, PARAMETER :: bigA=1.0+0.24*y*exp(-(4./y)**4)
    REAL, PARAMETER :: a=0.44*y-0.88
    REAL, PARAMETER :: bigB=0.183
    REAL, PARAMETER :: b=1.5
    REAL, PARAMETER :: bigC=0.019+0.107*y+0.19*exp(-(4./y)**4)
    REAL, PARAMETER :: c=2.4

    STOP 'B_TINKER: Be careful, redshift dependence is missing'
    
    b_Tinker=1.-bigA*(nu**a)/(nu**a+delta_c**a)+bigB*nu**b+bigC*nu**c
    
  END FUNCTION b_Tinker

  FUNCTION b2_nu(nu,hmod)

    !Bias function selection!
    IMPLICIT NONE
    REAL :: b2_nu
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       b2_nu=b2_ps(nu)
    ELSE IF(hmod%imf==2) THEN
       b2_nu=b2_st(nu)
    ELSE IF(hmod%imf==3) THEN
       STOP 'B2_NU: Error, second-order bias not specified for Tinker mass function'
    ELSE
       STOP 'B2_NU: Error, imf not specified correctly'
    END IF

  END FUNCTION b2_nu

  FUNCTION b2_ps(nu)

    !Press & Schechter (1974) second order bias
    IMPLICIT NONE
    REAL :: b2_ps
    REAL, INTENT(IN) :: nu
    REAL :: eps1, eps2, E1, E2

    REAL, PARAMETER :: a2=-17./21.
    REAL, PARAMETER :: p=0.0
    REAL, PARAMETER :: q=1.0
    REAL, PARAMETER :: dc=1.686 ! TODO: This is probably not consistent

    STOP 'B2_PS: Check this very carefully'
    !I just took the ST form and set p=0 and q=1

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2_ps=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2_ps

  FUNCTION b2_st(nu)

    !Sheth, Mo & Tormen (2001) second-order bias
    IMPLICIT NONE
    REAL :: b2_st
    REAL, INTENT(IN) :: nu
    REAL :: eps1, eps2, E1, E2

    !Notation follows from Cooray & Sheth (2002) pp 25-26

    REAL, PARAMETER :: a2=-17./21.
    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: dc=1.686 ! TODO: This is probably not consistent

    eps1=(q*nu**2-1.)/dc
    eps2=(q*nu**2)*(q*nu**2-3.)/dc**2
    E1=(2.*p)/(dc*(1.+(q*nu**2)**p))
    E2=((1.+2.*p)/dc+2.*eps1)*E1

    b2_st=2.*(1.+a2)*(eps1+E1)+eps2+E2

  END FUNCTION b2_st

  FUNCTION g_nu(nu,hmod)

    !Mass function
    IMPLICIT NONE
    REAL :: g_nu
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    IF(hmod%imf==1) THEN
       g_nu=g_ps(nu)
    ELSE IF(hmod%imf==2) THEN
       g_nu=g_st(nu)
    ELSE IF(hmod%imf==3) THEN
       g_nu=g_Tinker(nu)
    ELSE
       STOP 'G_NU: Error, imf specified incorrectly'
    END IF

  END FUNCTION g_nu

  FUNCTION g_ps(nu)

    !Press & Scheter (1974) mass function!
    IMPLICIT NONE
    REAL :: g_ps
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: A=sqrt(2./pi)

    g_ps=A*exp(-(nu**2)/2.)

  END FUNCTION g_ps

  FUNCTION g_st(nu)

    !Sheth & Tormen (1999) mass function, equation (10) in arXiv:9901122
    !Note I use nu=dc/sigma(M) and this Sheth & Tormen (1999) use nu=(dc/sigma)^2, which accounts for some small differences
    !Haloes defined with SO relative to mean matter density with SC Delta_v relation
    !A redshift dependent delta_c is used for barrier height, again from SC
    IMPLICIT NONE
    REAL :: g_st
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: A=0.21616

    g_st=A*(1.+((q*nu**2)**(-p)))*exp(-q*nu**2/2.)

  END FUNCTION g_st

  FUNCTION g_Tinker(nu)

    !Tinker et al. (2010; 1001.3162) mass function (also 2008; xxxx.xxxx)
    IMPLICIT NONE
    REAL :: g_Tinker
    REAL, INTENT(IN) :: nu
    REAL :: alpha, beta, gamma, phi, eta

    !Hard-coded at z=0.
    REAL, PARAMETER :: z=0.
    REAL, PARAMETER :: Dv=200.

    !Parameter arrays from Tinker (2010)
    INTEGER, PARAMETER :: n=9 !Number of entries in parameter lists
    REAL, PARAMETER :: Delta_v(n)=[200.,300.,400.,600.,800.,1200.,1600.,2400.,3200.]
    REAL, PARAMETER :: alpha0(n)=[0.368,0.363,0.385,0.389,0.393,0.365,0.379,0.355,0.327]
    REAL, PARAMETER :: beta0(n)=[0.589,0.585,0.544,0.543,0.564,0.623,0.637,0.673,0.702]
    REAL, PARAMETER :: gamma0(n)=[0.864,0.922,0.987,1.09,1.20,1.34,1.50,1.68,1.81]
    REAL, PARAMETER :: phi0(n)=[-0.729,-0.789,-0.910,-1.05,-1.20,-1.26,-1.45,-1.50,-1.49]
    REAL, PARAMETER :: eta0(n)=[-0.243,-0.261,-0.261,-0.273,-0.278,-0.301,-0.301,-0.319,-0.336]

    !Delta_v dependence
    alpha=find(Dv,Delta_v,alpha0,n,3,3,2)
    beta=find(Dv,Delta_v,beta0,n,3,3,2)
    gamma=find(Dv,Delta_v,gamma0,n,3,3,2)
    phi=find(Dv,Delta_v,phi0,n,3,3,2)
    eta=find(Dv,Delta_v,eta0,n,3,3,2)

    !Redshift dependence
    beta=beta*(1.+z)**0.20
    gamma=gamma**(1.+z)**(-0.01)
    phi=phi*(1.+z)**(-0.08)
    eta=eta*(1.+z)**0.27

    !The actual mass function
    g_Tinker=alpha*(1.+(beta*nu)**(-2.*phi))*nu**(2.*eta)*exp(-0.5*gamma*nu**2)
    
  END FUNCTION g_Tinker

  FUNCTION gb_nu(nu,hmod)

    !g(nu) times b(nu)
    IMPLICIT NONE
    REAL :: gb_nu
    REAL, INTENT(IN) :: nu
    TYPE(halomod), INTENT(INOUT) :: hmod

    gb_nu=g_nu(nu,hmod)*b_nu(nu,hmod)

  END FUNCTION gb_nu

  FUNCTION wk_isothermal(x)

    !The normlaised Fourier Transform of an isothermal profile
    IMPLICIT NONE
    REAL :: wk_isothermal
    REAL, INTENT(IN) :: x

    REAL, PARAMETER :: dx=1e-3

    !Taylor expansion used for low |x| to avoid cancellation problems

    IF(ABS(x)<ABS(dx)) THEN
       !Taylor series at low x
       wk_isothermal=1.-(x**2)/18.
    ELSE
       wk_isothermal=Si(x)/x
    END IF

  END FUNCTION wk_isothermal

  FUNCTION wk_isothermal_2(x,y)

    !The normlaised Fourier Transform of an isothemral profile from x -> y
    IMPLICIT NONE
    REAL :: wk_isothermal_2
    REAL, INTENT(IN) :: x, y

    wk_isothermal_2=(Si(x)-Si(y))/(x-y)

  END FUNCTION wk_isothermal_2

  FUNCTION halo_fraction(itype,m,hmod,cosm)

    !Mass fraction of a type within a halo
    IMPLICIT NONE
    REAL :: halo_fraction
    INTEGER, INTENT(IN) :: itype
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    If(itype==-1 .OR. itype==0) THEN
       halo_fraction=1.
    ELSE IF(itype==1) THEN
       halo_fraction=halo_CDM_fraction(m,hmod,cosm)
    ELSE IF(itype==2) THEN
       halo_fraction=halo_gas_fraction(m,hmod,cosm)
    ELSE IF(itype==3) THEN
       halo_fraction=halo_star_fraction(m,hmod,cosm)
    ELSE IF(itype==4) THEN
       halo_fraction=halo_boundgas_fraction(m,hmod,cosm)
    ELSE IF(itype==5) THEN
       halo_fraction=halo_freegas_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_FRACTION: Error, itype not specified correcntly'
    END IF

  END FUNCTION halo_fraction

  FUNCTION halo_CDM_fraction(m,hmod,cosm)

    !Mass fraction of a halo in CDM
    IMPLICIT NONE
    REAL :: halo_CDM_fraction
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap
    
    !To prevent compile-time warning
    crap=m
    crap=hmod%A

    !Always the universal value
    halo_CDM_fraction=cosm%om_c/cosm%om_m

  END FUNCTION halo_CDM_fraction

  FUNCTION halo_gas_fraction(m,hmod,cosm)

    !Mass fraction of a halo in gas
    IMPLICIT NONE
    REAL :: halo_gas_fraction
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    halo_gas_fraction=halo_boundgas_fraction(m,hmod,cosm)+halo_freegas_fraction(m,hmod,cosm)

  END FUNCTION halo_gas_fraction

  FUNCTION halo_boundgas_fraction(m,hmod,cosm)

    !Fraction of a halo in bound gas
    IMPLICIT NONE
    REAL :: halo_boundgas_fraction
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm    
    REAL :: m0, sigma, beta

    !Set the model
    !1 - Fedeli (2014a) bound gas model
    !2 - Schneider (2015) bound gas
    !3 - Universal baryon fraction
    INTEGER, PARAMETER :: imod=2 !Set the model

    IF(imod==1) THEN
       !From Fedeli (2014a)
       m0=1.e12
       sigma=3.
       IF(m<m0) THEN
          halo_boundgas_fraction=0.
       ELSE
          halo_boundgas_fraction=erf(log10(m/m0)/sigma)*cosm%om_b/cosm%om_m
       END IF
    ELSE IF(imod==2) THEN
       !From Schneider (2015)
       !m0=1.2d14
       M0=hmod%M0
       beta=0.6
       halo_boundgas_fraction=(cosm%om_b/cosm%om_m)/(1.+(M0/m)**beta)
    ELSE IF(imod==3) THEN
       !Universal baryon fraction model (account for stellar contribution)
       halo_boundgas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,hmod,cosm)
    ELSE
       STOP 'HALO_BOUNDGAS_FRACTION: Error, imod_boundfrac not specified correctly'
    END IF

  END FUNCTION halo_boundgas_fraction

  FUNCTION halo_freegas_fraction(m,hmod,cosm)

    !Mass fraction of a halo in free gas
    IMPLICIT NONE
    REAL :: halo_freegas_fraction
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !This is always all the gas that is not bound or in stars
    halo_freegas_fraction=cosm%om_b/cosm%om_m-halo_star_fraction(m,hmod,cosm)-halo_boundgas_fraction(m,hmod,cosm)
    IF(halo_freegas_fraction<0.) halo_freegas_fraction=0.

  END FUNCTION halo_freegas_fraction

  FUNCTION halo_star_fraction(m,hmod,cosm)

    !Mass fraction of a halo in stars
    IMPLICIT NONE
    REAL :: halo_star_fraction
    REAL, INTENT(IN) :: m
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: m0, sigma, A, min, crap

    !Set the model
    !1 - Fedeli (2014)
    !2 - Constant stellar fraction
    !3 - Fedeli (2014) but saturates at high halo mass
    !4 - No stars
    INTEGER, PARAMETER :: imod=3

    crap=cosm%A

    IF(imod==1 .OR. imod==3) THEN
       !Fedeli (2014)
       !A=0.02
       !IF(variation) A=param(5)
       A=hmod%Astar
       m0=hmod%Mstar
       sigma=hmod%sstar
       halo_star_fraction=A*exp(-((log10(m/m0))**2)/(2.*sigma**2))
       IF(imod==3) THEN
          !Suggested by Ian, the relation I have is for the central stellar mass
          !in reality this saturates for high-mass haloes (due to satellite contribution)
          min=0.01
          IF(halo_star_fraction<min .AND. m>m0) halo_star_fraction=min
       END IF
    ELSE IF(imod==2) THEN
       !Constant star fraction
       A=0.005
       halo_star_fraction=A
    ELSE IF(imod==4) THEN
       !No stars (actually, things crash with exactly zero stars)
       halo_star_fraction=1e-4
    ELSE
       STOP 'HALO_STAR_FRACTION: Error, imod_starfrac specified incorrectly'
    END IF

  END FUNCTION halo_star_fraction

  FUNCTION integrate_hmod(a,b,f,hmod,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
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

       !Fix the answer to zero if the integration limits are identical
       integrate_hmod=0.

    ELSE

       !Set the sum variable for the integration
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
             f1=f(a,hmod)
             f2=f(b,hmod)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,hmod)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod

  FUNCTION integrate_hmod_cosm(a,b,f,hmod,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
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

       !Fix the answer to zero if the integration limits are identical
       integrate_hmod_cosm=0.

    ELSE

       !Set the sum variable for the integration
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
             f1=f(a,hmod,cosm)
             f2=f(b,hmod,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,hmod,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod_cosm=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod_cosm

  FUNCTION integrate_hmod_cosm_exp(a,b,f,hmod,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
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

       !Fix the answer to zero if the integration limits are identical
       integrate_hmod_cosm_exp=0.

    ELSE

       !Set the sum variable for the integration
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
             f1=f(exp(a),hmod,cosm)*exp(a)
             f2=f(exp(b),hmod,cosm)*exp(b)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(exp(x),hmod,cosm)*exp(x)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_HMOD: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_hmod_cosm_exp=REAL(sum_new)

    END IF

  END FUNCTION integrate_hmod_cosm_exp

  FUNCTION integrate_scatter(c,dc,ih,k,z,m,rv,hmod,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_scatter
    REAL, INTENT(IN) :: c, dc, acc
    INTEGER, INTENT(IN) :: iorder
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: k, z, m, rv    
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
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

       !Fix the answer to zero if the integration limits are identical
       integrate_scatter=0.

    ELSE

       !Set the sum variable for the integration
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
             f1=scatter_integrand(a,c,dc,ih,k,z,m,rv,hmod,cosm)
             f2=scatter_integrand(b,c,dc,ih,k,z,m,rv,hmod,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=scatter_integrand(x,c,dc,ih,k,z,m,rv,hmod,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_SCATTER: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_SCATTER: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_scatter=REAL(sum_new)

    END IF

  END FUNCTION integrate_scatter

  FUNCTION scatter_integrand(c,mean_c,sigma_lnc,ih,k,z,m,rv,hmod,cosm)

    !Integrand for computing halo profiles with scatter
    IMPLICIT NONE
    REAL :: scatter_integrand
    REAL, INTENT(IN) :: c, mean_c, sigma_lnc
    INTEGER, INTENT(IN) :: ih(2)
    REAL, INTENT(IN) :: k, z, m, rv    
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: wk(2), pc, rs
    INTEGER :: j

    !Halo profiles
    DO j=1,2
       rs=rv/c
       wk(j)=win_type(.FALSE.,ih(j),1,k,z,m,rv,rs,hmod,cosm)
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

    REAL, PARAMETER :: seconds=2. ! How many seconds should each timing test take (approximate)

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
                ntime=CEILING(seconds/(t2-t1))
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
