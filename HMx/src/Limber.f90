MODULE Limber

  USE constants
  USE array_operations
  USE interpolate
  USE cosmology_functions

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: projection
  PUBLIC :: lensing
  PUBLIC :: maxdist
  PUBLIC :: calculate_Cell
  PUBLIC :: xcorr_type
  PUBLIC :: cell_contribution
  PUBLIC :: fill_projection_kernels
  PUBLIC :: get_nz
  PUBLIC :: set_xcorr_type
  PUBLIC :: write_projection_kernels
  PUBLIC :: write_xi
  PUBLIC :: xcorr
  PUBLIC :: calculate_xi
  PUBLIC :: write_cell
  PUBLIC :: k_ell
 
  ! Projection quantities that need to be calculated only once; these relate to the Limber integrals
  TYPE projection    
     REAL, ALLOCATABLE :: X(:), r_X(:)
     INTEGER :: nX
  END TYPE projection

  ! Quantities that are necessary for lensing specifically
  ! Possibly this could usefully be merged with the projection type
  TYPE lensing
     REAL, ALLOCATABLE :: q(:), r_q(:)
     REAL, ALLOCATABLE :: nz(:), z_nz(:)
     INTEGER :: nq, nnz
  END TYPE lensing
  
  ! P(k,a) look-up table parameters
  REAL, PARAMETER :: kmin_pka=1e-4   ! k' value for P(k,a) table; P(k<k',a)=0
  REAL, PARAMETER :: kmax_pka=1e2    ! k' value for P(k,a) table; P(k>k',a)=0

  ! xcorr - C(l) calculation
  REAL, PARAMETER :: kmin_xcorr=1e-3 ! Minimum k
  REAL, PARAMETER :: kmax_xcorr=1e1  ! Maximum k (some halo-model things go hectic for k>10)
  INTEGER, PARAMETER :: nk_xcorr=64  ! Number of k values (used to be 32)
  REAL, PARAMETER :: amin_xcorr=0.1 ! Minimum scale factor (problems with one-halo term if amin is less than 0.1 (CMB lensing?))
  REAL, PARAMETER :: amax_xcorr=1.0 ! Maximum scale factor
  INTEGER, PARAMETER :: na_xcorr=16 ! Number of scale factores
  LOGICAL, PARAMETER :: verbose_Limber=.FALSE. ! Verbosity
  !LOGICAL, PARAMETER :: verbose_cell=.FALSE. ! Verbosity
  LOGICAL, PARAMETER :: verbose_xi=.FALSE. ! Verbosity

  ! Maxdist
  REAL, PARAMETER :: dr_max=0.01 ! Small subtraction from maxdist to prevent numerical issues

  ! Limber integration parameters
  INTEGER, PARAMETER  :: n_Limber=1024 ! Number of samples to take in r for Limber integration  
  REAL, PARAMETER :: lcorr=0.5         ! 1/2 in LoVerde (2008) Limber correction k(r)=(l+1/2)/f_k(r)
  REAL, PARAMETER :: acc_Limber=1e-4   ! Accuracy parameter for Limber integration

  ! n(z)
  REAL, PARAMETER :: zmin_nz=0.  ! Minimum redshift for the analytic n(z) tables
  REAL, PARAMETER :: zmax_nz=2.5 ! Maximum redshift for the analytic n(z) tables
  INTEGER, PARAMETER :: n_nz=128 ! Number of entries in the analytic n(z) tables

  ! General kernel
  REAL, PARAMETER :: rmin_kernel=0.   ! Minimum r for table
  INTEGER, PARAMETER :: nx_kernel=128 ! Number of entires in look-up table

  ! Lensing kernel/efficiency
  REAL, PARAMETER :: rmin_lensing=0.      ! Minimum distance in integral    
  INTEGER, PARAMETER :: nX_lensing=128    ! Number of entries in X(r) table
  INTEGER, PARAMETER :: nq_efficiency=128 ! Number of entries in q(r) table

  ! Gravitational waves
  REAL, PARAMETER :: A_gwave=1.
  REAL, PARAMETER :: rmin_gwave=10.
  
CONTAINS

  FUNCTION xcorr_type(ix)

    ! Names for cross-correlation field types
    ! TODO: This is super ugly. Surely it should be combined with set_xcorr_type somehow?
    IMPLICIT NONE
    CHARACTER(len=256) :: xcorr_type
    INTEGER, INTENT(IN) :: ix

    xcorr_type=''
    IF(ix==1)  xcorr_type='RCSLenS lensing'
    IF(ix==2)  xcorr_type='Compton y'
    IF(ix==3)  xcorr_type='CMB lensing'
    IF(ix==4)  xcorr_type='CFHTLenS lensing'
    IF(ix==5)  xcorr_type='KiDS lensing (z = 0.1 -> 0.9)'
    IF(ix==6)  xcorr_type='KiDS lensing (z = 0.1 -> 0.3)'
    IF(ix==7)  xcorr_type='KiDS lensing (z = 0.3 -> 0.5)'
    IF(ix==8)  xcorr_type='KiDS lensing (z = 0.5 -> 0.7)'
    IF(ix==9)  xcorr_type='KiDS lensing (z = 0.7 -> 0.9)'
    IF(ix==10) xcorr_type='Gravitational waves'
    IF(ix==11) xcorr_type='KiDS 450 (z = 0.1 -> 0.9)'
    IF(ix==12) xcorr_type='KiDS 450 (z = 0.1 -> 0.5)'
    IF(ix==13) xcorr_type='KiDS 450 (z = 0.5 -> 0.9)'
    IF(ix==14) xcorr_type='KiDS 450 (z = 0.9 -> 3.5)'
    IF(xcorr_type=='') STOP 'XCORR_TYPE: Error, ix not specified correctly'
    
  END FUNCTION xcorr_type

  SUBROUTINE set_xcorr_type(ix,ip)

    ! Set the cross-correlation type
    ! TODO: This is super ugly. Surely it should be combined with xcorr_type somehow?
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    INTEGER, INTENT(OUT) :: ip(2)
    INTEGER :: i, j
    
    INTEGER, PARAMETER :: nx=14 ! Total number of cross-correlation fields

    ! Loop over two-components of xcorr
    DO i=1,2

       IF(ix(i)==-1) THEN
          WRITE(*,fmt='(A30,I3)') 'SET_XCORR_TYPE: Choose field: ', i
          WRITE(*,*) '========================='
          DO j=1,nx
             WRITE(*,fmt='(I3,A3,A30)') j, '- ', TRIM(xcorr_type(j))
          END DO
          READ(*,*) ix(i)
          WRITE(*,*) '========================='
          WRITE(*,*)
       END IF

       IF(ix(i)==2) THEN
          ! Compton y
          ip(i)=6 ! Profile type: 6: Pressure
       ELSE IF(ix(i)==10) THEN
          ! Gravitational waves
          ip(i)=-1 ! Profile type: -1: DMONLY
       ELSE
          ! Gravitational lensing
          ip(i)=0 ! Profile type: 0: Matter
       END IF
       
    END DO

  END SUBROUTINE set_xcorr_type

  SUBROUTINE xcorr(ix,mmin,mmax,ell,Cell,nl,hmod,cosm,verbose)

    ! Calculates the C(l) for the cross correlation of fields ix(1) and ix(2)
    ! Public-facing function
    ! TODO: Remove this explicit HMx dependence, it cannot be necessary, it is the only place it appears
    ! TODO: Maybe this needs to be moved anyway
    USE HMx 
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    REAL, INTENT(IN) :: mmin, mmax
    REAL, INTENT(IN) :: ell(nl)
    REAL, INTENT(OUT) :: Cell(nl)
    INTEGER, INTENT(IN) :: nl
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    REAL, ALLOCATABLE :: a(:), k(:), powa_lin(:,:), powa_2h(:,:), powa_1h(:,:), powa(:,:)
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

    ! Allocate power arrays
    ALLOCATE(powa(nk,na),powa_lin(nk,na),powa_2h(nk,na),powa_1h(nk,na))

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
    CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa,hmod,cosm,verbose,response=.FALSE.)

    ! Fill out the projection kernels
    CALL fill_projection_kernels(ix,proj,cosm)
    IF(verbose) CALL write_projection_kernels(proj,cosm)

    ! Set the range in comoving distance for the Limber integral
    r1=0.
    r2=maxdist(proj)

    ! Actually calculate the C(ell), but only for the full halo model part
    CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)

    ! Write to screen
    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE xcorr

  REAL FUNCTION k_ell(ell,a,cosm)

    ! Finds the k that corresponds to ell at the given a
    IMPLICIT NONE
    REAL, INTENT(IN) :: ell ! angular wave number
    REAL, INTENT(IN) :: a   ! scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm ! cosmology
    REAL :: r

    IF(a==1.) THEN
       ! This should really be infinite
       ! Stops a division by infinity
       k_ell=1e3
    ELSE
       r=comoving_distance(a,cosm)
       k_ell=(ell+lcorr)/f_k(r,cosm)
    END IF
    
  END FUNCTION k_ell

  SUBROUTINE write_nz(lens,output)

    ! Write out the n(z) to a file
    IMPLICIT NONE
    TYPE(lensing), INTENT(IN) :: lens
    CHARACTER(len=*), INTENT(IN) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,lens%nnz
       WRITE(7,*) lens%z_nz(i), lens%nz(i)
    END DO
    CLOSE(7)

  END SUBROUTINE write_nz

  SUBROUTINE fill_projection_kernels(ix,proj,cosm)

    ! Fill look-up tables for the two projection kerels X_ij
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix(2) ! Label for the type of projection kernel needed
    TYPE(projection), INTENT(OUT) :: proj(2) ! Output the projection kernel
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmological model
    INTEGER :: nk, i

    IF(ix(1)==ix(2)) THEN
       nk=1
    ELSE
       nk=2
    END IF

    ! Loop over the two kernels
    DO i=1,nk
       CALL fill_projection_kernel(ix(i),proj(i),cosm)
    END DO

    ! In case the autospectrum is being considered
    IF(nk==1) THEN
       proj(2)=proj(1)
    END IF

  END SUBROUTINE fill_projection_kernels

  SUBROUTINE fill_projection_kernel(ix,proj,cosm)

    ! Fill look-up table for a projection kernel X_i
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(projection), INTENT(OUT) :: proj
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(lensing) :: lens

    IF(ix==2 .OR. ix==10) THEN
       CALL fill_kernel(ix,proj,cosm)
    ELSE
       CALL fill_lensing_kernel(ix,proj,lens,cosm)
    END IF

  END SUBROUTINE fill_projection_kernel

  SUBROUTINE calculate_Cell(r1,r2,ell,Cell,nl,k,a,pow,nk,na,proj,cosm)

    ! Calculates C(l) using the Limber approximation
    ! Note that using Limber and flat-sky for sensible results limits lmin to ell~10
    IMPLICIT NONE
    REAL, INTENT(IN) :: r1, r2 ! Maximum and minimum comoving distances for integration
    REAL, INTENT(IN) :: ell(nl) ! Input array of desired ell values for C(ell)
    REAL, INTENT(OUT) :: Cell(nl) ! Output array for C(ell) values
    INTEGER, INTENT(IN) :: nl ! Number of ell values
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na) ! Input k, a and P(k,a) arrays
    INTEGER, INTENT(IN) :: nk, na ! Number of k and a values
    TYPE(projection), INTENT(IN) :: proj(2) ! Projection kernels for the Limber integration
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: logk(nk), loga(na), logpow(nk,na)
    INTEGER :: i, j

    ! Create log tables to speed up 2D find routine in find_pkz
    logk=log(k)
    loga=log(a)
    DO j=1,na
       logpow(:,j)=log((2.*pi**2)*pow(:,j)/k**3)
    END DO

    ! Write some useful things to the screen
    IF(verbose_Limber) THEN
       WRITE(*,*) 'CALCULATE_CELL: ell min:', REAL(ell(1))
       WRITE(*,*) 'CALCULATE_CELL: ell max:', REAL(ell(nl))
       WRITE(*,*) 'CALCULATE_CELL: number of ell:', nl
       WRITE(*,*) 'CALCULATE_CELL: number of k:', nk
       WRITE(*,*) 'CALCULATE_CELL: number of a:', na
       WRITE(*,*) 'CALCULATE_CELL: Minimum distance [Mpc/h]:', REAL(r1)
       WRITE(*,*) 'CALCULATE_CELL: Maximum distance [Mpc/h]:', REAL(r2)
    END IF

    ! Finally do the integration
    IF(verbose_Limber) WRITE(*,*) 'CALCULATE CELL: Doing calculation'
    DO i=1,nl
       Cell(i)=integrate_Limber(ell(i),r1,r2,logk,loga,logpow,nk,na,acc_Limber,3,proj,cosm)
    END DO
    IF(verbose_Limber) THEN
       WRITE(*,*) 'CALCULATE_CELL: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_Cell

  SUBROUTINE Cell_contribution(r1,r2,k,a,pow,nk,na,proj,cosm)

    ! Calculates the contribution to each ell of C(l) as a function of z, k, r
    ! Note that using Limber and flat-sky for sensible results limits lmin to ~10
    USE string_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: r1, r2 ! Integration range
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na) ! Input arrays of k, a and P(k,a)
    INTEGER, INTENT(IN) :: nk, na ! Number of entries in k and a
    TYPE(projection), INTENT(IN) :: proj(2) ! Projection kernels
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: logk(nk), loga(na), logpow(nk,na)
    INTEGER :: i, j, l
    CHARACTER(len=256) :: fbase, fext, outfile

    INTEGER, PARAMETER :: n=16 ! Number of ell values to take, from l=1 to l=2**(n-1); 2^15 ~ 32,000

    ! Create log tables to speed up 2D find routine in find_pka
    logk=log(k)
    loga=log(a)
    DO j=1,na
       logpow(:,j)=log((2.*pi**2)*pow(:,j)/k**3)
    END DO

    ! Now call the contribution subroutine
    fbase='data/Cell_contrib_ell_'
    fext='.dat'
    DO i=1,n
       l=2**(i-1) ! Set the l
       outfile=number_file(fbase,i,fext)
       CALL Limber_contribution(REAL(l),r1,r2,logk,loga,logpow,nk,na,proj,cosm,outfile)
    END DO

  END SUBROUTINE Cell_contribution

  SUBROUTINE write_Cell(ell,Cell,nl,output)

    ! Write C(l) to a file; writes l, C(l), l(l+1)C(l)/2pi
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nl
    REAL, INTENT(IN) :: ell(nl), Cell(nl)
    CHARACTER(len=256) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,nl
       WRITE(7,*) ell(i), Cell(i), ell(i)*(1.+ell(i))*Cell(i)/twopi
    END DO
    CLOSE(7)

  END SUBROUTINE write_Cell

  SUBROUTINE calculate_xi(th_tab,xi_tab,nth,l_tab,cl_tab,nl,lmax)

    ! Calcuate the correlation functions given a C(ell) table
    USE special_functions
    IMPLICIT NONE
    REAL, INTENT(IN) :: l_tab(nl), cl_tab(nl)
    REAL, INTENT(OUT) :: th_tab(nth), xi_tab(3,nth)
    INTEGER, INTENT(IN) :: nl, lmax, nth
    INTEGER :: i, j
    REAL :: logl(nl), logCl(nl)
    REAL :: theta, Cl, l, xi0, xi2, xi4

    ! Speed up find routine by doing logarithms in advance
    logl=log(l_tab)
    logCl=log(cl_tab)

    IF(verbose_xi) WRITE(*,*) 'CALCULATE_XI: Computing correlation functions via sum'
    DO i=1,nth

       ! Get theta value and convert from degrees to radians
       theta=th_tab(i)/rad2deg

       ! Set values to zero before summing
       xi0=0.
       xi2=0.
       xi4=0.

       ! Do the conversion from Cl to xi as a summation over integer l
       DO j=1,lmax

          l=REAL(j)
          Cl=exp(find(log(l),logl,logCl,nl,3,3,2))

          xi0=xi0+(2.*l+1.)*Cl*Bessel(0,l*theta) ! J0
          xi2=xi2+(2.*l+1.)*Cl*Bessel(2,l*theta) ! J2
          xi4=xi4+(2.*l+1.)*Cl*Bessel(4,l*theta) ! J4

       END DO

       ! Divide by correct pre-factor
       xi0=xi0/(4.*pi)
       xi2=xi2/(4.*pi)
       xi4=xi4/(4.*pi)

       ! Convert theta from radians to degrees
       theta=theta*rad2deg

       ! Populate tables
       th_tab(i)=theta
       xi_tab(1,i)=xi0
       xi_tab(2,i)=xi2
       xi_tab(3,i)=xi4

    END DO
    IF(verbose_xi) THEN
       WRITE(*,*) 'CALCULATE_XI: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_xi

  SUBROUTINE write_xi(th_tab,xi_tab,nth,output)

    IMPLICIT NONE
    REAL, INTENT(IN) :: th_tab(nth), xi_tab(3,nth)
    INTEGER, INTENT(IN) :: nth
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,nth
       WRITE(7,*) th_tab(i), xi_tab(1,i), xi_tab(2,i), xi_tab(3,i)
    END DO
    CLOSE(7)

  END SUBROUTINE write_xi

  FUNCTION maxdist(proj)!,cosm)

    !Calculates the maximum distance necessary for the lensing integration
    IMPLICIT NONE
    REAL :: maxdist
    TYPE(projection), INTENT(IN) :: proj(2)
    !TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rmax1, rmax2

    !Fix the maximum redshift and distance (which may fixed at source plane)
    rmax1=MAXVAL(proj(1)%r_x)
    rmax2=MAXVAL(proj(2)%r_x)
    
    !Subtract a small distance here because of rounding errors in recalculating zmax
    maxdist=MIN(rmax1,rmax2)-dr_max
    
  END FUNCTION maxdist

  SUBROUTINE write_projection_kernels(proj,cosm)

    IMPLICIT NONE
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=256) :: output

    output=TRIM('projection/kernel1.dat')
    CALL write_projection_kernel(proj(1),cosm,output)

    output=TRIM('projection/kernel2.dat')
    CALL write_projection_kernel(proj(2),cosm,output)

  END SUBROUTINE write_projection_kernels

  SUBROUTINE write_projection_kernel(proj,cosm,output)

    IMPLICIT NONE
    TYPE(projection), INTENT(IN) :: proj
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER :: i
    REAL :: r, z
    
    !LOGICAL, PARAMETER :: verbose=.FALSE. ! Verbosity

    !Kernel 1
    IF(verbose_Limber) WRITE(*,*) 'WRITE_PROJECTION_KERNEL: Writing out kernel: ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,proj%nX    
       r=proj%r_X(i)
       z=redshift_r(r,cosm)
       WRITE(7,*) r, z, proj%X(i)
    END DO
    CLOSE(7)
    IF(verbose_Limber) THEN
       WRITE(*,*) 'WRITE_PROJECTION_KERNEL: Writing done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_projection_kernel

  SUBROUTINE fill_lensing_kernel(ix,proj,lens,cosm)

    !Fill the lensing projection kernel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing) :: lens
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(projection), INTENT(OUT) :: proj
    REAL :: zmin, zmax, rmax, amax, r
    CHARACTER(len=256) :: output
    INTEGER :: i, nX

    IF(ix==2 .OR. ix==10) STOP 'FILL_LENSING_KERNEL: Error, trying to do this for a non-lensing ix'

    !Choose either n(z) or fixed z_s
    IF(ix==3) THEN      
       zmin=0.
       zmax=cosm%z_cmb
       IF(verbose_Limber) WRITE(*,*) 'FILL_LENSING_KERNEL: Source plane redshift:', REAL(zmax)
    ELSE
       CALL get_nz(ix,lens)
       zmin=lens%z_nz(1)
       zmax=lens%z_nz(lens%nnz)
       output='data/nz.dat'
       CALL write_nz(lens,output)
    END IF

    !Get the distance range for the lensing kernel
    amax=scale_factor_z(zmax)
    rmax=comoving_distance(amax,cosm)
    IF(verbose_Limber) THEN
       WRITE(*,*) 'FILL_LENSING_KERNEL: minimum r [Mpc/h]:', REAL(rmin_lensing)
       WRITE(*,*) 'FILL_LENSING_KERNEL: maximum r [Mpc/h]:', REAL(rmax)
       WRITE(*,*) 'FILL_LENSING_KERNEL: minimum z:', REAL(zmin)
       WRITE(*,*) 'FILL_LENSING_KERNEL: maximum z:', REAL(zmax)
       WRITE(*,*)
    END IF

    !Fill the q(r) table
    CALL fill_lensing_efficiency(ix,rmin_lensing,rmax,zmax,lens,cosm)

    !Write the q(r) table to file
    output='data/lensing_efficiency.dat'
    CALL write_lensing_efficiency(lens,cosm,output)

    !Assign arrays for the projection function
    nX=nX_lensing
    proj%nx=nx
    CALL fill_array(rmin_lensing,rmax,proj%r_x,nx)
    IF(verbose_Limber) WRITE(*,*) 'FILL_LENSING_KERNEL: number of points:', nx
    IF(ALLOCATED(proj%x)) DEALLOCATE(proj%x)
    ALLOCATE(proj%x(proj%nx))

    DO i=1,nx
       !Get r and fill X(r)
       r=proj%r_x(i)
       proj%x(i)=lensing_kernel(r,lens,cosm)  
       !Enforce that the kernel must not be negative (is this necessary?)
       IF(proj%x(i)<0.) proj%x(i)=0.
    END DO
    IF(verbose_Limber) THEN
       WRITE(*,*) 'FILL_LENSING_KERNEL: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE fill_lensing_kernel

  SUBROUTINE fill_lensing_efficiency(ix,rmin,rmax,zmax,lens,cosm)

    !Fill the table for lensing efficiency
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    REAL, INTENT(IN) :: rmin, rmax, zmax
    TYPE(lensing), INTENT(INOUT) :: lens
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r, z
    INTEGER :: i

    !Fill the r vs. q(r) tables
    lens%nq=nq_efficiency
    IF(verbose_Limber) WRITE(*,*) 'FILL_EFFICIENCY: number of points:', lens%nq
    CALL fill_array(rmin,rmax,lens%r_q,lens%nq)
    IF(ALLOCATED(lens%q)) DEALLOCATE(lens%q)
    ALLOCATE(lens%q(lens%nq))

    DO i=1,lens%nq       
       r=lens%r_q(i)
       z=redshift_r(r,cosm)
       IF(r==0.) THEN
          !To avoid division by zero
          lens%q(i)=1.
       ELSE
          IF(ix==3) THEN
             !q(r) for a fixed source plane
             lens%q(i)=f_k(rmax-r,cosm)/f_k(rmax,cosm)
          ELSE
             !q(r) for a n(z) distribution 
             lens%q(i)=integrate_q(r,z,zmax,acc_Limber,3,lens,cosm)
          END IF
       END IF
    END DO
    IF(verbose_Limber) THEN
       WRITE(*,*) 'FILL_EFFICIENCY: Done writing'
       WRITE(*,*)
    END IF

  END SUBROUTINE fill_lensing_efficiency

  FUNCTION q_r(r,lens)

    !Interpolation function for q(r)
    IMPLICIT NONE
    REAL :: q_r
    REAL, INTENT(IN) :: r
    TYPE(lensing), INTENT(IN) :: lens

    q_r=find(r,lens%r_q,lens%q,lens%nq,3,3,2)

  END FUNCTION q_r

  SUBROUTINE write_lensing_efficiency(lens,cosm,output)

    !Write lensing efficiency q(r) function to a file
    IMPLICIT NONE
    TYPE(lensing), INTENT(INOUT) :: lens
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=256), INTENT(IN) :: output
    REAL :: r, z, q
    INTEGER :: i

    WRITE(*,*) 'WRITE_EFFICIENCY: Writing q(r): ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,lens%nq
       r=lens%r_q(i)
       z=redshift_r(r,cosm)
       q=lens%q(i)
       WRITE(7,*) r, z, q
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_EFFICIENCY: Done'
    WRITE(*,*)

  END SUBROUTINE write_lensing_efficiency

  FUNCTION lensing_kernel(r,lens,cosm)

    !The lensing projection kernel
    IMPLICIT NONE
    REAL :: lensing_kernel
    REAL, INTENT(IN) :: r
    TYPE(lensing) :: lens
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: z, q

    !Get z(r)
    z=redshift_r(r,cosm)

    !Get the lensing efficiency
    q=q_r(r,lens)

    !This is then the projection kernel (X_kappa)
    lensing_kernel=(1.+z)*f_k(r,cosm)*q
    lensing_kernel=lensing_kernel*1.5*cosm%om_m/(Hdist**2)

  END FUNCTION lensing_kernel

  SUBROUTINE fill_kernel(ix,proj,cosm)

    ! Fill a table of kernel values
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(projection), INTENT(OUT) :: proj
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, nX
    REAL :: rmax, r

    nX=nx_kernel
    proj%nx=nx

    ! Get the distance range for the projection function
    ! Use the same as that for the distance calculation
    ! Assign arrays for the kernel function
    rmax=MAXVAL(cosm%r)
    CALL fill_array(rmin_kernel,rmax,proj%r_X,proj%nX)
    WRITE(*,*) 'FILL_KERNEL: minimum r [Mpc/h]:', REAL(rmin_kernel)
    WRITE(*,*) 'FILL_KERNEL: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'FILL_KERNEL: number of points:', nX

    IF(ALLOCATED(proj%X)) DEALLOCATE(proj%X)
    ALLOCATE(proj%X(nX))

    ! Now fill the kernels
    DO i=1,nX
       r=proj%r_x(i)
       IF(ix==2) THEN
          proj%x(i)=y_kernel(r,cosm)
       ELSE IF(ix==10) THEN
          proj%x(i)=gwave_kernel(r,cosm)
       ELSE
          STOP 'FILL_KERNEL: Error, ix not specified correctly'
       END IF
    END DO

    WRITE(*,*) 'FILL_KERNEL: Done:'
    WRITE(*,*)

  END SUBROUTINE fill_kernel

  FUNCTION y_kernel(r,cosm)

    ! The Compton-y projection kernel
    IMPLICIT NONE
    REAL :: y_kernel
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: z, a

    ! Get the scale factor
    z=redshift_r(r,cosm)
    a=scale_factor_z(z)

    ! Make the kernel and do some unit conversions
    y_kernel=yfac                     ! yfac = sigma_T / m_e c^2 [kg^-1 s^2]
    y_kernel=y_kernel*Mpc/cosm%h      ! NEW: Add Mpc/h units from the dr in the integral (h is new)
    y_kernel=y_kernel/a**2            ! NEW: These come from 'a^-3' for pressure multiplied by 'a' for comoving distance
    y_kernel=y_kernel*eV*(0.01)**(-3) ! Convert units of pressure spectrum from [eV/cm^3] to [J/m^3]

  END FUNCTION y_kernel

  REAL FUNCTION gwave_kernel(r,cosm)

    ! Projection kernel for gravitational waves (~ 1/r)
    IMPLICIT NONE
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    ! To stop compile-time warnings
    crap=cosm%Om_m

    IF(r<rmin_gwave) THEN
       gwave_kernel=A_gwave/rmin_gwave
    ELSE
       gwave_kernel=A_gwave/r
    END IF

  END FUNCTION gwave_kernel

  SUBROUTINE get_nz(ix,lens)

    ! The the n(z) function for lensing
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens

    IF(ix==1 .OR. ix==4) THEN
       CALL fill_analytic_nz_table(ix,lens)
    ELSE
       CALL fill_nz_table(ix,lens)
    END IF

    IF(verbose_Limber) THEN
       WRITE(*,*) 'GET_NZ: zmin:', lens%z_nz(1)
       WRITE(*,*) 'GET_NZ: zmax:', lens%z_nz(lens%nnz)
       WRITE(*,*) 'GET_NZ: nz:', lens%nnz
       WRITE(*,*)
    END IF

  END SUBROUTINE get_nz

  SUBROUTINE fill_analytic_nz_table(ix,lens)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens
    INTEGER :: i

    ! From analytical function
    lens%nnz=n_nz
    IF(ALLOCATED(lens%z_nz)) DEALLOCATE(lens%z_nz)
    IF(ALLOCATED(lens%nz))   DEALLOCATE(lens%nz)
    ALLOCATE(lens%z_nz(lens%nnz),lens%nz(lens%nnz))

    ! Fill the look-up tables
    CALL fill_array(zmin_nz,zmax_nz,lens%z_nz,lens%nnz)
    DO i=1,n_nz
       lens%nz(i)=nz_lensing(lens%z_nz(i),ix)
    END DO

  END SUBROUTINE fill_analytic_nz_table

  SUBROUTINE fill_nz_table(ix,lens)

    USE file_info
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens
    INTEGER :: i
    REAL :: spam
    CHARACTER(len=256) :: input

    ! Get file name
    IF(ix==5) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.1-0.9_MEAD.txt'
    ELSE IF(ix==6) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.1-0.3.txt'
    ELSE IF(ix==7) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.3-0.5.txt'
    ELSE IF(ix==8) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.5-0.7.txt'
    ELSE IF(ix==9) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS_z0.7-0.9.txt'
    ELSE IF(ix==11 .OR. ix==12 .OR. ix==13 .OR. ix==14) THEN
       input='/Users/Mead/Physics/KiDS/nz/KiDS-450_fat_bin_nofz.txt'
    ELSE
       STOP 'GET_NZ: ix not specified correctly'
    END IF
    WRITE(*,*) 'GET_NZ: Input file:', TRIM(input)

    ! Allocate arrays
    lens%nnz=count_number_of_lines(input)
    IF(ALLOCATED(lens%z_nz)) DEALLOCATE(lens%z_nz)
    IF(ALLOCATED(lens%nz))   DEALLOCATE(lens%nz)
    ALLOCATE(lens%z_nz(lens%nnz),lens%nz(lens%nnz))

    ! Read in n(z) table
    OPEN(7,file=input)
    DO i=1,lens%nnz
       IF(ix==5 .OR. ix==6 .OR. ix==7 .OR. ix==8 .OR. ix==9) THEN
          READ(7,*) lens%z_nz(i), lens%nz(i) ! Second column
       ELSE IF(ix==11) THEN
          READ(7,*) lens%z_nz(i), lens%nz(i) ! Second column (z = 0.1 -> 0.9)
       ELSE IF(ix==12) THEN
          READ(7,*) lens%z_nz(i), spam, lens%nz(i) ! Third column (z = 0.1 -> 0.5)
       ELSE IF(ix==13) THEN
          READ(7,*) lens%z_nz(i), spam, spam, lens%nz(i) ! Fourth column (z = 0.5 -> 0.9)
       ELSE IF(ix==14) THEN
          READ(7,*) lens%z_nz(i), spam, spam, spam, lens%nz(i) ! Fifth column (z = 0.9 -> 3.5)
       ELSE
          STOP 'GET_NZ: ix not specified correctly'
       END IF
    END DO
    CLOSE(7)

    ! Do this because the KiDS-450 files contain the lower left edge of histograms
    ! The bin sizes are 0.05 in z, so need to add 0.05/2 = 0.025
    IF(ix==11 .OR. ix==12 .OR. ix==13 .OR. ix==14) THEN       
       lens%z_nz=lens%z_nz+0.025
    END IF

  END SUBROUTINE fill_nz_table

  FUNCTION nz_lensing(z,ix)

    ! Analytical n(z) for different surveys
    IMPLICIT NONE
    REAL :: nz_lensing
    REAL, INTENT(IN) :: z
    INTEGER, INTENT(IN) :: ix
    REAL :: a, b, c, d, e, f, g, h, i
    REAL :: n1, n2, n3
    REAL :: z1, z2

    IF(ix==1) THEN
       ! RCSLenS
       a=2.94
       b=-0.44
       c=1.03
       d=1.58
       e=0.40
       f=0.25
       g=0.38
       h=0.81
       i=0.12
       n1=a*z*exp(-(z-b)**2/c**2)
       n2=d*z*exp(-(z-e)**2/f**2)
       n3=g*z*exp(-(z-h)**2/i**2)
       nz_lensing=n1+n2+n3
    ELSE IF(ix==4) THEN
       ! CFHTLenS
       z1=0.7 ! Not a free parameter in Van Waerbeke et al. (2013)
       z2=1.2 ! Not a free parameter in Van Waerbeke et al. (2013)
       a=1.50
       b=0.32
       c=0.20
       d=0.46
       nz_lensing=a*exp(-((z-z1)/b)**2)+c*exp(-((z-z2)/d)**2)
    ELSE
       STOP 'NZ_LENSING: ix specified incorrectly'
    END IF

  END FUNCTION nz_lensing

  FUNCTION integrate_q(r,a,b,acc,iorder,lens,cosm)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_q
    REAL, INTENT(IN) :: a, b, r, acc
    INTEGER, INTENT(IN) :: iorder
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(lensing), INTENT(IN) :: lens
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5  ! Standard integration parameters
    INTEGER, PARAMETER :: jmax=30 ! Standard integration parameters

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_q=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=q_integrand(a,r,lens,cosm)
             f2=q_integrand(b,r,lens,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=q_integrand(x,r,lens,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_Q: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             ! jmin avoids spurious early convergence
             integrate_q=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             integrate_q=0.d0
             STOP 'INTEGRATE_Q: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             integrate_q=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_q

  FUNCTION q_integrand(z,r,lens,cosm)

    ! The lensing efficiency integrand, which is a function of z
    ! z is integrated over while r is just a parameter
    ! This is only called for n(z)
    IMPLICIT NONE
    REAL :: q_integrand
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(lensing), INTENT(IN) :: lens
    REAL :: rdash, nz, a

    a=scale_factor_z(z)

    IF(z==0.) THEN
       q_integrand=0.
    ELSE
       ! Find the r'(z) variable that is integrated over     
       rdash=comoving_distance(a,cosm)
       ! Find the n(z)
       nz=find(z,lens%z_nz,lens%nz,lens%nnz,3,3,2)
       ! This is then the integrand
       q_integrand=nz*f_k(rdash-r,cosm)/f_k(rdash,cosm)
    END IF

  END FUNCTION q_integrand

  FUNCTION integrate_Limber(l,a,b,logktab,logatab,logptab,nk,na,acc,iorder,proj,cosm)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_Limber
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    REAL, INTENT(IN) :: l
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)
    INTEGER, INTENT(IN) :: nk, na
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5  ! Standard integration parameters
    INTEGER, PARAMETER :: jmax=25 ! Standard integration parameters

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate_Limber=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=Limber_integrand(a,l,logktab,logatab,logptab,nk,na,proj,cosm)
             f2=Limber_integrand(b,l,logktab,logatab,logptab,nk,na,proj,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=Limber_integrand(x,l,logktab,logatab,logptab,nk,na,proj,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE_LIMBER: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             ! jmin avoids spurious early convergence
             integrate_Limber=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             integrate_Limber=0.d0
             STOP 'INTEGRATE_LIMBER: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             integrate_Limber=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_Limber

  FUNCTION Limber_integrand(r,l,logktab,logatab,logptab,nk,na,proj,cosm)

    ! The integrand for the Limber integral
    IMPLICIT NONE
    REAL :: Limber_integrand
    REAL, INTENT(IN) :: r, l
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)
    INTEGER, INTENT(IN) :: nk, na
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(projection), INTENT(IN) :: proj(2)
    REAL :: z, a, k, X(2)
    INTEGER :: i

    IF(r==0.) THEN

       Limber_integrand=0.

    ELSE

       ! Get the two kernels
       DO i=1,2
          X(i)=find(r,proj(i)%r_X,proj(i)%X,proj(i)%nX,3,3,2)
          !X(i)=exp(find(log(r),log(proj(i)%r_X),log(proj(i)%X),proj(i)%nX,3,3,2)) ! Barfed with this
       END DO

       ! Get variables r, z(r) and k(r) for P(k,z)
       z=redshift_r(r,cosm)
       a=scale_factor_z(z)
       k=(l+lcorr)/f_k(r,cosm) ! LoVerde et al. (2008) Limber correction

       ! Construct the integrand
       Limber_integrand=X(1)*X(2)*find_pka(k,a,logktab,logatab,logptab,nk,na)/f_k(r,cosm)**2

    END IF

  END FUNCTION Limber_integrand

  SUBROUTINE Limber_contribution(l,r1,r2,logktab,logatab,logptab,nk,na,proj,cosm,outfile)

    USE array_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: l, r1, r2
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)
    INTEGER, INTENT(IN) :: nk, na
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(projection), INTENT(IN) :: proj(2)
    CHARACTER(len=256), INTENT(IN) :: outfile
    REAL :: k, r, z, total, int, a
    INTEGER :: i

    ! Calculate the integral for this value of ell
    total=integrate_Limber(l,r1,r2,logktab,logatab,logptab,nk,na,acc_Limber,3,proj,cosm)

    ! Now split up the contributions
    ! You need the Jacobian and to remember that the contribution is split in ln(k), ln(z) and ln(R)
    ! This means factors of:
    ! k - f_k(r)/f'_k(r) (=r if flat)
    ! z - z/H(z)
    ! r - r
    OPEN(7,file=TRIM(outfile))
    DO i=1,n_Limber
       r=progression(r1,r2,i,n_Limber)
       IF(r==0.) THEN
          ! Avoid trouble for r = 0 exactly
          ! Maybe should do some sort of awful Taylor expansion here
          CYCLE
       ELSE
          k=(l+lcorr)/f_k(r,cosm)
          z=redshift_r(r,cosm)
          a=scale_factor_z(z)
          int=Limber_integrand(r,l,logktab,logatab,logptab,nk,na,proj,cosm)
          WRITE(7,*) k, int*f_k(r,cosm)/(fdash_k(r,cosm)*total), z, int*z/(sqrt(Hubble2(a,cosm))*total), r, int*r/total
       END IF
    END DO
    CLOSE(7)

  END SUBROUTINE Limber_contribution

  REAL FUNCTION find_pka(k,a,logktab,logatab,logptab,nk,na)

    ! Looks up the power as a 2D function of k and a
    ! Note that it cuts P(k,a) off above and below certain wavenumbers defined in the header (kmin_pka, kmax_pka)
    ! It will interpolate in log(a) outside range
    ! It will interpolate in log(k) outside range of ktab until kmin_pka/kmax_pka
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, a ! Input desired values of k and a
    INTEGER, INTENT(IN) :: nk, na ! Number of entried of k and a in arrays
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na) ! Arrays of log(k), log(a) and log(P(k,a))

    ! Get the power
    IF(k<kmin_pka .OR. k>kmax_pka) THEN
       find_pka=0.
    ELSE
       find_pka=exp(find2d(log(k),logktab,log(a),logatab,logptab,nk,na,3,3,1))
    END IF

  END FUNCTION find_pka
  
END MODULE Limber
