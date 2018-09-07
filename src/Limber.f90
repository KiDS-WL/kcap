MODULE Limber

  USE constants
  USE array_operations
  USE interpolate
  USE cosmology_functions

  IMPLICIT NONE

  REAL, PARAMETER :: lcorr=0.5
  REAL, PARAMETER :: acc_Limber=1e-4

  !Projection quantities that need to be calculated only once
  !These relate to the Limber integrals
  TYPE projection    
     REAL, ALLOCATABLE :: X(:), r_X(:)
     INTEGER :: nX
  END TYPE projection

  !Quantities that are necessary for lensing specifically
  !Possibly this could usefully be merged with the projection type
  TYPE lensing
     REAL, ALLOCATABLE :: q(:), r_q(:)
     REAL, ALLOCATABLE :: nz(:), z_nz(:)
     INTEGER :: nq, nnz
  END TYPE lensing
  
CONTAINS

  FUNCTION xcorr_type(ix)

    !Names for cross-correlation field types
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
    IF(xcorr_type=='') STOP 'XCORR_TYPE: Error, ix not specified correctly'
    
  END FUNCTION xcorr_type

  SUBROUTINE set_xcorr_type(ix,ip)

    !Set the cross-correlation type
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    INTEGER, INTENT(OUT) :: ip(2)
    INTEGER :: i, j

    DO i=1,2

       IF(ix(i)==-1) THEN
          WRITE(*,fmt='(A30,I3)') 'SET_XCORR_TYPE: Choose field: ', i
          WRITE(*,*) '========================='
          DO j=1,10
             WRITE(*,fmt='(I3,A3,A30)') j, '- ', TRIM(xcorr_type(j))
          END DO
          READ(*,*) ix(i)
          WRITE(*,*) '========================='
          WRITE(*,*)
       END IF

       IF(ix(i)==2) THEN
          !Compton y
          ip(i)=6 !Profile type: 6 - Pressure
       ELSE IF(ix(i)==10) THEN
          !Gravitational waves
          ip(i)=-1 !Profile type: -1 DMONLY
       ELSE
          !Gravitational lensing (should be set to 0 eventually)
          ip(i)=-1 !Profile type: 0 - Matter
       END IF
       
    END DO

  END SUBROUTINE set_xcorr_type

  SUBROUTINE xcorr(ix,mmin,mmax,ell,Cell,nl,hmod,cosm,verbose)

    !Calculates the C(l) for the cross correlation of fields ix(1) and ix(2)
    USE HMx !TODO: Remove this explicit HMx dependence, it cannot be necessary
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ix(2)
    INTEGER, INTENT(IN) :: nl
    REAL, INTENT(IN) :: ell(nl), mmin, mmax
    REAL, INTENT(OUT) :: Cell(nl)
    TYPE(halomod), INTENT(INOUT) :: hmod
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, ALLOCATABLE :: a(:), k(:), powa_lin(:,:), powa_2h(:,:), powa_1h(:,:), powa(:,:)
    !TYPE(tables) :: lut
    !TYPE(lensing) :: lens
    TYPE(projection) :: proj(2)
    LOGICAL, INTENT(IN) :: verbose
    REAL :: kmin, kmax, amin, amax, lmin, lmax
    INTEGER :: nk, na, ip(2)
    REAL :: r1, r2

    !Set the k range
    kmin=1e-3
    kmax=1e1
    nk=32
    CALL fill_array(log(kmin),log(kmax),k,nk)
    k=exp(k)   

    !Set the a range
    amin=0.1 !scale_factor(cosm%z_cmb) !Problems with one-halo term if amin is less than 0.1
    amax=1.
    na=16
    CALL fill_array(amin,amax,a,na)

    lmin=ell(1)
    lmax=ell(nl)

    !Allocate power arrays
    ALLOCATE(powa(nk,na),powa_lin(nk,na),powa_2h(nk,na),powa_1h(nk,na))

    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Cross-correlation information'
       WRITE(*,*) 'XCORR: P(k) minimum k [h/Mpc]:', REAL(kmin)
       WRITE(*,*) 'XCORR: P(k) maximum k [h/Mpc]:', REAL(kmax)
       WRITE(*,*) 'XCORR: Number of k:', nk
       WRITE(*,*) 'XCORR: minimum a:', REAL(amin)
       WRITE(*,*) 'XCORR: maximum a:', REAL(amax)
       WRITE(*,*) 'XCORR: number of a:', na
       WRITE(*,*) 'XCORR: minimum ell:', REAL(lmin)
       WRITE(*,*) 'XCORR: maximum ell:', REAL(lmax)
       WRITE(*,*) 'XCORR: number of ell:', nl
       WRITE(*,*)
    END IF

    CALL set_xcorr_type(ix,ip)

    CALL calculate_HMx(ip,mmin,mmax,k,nk,a,na,powa_lin,powa_2h,powa_1h,powa,hmod,cosm,verbose,response=.FALSE.)

    !Fill out the projection kernels
    CALL fill_projection_kernels(ix,proj,cosm)
    IF(verbose) CALL write_projection_kernels(proj,cosm)

    !Set the distance range for the Limber integral
    r1=0.
    r2=maxdist(proj)

    !Actually calculate the C(ell), but only for the full halo model part
    CALL calculate_Cell(r1,r2,ell,Cell,nl,k,a,powa,nk,na,proj,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'XCORR: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE xcorr

  SUBROUTINE write_nz(lens,output)

    IMPLICIT NONE
    TYPE(lensing), INTENT(IN) :: lens
    CHARACTER(len=256), INTENT(IN) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,lens%nnz
       WRITE(7,*) lens%z_nz(i), lens%nz(i)
    END DO
    CLOSE(7)

  END SUBROUTINE write_nz

  SUBROUTINE fill_projection_kernels(ix,proj,cosm)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    TYPE(projection) :: proj(2)
    INTEGER :: nk, i

    IF(ix(1)==ix(2)) THEN
       nk=1
    ELSE
       nk=2
    END IF

    !Loop over the two kernels
    DO i=1,nk
       CALL fill_projection_kernel(ix(i),proj(i),cosm)
    END DO

    !In case the autospectrum is being considered
    IF(nk==1) THEN
       proj(2)=proj(1)
    END IF

  END SUBROUTINE fill_projection_kernels

  SUBROUTINE fill_projection_kernel(ix,proj,cosm)

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

    !Calculates C(l) using the Limber approximation
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nl, nk, na
    REAL, INTENT(IN) :: ell(nl)
    REAL, INTENT(OUT) :: Cell(nl)
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    REAL, INTENT(IN) :: r1, r2
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: logk(nk), loga(na), logpow(nk,na)
    INTEGER :: i, j

    LOGICAL, PARAMETER :: verbose=.FALSE.

    !Note that using Limber and flat-sky for sensible results limits lmin to ell~10

    !Create log tables to speed up 2D find routine in find_pkz
    logk=log(k)
    loga=log(a)
    DO j=1,na
       logpow(:,j)=log((2.*pi**2)*pow(:,j)/k**3)
    END DO

    !Write some useful things to the screen
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_CELL: ell min:', REAL(ell(1))
       WRITE(*,*) 'CALCULATE_CELL: ell max:', REAL(ell(nl))
       WRITE(*,*) 'CALCULATE_CELL: number of ell:', nl
       WRITE(*,*) 'CALCULATE_CELL: number of k:', nk
       WRITE(*,*) 'CALCULATE_CELL: number of a:', na
       WRITE(*,*) 'CALCULATE_CELL: Minimum distance [Mpc/h]:', REAL(r1)
       WRITE(*,*) 'CALCULATE_CELL: Maximum distance [Mpc/h]:', REAL(r2)
    END IF

    !Finally do the integration
    IF(verbose) WRITE(*,*) 'CALCULATE CELL: Doing calculation'
    DO i=1,nl
       Cell(i)=integrate_Limber(ell(i),r1,r2,logk,loga,logpow,nk,na,acc_Limber,3,proj,cosm)
    END DO
    IF(verbose) THEN
       WRITE(*,*) 'CALCULATE_CELL: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE calculate_Cell

  SUBROUTINE Cell_contribution(r1,r2,k,a,pow,nk,na,proj,cosm)

    !Calculates C(l) using the Limber approximation
    USE string_operations
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k(nk), a(na), pow(nk,na)
    REAL, INTENT(IN) :: r1, r2
    TYPE(projection), INTENT(IN) :: proj(2)
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: logk(nk), loga(na), logpow(nk,na)
    INTEGER :: i, j, l
    CHARACTER(len=256) :: fbase, fext, outfile

    INTEGER, PARAMETER :: n=16 !Number of ell values to take, from ell=1 to ell=2**(n-1)

    !Note that using Limber and flat-sky for sensible results limits lmin to ~10

    !Create log tables to speed up 2D find routine in find_pkz
    logk=log(k)
    loga=log(a)
    DO j=1,na
       logpow(:,j)=log((2.*pi**2)*pow(:,j)/k**3)
    END DO

    !Now call the contribution subroutine
    fbase='Limber/Cell_contrib_ell_'
    fext='.dat'
    DO i=1,n
       l=2**(i-1)
       outfile=number_file(fbase,i,fext)
       CALL Limber_contribution(REAL(l),r1,r2,logk,loga,logpow,nk,na,proj,cosm,outfile)
    END DO

  END SUBROUTINE Cell_contribution

  SUBROUTINE write_Cell(ell,Cell,nl,output)

    !Write C(ell) to a file
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nl
    REAL, INTENT(IN) :: ell(nl), Cell(nl)
    CHARACTER(len=256) :: output
    INTEGER :: i

    OPEN(7,file=output)
    DO i=1,nl
       WRITE(7,*) ell(i), Cell(i), ell(i)*(1.+ell(i))*Cell(i)/(2.*pi)
    END DO
    CLOSE(7)

  END SUBROUTINE write_Cell

  SUBROUTINE calculate_xi(th_tab,xi_tab,nth,l_tab,cl_tab,nl,lmax)

    !Calcuate the correlation functions given a C(ell) table
    USE special_functions
    IMPLICIT NONE
    REAL, INTENT(IN) :: l_tab(nl), cl_tab(nl)
    REAL, INTENT(OUT) :: th_tab(nth), xi_tab(3,nth)
    INTEGER, INTENT(IN) :: nl, lmax, nth
    INTEGER :: i, j
    REAL :: logl(nl), logCl(nl)
    REAL :: theta, Cl, l, xi0, xi2, xi4

    !Speed up find routine by doing logarithms in advance
    logl=log(l_tab)
    logCl=log(cl_tab)

    !WRITE(*,*) 'CALCULATE_XI: Computing correlation functions via sum'
    DO i=1,nth

       !Get theta value and convert from degrees to radians
       theta=th_tab(i)/rad2deg

       !Set values to zero before summing
       xi0=0.
       xi2=0.
       xi4=0.

       !Do the conversion from Cl to xi as a summation over integer ell
       DO j=1,lmax

          l=REAL(j)
          Cl=exp(find(log(l),logl,logCl,nl,3,3,2))

          xi0=xi0+(2.*l+1.)*Cl*Bessel(0,l*theta)
          xi2=xi2+(2.*l+1.)*Cl*Bessel(2,l*theta)
          xi4=xi4+(2.*l+1.)*Cl*Bessel(4,l*theta)

       END DO

       !Divide by correct pre-factor
       xi0=xi0/(4.*pi)
       xi2=xi2/(4.*pi)
       xi4=xi4/(4.*pi)

       !Convert theta from radians to degrees
       theta=theta*rad2deg

       !Populate tables
       th_tab(i)=theta
       xi_tab(1,i)=xi0
       xi_tab(2,i)=xi2
       xi_tab(3,i)=xi4

    END DO
    !WRITE(*,*) 'CALCULATE_XI: Done'
    !WRITE(*,*)

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

    REAL, PARAMETER :: dr=0.01

    !Fix the maximum redshift and distance (which may fixed at source plane)
    rmax1=MAXVAL(proj(1)%r_x)
    rmax2=MAXVAL(proj(2)%r_x)
    
    !Subtract a small distance here because of rounding errors in recalculating zmax
    maxdist=MIN(rmax1,rmax2)-dr
    
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
    
    LOGICAL, PARAMETER :: verbose=.FALSE.

    !Kernel 1
    IF(verbose) WRITE(*,*) 'WRITE_PROJECTION_KERNEL: Writing out kernel: ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,proj%nX    
       r=proj%r_X(i)
       z=redshift_r(r,cosm)
       WRITE(7,*) r, z, proj%X(i)
    END DO
    CLOSE(7)
    IF(verbose) THEN
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
    INTEGER :: i

    !Parameters
    REAL, PARAMETER :: rmin=0. !Minimum distance in integral    
    INTEGER, PARAMETER :: nx=128 !Number of entries in X(r) table

    IF(ix==2 .OR. ix==10) STOP 'FILL_LENSING_KERNEL: Error, trying to do this for a non-lensing ix'

    !Choose either n(z) or fixed z_s
    IF(ix==3) THEN      
       zmin=0.
       zmax=cosm%z_cmb
       WRITE(*,*) 'FILL_LENSING_KERNEL: Source plane redshift:', REAL(zmax)
    ELSE
       CALL get_nz(ix,lens)
       zmin=lens%z_nz(1)
       zmax=lens%z_nz(lens%nnz)
       output='lensing/nz.dat'
       CALL write_nz(lens,output)
    END IF

    !Get the distance range for the lensing kernel
    amax=scale_factor_z(zmax)
    rmax=comoving_distance(amax,cosm)
    WRITE(*,*) 'FILL_LENSING_KERNEL: minimum r [Mpc/h]:', REAL(rmin)
    WRITE(*,*) 'FILL_LENSING_KERNEL: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'FILL_LENSING_KERNEL: minimum z:', REAL(zmin)
    WRITE(*,*) 'FILL_LENSING_KERNEL: maximum z:', REAL(zmax)
    WRITE(*,*)

    !Fill the q(r) table
    CALL fill_lensing_efficiency(ix,rmin,rmax,zmax,lens,cosm)

    !Write the q(r) table to file
    output='lensing/efficiency.dat'
    CALL write_lensing_efficiency(lens,cosm,output)

    !Assign arrays for the projection function
    proj%nx=nx
    CALL fill_array(rmin,rmax,proj%r_x,nx)
    WRITE(*,*) 'FILL_LENSING_KERNEL: number of points:', nx
    IF(ALLOCATED(proj%x)) DEALLOCATE(proj%x)
    ALLOCATE(proj%x(proj%nx))

    DO i=1,nx
       !Get r and fill X(r)
       r=proj%r_x(i)
       proj%x(i)=lensing_kernel(r,lens,cosm)  
       !Enforce that the kernel must not be negative (is this necessary?)
       IF(proj%x(i)<0.) proj%x(i)=0.
    END DO
    WRITE(*,*) 'FILL_LENSING_KERNEL: Done'
    WRITE(*,*)

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

    INTEGER, PARAMETER :: nq=128 !Number of entries in q(r) table

    !Fill the r vs. q(r) tables
    lens%nq=nq
    WRITE(*,*) 'FILL_EFFICIENCY: number of points:', lens%nq
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
    WRITE(*,*) 'FILL_EFFICIENCY: Done writing'
    WRITE(*,*)

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

    !Fill a table of kernel values
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(projection), INTENT(OUT) :: proj
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: rmax, r

    REAL, PARAMETER :: rmin=0. !Minimum r for table
    INTEGER, PARAMETER :: nx=128 !Entires in look-up table

    proj%nx=nx

    !Get the distance range for the projection function
    !Use the same as that for the distance calculation
    !Assign arrays for the kernel function
    rmax=MAXVAL(cosm%r)
    CALL fill_array(rmin,rmax,proj%r_X,proj%nX)
    WRITE(*,*) 'FILL_KERNEL: minimum r [Mpc/h]:', REAL(rmin)
    WRITE(*,*) 'FILL_KERNEL: maximum r [Mpc/h]:', REAL(rmax)
    WRITE(*,*) 'FILL_KERNEL: number of points:', nx

    IF(ALLOCATED(proj%X)) DEALLOCATE(proj%X)
    ALLOCATE(proj%X(nX))

    !Now fill the kernels
    DO i=1,nX
       r=proj%r_x(i)
       IF(ix==2) THEN
          proj%x(i)=y_kernel(r,cosm)
       ELSE IF(ix==10) THEN
          proj%x(i)=gravity_kernel(r,cosm)
       ELSE
          STOP 'FILL_KERNEL: Error, ix not specified correctly'
       END IF
    END DO

    WRITE(*,*) 'FILL_KERNEL: Done:'
    WRITE(*,*)

  END SUBROUTINE fill_kernel

  FUNCTION y_kernel(r,cosm)

    !The Compton-y projeciton kernel
    IMPLICIT NONE
    REAL :: y_kernel
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: z, a
    REAL :: crap

    z=redshift_r(r,cosm)
    a=scale_factor_z(z)

    !To stop compile-time warnings
    crap=cosm%om_m
    crap=r 

    y_kernel=yfac*mpc !Convert some units; note that there is no factor of 'a'
    y_kernel=y_kernel*a
    y_kernel=y_kernel*eV*cm**(-3) !Convert from eV cm^-3 to J m^-3

  END FUNCTION y_kernel

  FUNCTION gravity_kernel(r,cosm)

    !Projection kernel for gravitational waves (~ 1/r)
    IMPLICIT NONE
    REAL :: gravity_kernel
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap

    REAL, PARAMETER :: A=1.
    REAL, PARAMETER :: rmin=10.

    !To stop compile-time warnings
    crap=cosm%om_m

    IF(r<rmin) THEN
       gravity_kernel=A/rmin
    ELSE
       gravity_kernel=A/r
    END IF

  END FUNCTION gravity_kernel

  SUBROUTINE get_nz(ix,lens)

    !The the n(z) function for lensing
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens

    !CHARACTER(len=256) :: names(7)
    !names(1)='1 - RCSLenS'
    !names(2)='2 - KiDS (z = 0.1 -> 0.9)'
    !names(3)='3 - KiDS (z = 0.1 -> 0.3)'
    !names(4)='4 - KiDS (z = 0.3 -> 0.5)'
    !names(5)='5 - KiDS (z = 0.5 -> 0.7)'
    !names(6)='6 - KiDS (z = 0.7 -> 0.9)'
    !names(7)='7 - CFHTLenS (Van Waerbeke 2013)'

    !IF(inz==-1) THEN
    !   WRITE(*,*) 'GET_NZ: Choose n(z)'
    !   WRITE(*,*) '==================='
    !   DO i=1,SIZE(names)
    !      WRITE(*,*) TRIM(names(i))
    !   END DO
    !   READ(*,*) inz
    !   WRITE(*,*) '==================='
    !   WRITE(*,*)
    !END IF

    IF(ix==1 .OR. ix==4) THEN
       CALL fill_analytic_nz_table(ix,lens)
    ELSE
       CALL fill_nz_table(ix,lens)
    END IF

    !WRITE(*,*) 'GET_NZ: ', TRIM(names(inz))
    WRITE(*,*) 'GET_NZ: zmin:', lens%z_nz(1)
    WRITE(*,*) 'GET_NZ: zmax:', lens%z_nz(lens%nnz)
    WRITE(*,*) 'GET_NZ: nz:', lens%nnz
    WRITE(*,*)

  END SUBROUTINE get_nz

  SUBROUTINE fill_analytic_nz_table(ix,lens)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens
    INTEGER :: i

    REAL, PARAMETER :: zmin=0.
    REAL, PARAMETER :: zmax=2.5
    INTEGER, PARAMETER :: n=128

    !From analytical function
    lens%nnz=n
    IF(ALLOCATED(lens%z_nz)) DEALLOCATE(lens%z_nz)
    IF(ALLOCATED(lens%nz)) DEALLOCATE(lens%nz)
    ALLOCATE(lens%z_nz(lens%nnz),lens%nz(lens%nnz))

    !Fill the look-up tables
    CALL fill_array(zmin,zmax,lens%z_nz,lens%nnz)
    DO i=1,n
       lens%nz(i)=nz_lensing(lens%z_nz(i),ix)
    END DO

  END SUBROUTINE fill_analytic_nz_table

  SUBROUTINE fill_nz_table(ix,lens)

    USE file_info
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix
    TYPE(lensing), INTENT(INOUT) :: lens
    INTEGER :: i
    CHARACTER(len=256) :: input

    !Get file name
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
    ELSE
       STOP 'GET_NZ: ix not specified correctly'
    END IF
    WRITE(*,*) 'GET_NZ: Input file:', TRIM(input)

    !Allocate arrays
    lens%nnz=count_number_of_lines(input)
    IF(ALLOCATED(lens%z_nz)) DEALLOCATE(lens%z_nz)
    IF(ALLOCATED(lens%nz))   DEALLOCATE(lens%nz)
    ALLOCATE(lens%z_nz(lens%nnz),lens%nz(lens%nnz))

    !Read in n(z) table
    OPEN(7,file=input)
    DO i=1,lens%nnz
       READ(7,*) lens%z_nz(i), lens%nz(i)
    END DO
    CLOSE(7)

  END SUBROUTINE fill_nz_table

  FUNCTION nz_lensing(z,ix)

    IMPLICIT NONE
    REAL :: nz_lensing
    REAL, INTENT(IN) :: z
    INTEGER, INTENT(IN) :: ix
    REAL :: a, b, c, d, e, f, g, h, i
    REAL :: n1, n2, n3
    REAL :: z1, z2

    IF(ix==1) THEN
       !RCSLenS
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
       !CFHTLenS
       z1=0.7 !Not a free parameter in Van Waerbeke 2013
       z2=1.2 !Not a free parameter in Van Waerbeke 2013
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

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
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

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate_q=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=q_integrand(a,r,lens,cosm)
             f2=q_integrand(b,r,lens,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*DBLE(i-1)/DBLE(n-1)
                x=progression(a,b,i,n)
                fx=q_integrand(x,r,lens,cosm)
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
                STOP 'INTEGRATE_Q: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             integrate_q=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             integrate_q=0.d0
             STOP 'INTEGRATE_Q: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             integrate_q=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_q

  FUNCTION q_integrand(z,r,lens,cosm)

    !The lensing efficiency integrand, which is a function of z
    !z is integrated over while r is just a parameter
    !This is only called for n(z)
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
       !Find the r'(z) variable that is integrated over
       !rdash=find(z,cosm%z_r,cosm%r,cosm%nr,3,3,2)       
       rdash=comoving_distance(a,cosm)
       !Find the n(z)
       nz=find(z,lens%z_nz,lens%nz,lens%nnz,3,3,2)
       !This is then the integrand
       q_integrand=nz*f_k(rdash-r,cosm)/f_k(rdash,cosm)
    END IF

  END FUNCTION q_integrand

  FUNCTION integrate_Limber(l,a,b,logktab,logatab,logptab,nk,na,acc,iorder,proj,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
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

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=25

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate_Limber=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          !WRITE(*,*) j, n

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=Limber_integrand(a,l,logktab,logatab,logptab,nk,na,proj,cosm)
             f2=Limber_integrand(b,l,logktab,logatab,logptab,nk,na,proj,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                !x=a+(b-a)*DBLE(i-1)/DBLE(n-1)
                x=progression(a,b,i,n)
                fx=Limber_integrand(x,l,logktab,logatab,logptab,nk,na,proj,cosm)
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
                STOP 'INTEGRATE_LIMBER: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             integrate_Limber=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             integrate_Limber=0.d0
             STOP 'INTEGRATE_LIMBER: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             integrate_Limber=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION integrate_Limber

  FUNCTION Limber_integrand(r,l,logktab,logatab,logptab,nk,na,proj,cosm)

    !The integrand for the Limber integral
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

       !Get the two kernels
       DO i=1,2
          X(i)=find(r,proj(i)%r_X,proj(i)%X,proj(i)%nX,3,3,2)
       END DO
       !x2=find(r,proj%r_x2,proj%x2,proj%nx2,3,3,2)
       !x1=exp(find(log(r),log(proj%r_x1),log(proj%x1),proj%nx1,3,3,2)) !Barfed with this
       !x2=exp(find(log(r),log(proj%r_x2),log(proj%x2),proj%nx2,3,3,2)) !Barfed with this

       !Get variables r, z(r) and k(r) for P(k,z)
       z=redshift_r(r,cosm)
       a=scale_factor_z(z)
       k=(l+lcorr)/f_k(r,cosm) !LoVerde et al. (2008) Limber correction

       !Construct the integrand
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

    INTEGER, PARAMETER  :: n=1024 !Number of samples to take in r

    !Calculate the integral for this value of ell
    total=integrate_Limber(l,r1,r2,logktab,logatab,logptab,nk,na,acc_Limber,3,proj,cosm)

    !Now split up the contributions
    !You need the Jacobian and to remember that the contribution is split in ln(k), ln(z) and ln(R)
    !This means factors of:
    !k - f_k(r)/f'_k(r) (=r if flat)
    !z - z/H(z)
    !r - r
    OPEN(7,file=TRIM(outfile))
    DO i=1,n
       r=progression(r1,r2,i,n)
       IF(r==0.) THEN
          CYCLE
       ELSE
          k=(l+lcorr)/f_k(r,cosm)
          !z=find(r,cosm%r,cosm%z_r,cosm%nr,3,3,2)
          z=redshift_r(r,cosm)
          a=scale_factor_z(z)
          int=Limber_integrand(r,l,logktab,logatab,logptab,nk,na,proj,cosm)
          WRITE(7,*) k, int*f_k(r,cosm)/(fdash_k(r,cosm)*total), z, int*z/(sqrt(Hubble2(a,cosm))*total), r, int*r/total
       END IF
    END DO
    CLOSE(7)

  END SUBROUTINE Limber_contribution

  FUNCTION find_pka(k,a,logktab,logatab,logptab,nk,na)

    !Looks up the power as a 2D function of k and a
    IMPLICIT NONE
    REAL :: find_pka
    INTEGER, INTENT(IN) :: nk, na
    REAL, INTENT(IN) :: k, a
    REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk,na)

    REAL, PARAMETER :: kmin=1e-3 !kmin value
    REAL, PARAMETER :: kmax=1e2 !kmax value

    IF(k<=kmin .OR. k>=kmax) THEN
       find_pka=0.
    ELSE
       find_pka=exp(find2d(log(k),logktab,log(a),logatab,logptab,nk,na,3,3,1))
    END IF

    !Convert from Delta^2 -> P(k) - with dimensions of (Mpc/h)^3
    !find_pkz=(2.*pi**2)*find_pkz/k**3

  END FUNCTION find_pka
  
END MODULE Limber
