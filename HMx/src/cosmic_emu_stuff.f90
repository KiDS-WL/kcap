MODULE cosmic_emu_stuff

  USE cosmology_functions

CONTAINS

  SUBROUTINE get_FrankenEmu_power(k,P,n,z,cosm,rebin)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), P(:)
    INTEGER, INTENT(OUT) :: n
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    LOGICAL, INTENT(IN) :: rebin
    CHARACTER(len=256) :: output
    INTEGER :: i, nk
    REAL, ALLOCATABLE :: k2(:), P2(:)
    REAL :: kmin, kmax
    
    CALL SYSTEM('rm emu_params.txt')
    CALL SYSTEM('rm emu_power.dat')

    OPEN(7,file='emu_params.txt')
    WRITE(7,fmt='(A20,7F10.5)') 'emu_power.dat', cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%n, -cosm%w, cosm%sig8, cosm%h, z
    CLOSE(7)

    CALL SYSTEM('/Users/Mead/Physics/FrankenEmu_Cl/emu.exe emu_params.txt')! > /dev/null')

    output='emu_power.dat'

    n=file_length(output,verbose=.FALSE.)
    n=n-5
    IF(ALLOCATED(k)) DEALLOCATE(k)
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE(k(n),P(n))

    WRITE(*,*) 'GET_FRANKENEMU_POWER: z:', z
    WRITE(*,*) 'GET_FRANKENEMU_POWER: P(k) file length:', n
    WRITE(*,*) 'GET_FRANKENEMU_POWER: Reading in P(k)'

    OPEN(7,file=output)
    READ(7,*)
    READ(7,*)
    READ(7,*)
    READ(7,*)
    READ(7,*)
    DO i=1,n
       READ(7,*) k(i), P(i)
    END DO
    CLOSE(7)

    !Convert P(k) to Delta^2(k)
    !P=P*(k**3)*4.*pi/(2.*pi)**3

    !Convert k to k/h
    k=k/cosm%h

    ! Rebin on a log-linear axis
    IF(rebin) THEN
       kmin=1e-2
       kmax=7.
       nk=128
       CALL fill_array(log(kmin),log(kmax),k2,nk)
       k2=exp(k2)
       ALLOCATE(P2(nk))
       CALL interpolate_array(log(k),log(P),n,log(k2),P2,nk,3,3,2)
       P2=exp(P2)
       DEALLOCATE(k,P)
       ALLOCATE(k(nk),P(nk))
       k=k2
       P=P2
       n=nk
       DEALLOCATE(k2,P2)
    END IF

    WRITE(*,*) 'GET_FRANKENEMU_POWER: Done'
    WRITE(*,*)

  END SUBROUTINE get_FrankenEmu_power

  SUBROUTINE get_Mira_Titan_power(k,P,n,z,cosm,rebin)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), P(:)
    INTEGER, INTENT(OUT) :: n
    REAL, INTENT(IN) :: z    
    TYPE(cosmology), INTENT(IN) :: cosm
    LOGICAL, INTENT(IN) :: rebin
    CHARACTER(len=256) :: output
    INTEGER :: i, nk
    REAL, ALLOCATABLE :: k2(:), P2(:)
    REAL :: kmin, kmax

    CALL SYSTEM('rm xstar.dat')
    CALL SYSTEM('rm EMU0.txt')

    OPEN(7,file='xstar.dat')
    WRITE(7,*) (cosm%Om_m*cosm%h**2), (cosm%Om_b*cosm%h**2), cosm%sig8, cosm%h, cosm%n, cosm%w, cosm%wa, (cosm%om_nu*cosm%h**2), z
    CLOSE(7)

    CALL SYSTEM('/Users/Mead/Physics/MiraTitan/P_tot/emu.exe')
    output='EMU0.txt'

    n=file_length(output,verbose=.FALSE.)
    IF(ALLOCATED(k)) DEALLOCATE(k)
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE(k(n),P(n))

    WRITE(*,*) 'GET_MIRA_TITAN_POWER: z:', z
    WRITE(*,*) 'GET_MIRA_TITAN_POWER: P(k) file length:', n
    WRITE(*,*) 'GET_MIRA_TITAN_POWER: Reading in P(k)'

    OPEN(7,file=output)
    DO i=1,n
       READ(7,*) k(i), P(i)
    END DO
    CLOSE(7)

    !Convert P(k) to Delta^2(k)
    P=P*(k**3)*4.*pi/(2.*pi)**3

    !Convert k to k/h
    k=k/cosm%h

    ! Rebin on a log-linear axis
    IF(rebin) THEN
       kmin=1e-2
       kmax=7.
       nk=128
       CALL fill_array(log(kmin),log(kmax),k2,nk)
       k2=exp(k2)
       ALLOCATE(P2(nk))
       CALL interpolate_array(log(k),log(P),n,log(k2),P2,nk,3,3,2)
       P2=exp(P2)
       DEALLOCATE(k,P)
       ALLOCATE(k(nk),P(nk))
       k=k2
       P=P2
       n=nk
       DEALLOCATE(k2,P2)
    END IF

    WRITE(*,*) 'GET_MIRA_TITAN_POWER: Done'
    WRITE(*,*)

  END SUBROUTINE get_Mira_Titan_power
  
END MODULE cosmic_emu_stuff
