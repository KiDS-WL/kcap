MODULE cosmic_emu_stuff

  USE cosmology_functions

CONTAINS

  SUBROUTINE get_cosmic_emu_power(z,output,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    CHARACTER(len=256) :: output

    STOP 'GET_COSMIC_EMU_POWER: Not finished yet'

    !OPEN(7,file=TRIM(stem)//'FrankenEmu/random.txt')
    WRITE(7,fmt='(A20,7F10.5)') output, (cosm%om_m*cosm%h**2.), (cosm%om_b*cosm%h**2.), cosm%ns, -cosm%w, cosm%sig8, cosm%h, z
    CLOSE(7)

    CALL SYSTEM('/Users/Mead/Physics/FrankenEmu_CL/emu.exe ./../FrankenEmu/random.txt > /dev/null')

  END SUBROUTINE get_cosmic_emu_power

  SUBROUTINE get_Mira_Titan_power(k,P,n,z,cosm,rebin)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), P(:)
    REAL, INTENT(IN) :: z
    INTEGER, INTENT(OUT) :: n
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

    n=file_length(output,.FALSE.)
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
