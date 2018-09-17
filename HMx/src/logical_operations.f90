MODULE logical_operations

CONTAINS

  LOGICAL FUNCTION positive(x)

    ! Logical function that returns .TRUE. if x>=0.
    ! If x=0. then it will return true
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    IF(x<0.) THEN
       positive=.FALSE.
    ELSE
       positive=.TRUE.
    END IF
    
  END FUNCTION positive

  LOGICAL FUNCTION even(i)

    ! Tests for i being even, returns true if even
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i

    IF(odd(i)) THEN
       even=.FALSE.
    ELSE
       even=.TRUE.
    END IF
    
  END FUNCTION even

  LOGICAL FUNCTION odd(i)

    ! Tests for i being odd, returns true if odd
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i  

    IF(mod(i,2)==0) THEN
       odd=.FALSE.
    ELSE
       odd=.TRUE.
    END IF

  END FUNCTION odd

  LOGICAL FUNCTION requal(x,y)

    ! Tests if two REAL numbers are within some tolerance of each other
    ! If they are, then they should be considered equal
    ! Adapted from https://stackoverflow.com/questions/4915462/how-should-i-do-floating-point-comparison/4915891#4915891
    IMPLICIT NONE
    REAL, INTENT(IN) :: x, y    
    REAL :: absx, absy, diff
    REAL, PARAMETER :: eps=1e-4

    STOP 'REQUAL: Test this'

    absx=ABS(x)
    absy=ABS(y)
    diff=ABS(x-y)

    IF(x==y) THEN
       requal=.TRUE.
    ELSE IF(x==0. .OR. y==0. .OR. diff<TINY(x)) THEN
       IF(diff<eps*TINY(x)) THEN
          requal=.TRUE.
       ELSE
          requal=.FALSE.
       END IF
    ELSE
       IF(diff/(absx+absy)<eps) THEN
          requal=.TRUE.
       ELSE
          requal=.FALSE.
       END IF
    END IF
    
  END FUNCTION requal

END MODULE logical_operations
