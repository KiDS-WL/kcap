MODULE logical_operations

CONTAINS

  FUNCTION positive(x)

    !Logical function that returns .TRUE. if x>=0.
    IMPLICIT NONE
    LOGICAL :: positive
    REAL, INTENT(IN) :: x

    IF(x<0.) THEN
       positive=.FALSE.
    ELSE
       positive=.TRUE.
    END IF
    
  END FUNCTION positive

  FUNCTION odd(i)

    IMPLICIT NONE
    LOGICAL :: odd
    INTEGER :: i

    !Tests for i being odd or even returns true if odd

    IF(mod(i,2)==0) THEN
       odd=.FALSE.
    ELSE
       odd=.TRUE.
    END IF

  END FUNCTION odd

END MODULE logical_operations
