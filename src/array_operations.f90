MODULE array_operations

CONTAINS

  INTEGER FUNCTION array_position(x,a,n)

    ! Returns the location in the array of value x
    ! If x is not in array then returns zero
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
    INTEGER, INTENT(IN) :: a(n) ! Array to check 
    INTEGER, INTENT(IN) :: n    ! Size of array
    INTEGER :: i

    array_position=0
    
    DO i=1,n
       IF(a(i)==x) THEN
          array_position=i
          EXIT
       END IF
    END DO
    
  END FUNCTION array_position

  FUNCTION sum_double(a,n)

    ! Sum using double precision, which is necessary for many array elements
    IMPLICIT NONE
    REAL :: sum_double
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION :: sum
    INTEGER :: i
    
    sum=0.d0

    DO i=1,n
       sum=sum+a(i)
    END DO

    sum_double=REAL(sum)

  END FUNCTION sum_double

  SUBROUTINE amputate(arr,n_old,n_new)

    ! Chop an array down to a smaller size
    ! TODO: Retire
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    REAL, ALLOCATABLE :: hold(:)
    INTEGER, INTENT(IN) :: n_new
    INTEGER, INTENT(IN) :: n_old
    INTEGER :: i    

    IF(n_old<n_new) STOP 'AMPUTATE: Error, new array should be smaller than the old one'

    ALLOCATE(hold(n_old))
    hold=arr
    DEALLOCATE(arr)
    ALLOCATE(arr(n_new))
    
    DO i=1,n_new
       arr(i)=hold(i)
    END DO
    
    DEALLOCATE(hold)

  END SUBROUTINE amputate

  SUBROUTINE amputate_general(a,n,m,i1,i2)

    ! Chop an array of size a(n) down to a smaller size demarked by indices i1, i2
    ! If i1=1 and i2=n then this does nothing
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: m
    INTEGER, INTENT(IN) :: i1, i2
    REAL, ALLOCATABLE :: b(:)
    INTEGER :: i

    IF(i2<i1) THEN
       STOP 'AMPUTATE: Error, i2 should be greater than i1'
    END IF

    m=i2-i1+1
    IF(n<m) THEN
       STOP 'AMPUTATE: Error, new array should be smaller than the old one'
    END IF

    ! Store input array and then deallocate
    ALLOCATE(b(n))
    b=a
    DEALLOCATE(a)

    ! Allocate new output array
    ALLOCATE(a(m))

    ! Fill new output array
    DO i=1,m
       a(i)=b(i+i1-1)
    END DO

    ! Deallocate holding array
    DEALLOCATE(b)

  END SUBROUTINE amputate_general

  SUBROUTINE reduce(arr1,n1,arr2,n2)

    ! Reduces the size of array1 to the size of array2
    ! This will not preserve the spacing of entries in array1 and might be a terrible idea in many cases
    IMPLICIT NONE
    REAL, INTENT(IN) :: arr1(n1)
    REAL, INTENT(OUT) :: arr2(n2)
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER :: i, j

    DO i=1,n2
       j=1+CEILING(REAL((n1-1)*(i-1))/REAL(n2-1))
       arr2(i)=arr1(j)
    END DO

  END SUBROUTINE reduce

  SUBROUTINE reduceto(arr1,n)

    ! Reduces the array from whatever size to size 'n'
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: arr1(:)
    INTEGER, INTENT(IN) :: n
    REAL, ALLOCATABLE :: hold(:)
    INTEGER :: i, j

    ALLOCATE(hold(n))

    DO i=1,n
       j=1+CEILING(REAL((n-1)*(i-1))/REAL(n-1))
       hold(i)=arr1(j)
    END DO

    DEALLOCATE(arr1)
    ALLOCATE(arr1(n))

    arr1=hold

    DEALLOCATE(hold)

  END SUBROUTINE reduceto

  SUBROUTINE reverse(arry,n)

    ! Reverses the contents of arry
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: arry(n)
    INTEGER :: i
    REAL :: hold(n) 

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

  END SUBROUTINE reverse

  FUNCTION splay(a,n1,n2,n3)

    ! This splays out a 3d array 'a' into a 1d array 'b' of the same size (n1*n2*n3)
    IMPLICIT NONE
    REAL :: splay(n1*n2*n3)
    REAL, INTENT(IN) :: a(n1,n2,n3)
    INTEGER, INTENT(IN) :: n1, n2, n3
    INTEGER :: i, j, k, ii

    ! Set sum integer to zero
    ii=0

    DO i=1,n1
       DO j=1,n2
          DO k=1,n3             
             ii=ii+1
             splay(ii)=a(i,j,k)
          END DO
       END DO
    END DO

  END FUNCTION splay

  SUBROUTINE binning(a,a1,a2,n,b,c,ilog)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, ilog
    REAL, INTENT(IN) :: a(:)
    REAL, INTENT(OUT) :: b(n), c(n)
    REAL :: a1, a2, min, max
    REAL, ALLOCATABLE :: binlim(:)
    INTEGER :: i, j

    min=a1
    max=a2

    WRITE(*,*) 'Binning'
    WRITE(*,*) 'Min:', min
    WRITE(*,*) 'Max:', max

    IF(ilog==1) THEN
       min=log10(min)
       max=log10(max)
    END IF

    ! This sets the limits for the bins!
    CALL fill_array(min,max,binlim,n+1)

    ! This sets the centre value for each bin!
    DO i=1,n
       b(i)=(binlim(i)+binlim(i+1))/2.
    END DO

    IF(ilog==1) THEN
       binlim=10.**binlim
       b=10.**b
    END IF

    c=0.

    DO i=1,SIZE(a)
       DO j=1,n
          IF(a(i)>=binlim(j) .AND. a(i)<=binlim(j+1)) THEN
             c(j)=c(j)+1.
          END IF
       END DO
    END DO

    WRITE(*,*) 'Binning complete'
    WRITE(*,*)

  END SUBROUTINE binning

  SUBROUTINE merge_arrays(a,na,b,nb,c,nc)

    ! Takes arrays a and b and merges them together to make c with length SIZE(a)+SIZE(b)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a(na), b(nb)
    REAL, ALLOCATABLE, INTENT(OUT) :: c(:)
    INTEGER, INTENT(IN) :: na, nb
    INTEGER, INTENT(OUT) :: nc
    INTEGER :: i
    
    nc=na+nb

    IF(ALLOCATED(c)) DEALLOCATE(c)
    ALLOCATE(c(nc))

    DO i=1,na
       c(i)=a(i)
    END DO

    DO i=1,nb
       c(i+na)=b(i)
    END DO

  END SUBROUTINE merge_arrays

  FUNCTION concatenate_arrays(a,na,b,nb)

    ! Concatenate arrays a and b to form new array with length SIZE(a)+SIZE(b)
    IMPLICIT NONE
    REAL :: concatenate_arrays(na+nb)
    REAL, INTENT(IN) :: a(na), b(nb)
    INTEGER, INTENT(IN) :: na, nb
    INTEGER :: i

    DO i=1,na
       concatenate_arrays(i)=a(i)
    END DO

    DO i=1,nb
       concatenate_arrays(i+na)=b(i)
    END DO

  END FUNCTION concatenate_arrays

  SUBROUTINE fill_array(min,max,arr,n)

    ! Fills array 'arr' in equally spaced intervals
    ! TODO: I'm not sure if inputting an array like this is okay
    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL, ALLOCATABLE, INTENT(OUT) :: arr(:)
    INTEGER, INTENT(IN) :: n

    ! Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=progression(min,max,i,n)
       END DO
    END IF

  END SUBROUTINE fill_array

  SUBROUTINE fill_array_double(min,max,arr,n)

    ! Fills array 'arr' in equally spaced intervals
    ! TODO: I'm not sure if inputting an array like this is okay
    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    INTEGER, INTENT(IN) :: n

    ! Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.d0

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=progression8(min,max,i,n)
       END DO
    END IF

  END SUBROUTINE fill_array_double

  FUNCTION progression(xmin,xmax,i,n)

    IMPLICIT NONE
    REAL :: progression
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    IF(n==1) THEN
       progression=xmin
    ELSE
       progression=xmin+(xmax-xmin)*REAL(i-1)/REAL(n-1)
    END IF
    
  END FUNCTION progression

  FUNCTION progression_log(xmin,xmax,i,n)

    IMPLICIT NONE
    REAL :: progression_log
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression_log=exp(progression(log(xmin),log(xmax),i,n))
    
  END FUNCTION progression_log

  FUNCTION progression8(xmin,xmax,i,n)

    IMPLICIT NONE
    DOUBLE PRECISION :: progression8
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression8=xmin+(xmax-xmin)*DBLE(i-1)/DBLE(n-1)
    
  END FUNCTION progression8

  FUNCTION maximum(x,y,n)

    USE fix_polynomial
    
    ! From an array y(x) finds the x location of the first maximum
    IMPLICIT NONE
    REAL :: maximum
    REAL, INTENT(IN) :: x(n), y(n)
    INTEGER, INTENT(IN) :: n
    REAL :: x1, x2, x3, y1, y2, y3, a, b, c
    INTEGER :: i

    ! Need this to stop a compile-time warning
    maximum=0.

    DO i=1,n-1
       
       IF(y(i+1)<y(i)) THEN

          ! Get the x positions
          x1=x(i-1)
          x2=x(i)
          x3=x(i+1)

          ! Get the y values
          y1=y(i-1)
          y2=y(i)
          y3=y(i+1)

          ! Fix a quadratic around the maximum
          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          ! Read off the maximum x from the parabola
          maximum=-b/(2.*a)

          ! Exit the loop
          EXIT
          
       ELSE IF(i<n-1) THEN
          
          CYCLE
          
       ELSE
          
          STOP 'MAXIMUM: Error, array does not have a maximum'
          
       END IF
       
    END DO

  END FUNCTION maximum

  SUBROUTINE mask(okay,m,n,min,max)

    ! Flags objects that make the cut as 'okay'
    ! Can be applied to any scalar array, not just mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m(n), min, max
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(OUT) :: okay(n)
    INTEGER :: i, o

    WRITE(*,*) 'MASK: Imposing property cut'    
    WRITE(*,*) 'MASK: Minimum value:', min
    WRITE(*,*) 'MASK: Maximum value:', max
    WRITE(*,*) 'MASK: Original number of objects', n

    okay=.FALSE.

    DO i=1,n
       IF(m(i)>=min .AND. m(i)<=max) okay(i)=.TRUE.
    END DO

    o=COUNT(okay)

    WRITE(*,*) 'MASK: Final number of objects:', o
    WRITE(*,*) 'MASK: Fraction remaining:', REAL(o)/REAL(n)
    WRITE(*,*) 'MASK: Fraction culled:', 1.-REAL(o)/REAL(n)
    WRITE(*,*) 'MASK: Done'
    WRITE(*,*)

  END SUBROUTINE mask

END MODULE array_operations
