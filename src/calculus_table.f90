MODULE calculus_table

CONTAINS

  FUNCTION derivative_table(x,xin,yin,n,iorder,imeth)

    USE table_integer
    USE fix_polynomial
    USE array_operations

    ! Given two arrays x and y such that y=y(x) this uses interpolation to calculate the derivative y'(x_i) at position x_i
    IMPLICIT NONE
    REAL :: derivative_table
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xin(n), yin(n)
    REAL ::  xtab(n), ytab(n)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i
    INTEGER, INTENT(IN) :: imeth, iorder

    ! This version interpolates if the value is off either end of the array!
    ! Care should be chosen to insert x, xtab, ytab as log if this might give better!
    ! Results from the interpolation!

    ! imeth = 1 => find x in xtab by crudely searching
    ! imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    ! imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    ! iorder = 1 => linear interpolation
    ! iorder = 2 => quadratic interpolation
    ! iorder = 3 => cubic interpolation

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       ! Reverse the arrays in this case
       CALL reverse(xtab,n)
       CALL reverse(ytab,n)
    END IF

    IF(iorder==1) THEN

       IF(n<2) STOP 'DERIVATIVE_TABLE: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x2=xtab(2)
          x1=xtab(1)

          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-1)) THEN

          x2=xtab(n)
          x1=xtab(n-1)

          y2=ytab(n)
          y1=ytab(n-1)

       ELSE

          i=select_table_integer(x,xtab,n,imeth)

          x2=xtab(i+1)
          x1=xtab(i)

          y2=ytab(i+1)
          y1=ytab(i)

       END IF

       CALL fix_line(a,b,x1,y1,x2,y2)
       derivative_table=a

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'DERIVATIVE_TABLE: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x3=xtab(3)
             x2=xtab(2)
             x1=xtab(1)

             y3=ytab(3)
             y2=ytab(2)
             y1=ytab(1)

          ELSE IF (x>=xtab(n-1)) THEN

             x3=xtab(n)
             x2=xtab(n-1)
             x1=xtab(n-2)

             y3=ytab(n)
             y2=ytab(n-1)
             y1=ytab(n-2)

          END IF

          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          derivative_table=2.*a*x+b

       ELSE

          i=select_table_integer(x,xtab,n,imeth)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          ! In this case take the average of two separate quadratic spline values

          derivative_table=0.

          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          derivative_table=derivative_table+(2.*a*x+b)/2.

          CALL fix_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          derivative_table=derivative_table+(2.*a*x+b)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'DERIVATIVE_TABLE: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x4=xtab(4)
          x3=xtab(3)
          x2=xtab(2)
          x1=xtab(1)

          y4=ytab(4)
          y3=ytab(3)
          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-2)) THEN

          x4=xtab(n)
          x3=xtab(n-1)
          x2=xtab(n-2)
          x1=xtab(n-3)

          y4=ytab(n)
          y3=ytab(n-1)
          y2=ytab(n-2)
          y1=ytab(n-3)

       ELSE

          i=select_table_integer(x,xtab,n,imeth)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       derivative_table=3.*a*(x**2.)+2.*b*x+c

    ELSE

       STOP 'DERIVATIVE_TABLE: Error, order not specified correctly'

    END IF

  END FUNCTION derivative_table

  FUNCTION integrate_table(x,y,n,n1,n2,iorder)

    USE fix_polynomial

    ! Integrates tables y(x)dx
    IMPLICIT NONE
    REAL :: integrate_table
    INTEGER, INTENT(IN) :: n, n1, n2
    REAL, INTENT(IN) :: x(n), y(n)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    DOUBLE PRECISION :: sum
    INTEGER :: i, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    sum=0.d0

    ! I think if n1=n2 then the result will just be zero anyway
    ! IF(n2<=n1) STOP 'INTEGRATE_TABLE: Error n2 must be greater than n1'

    IF(iorder==1) THEN

       ! Sums over all Trapezia (a+b)*h/2
       DO i=n1,n2-1
          a=y(i+1)
          b=y(i)
          h=x(i+1)-x(i)
          sum=sum+(a+b)*h/2.d0
       END DO

    ELSE IF(iorder==2) THEN

       DO i=n1,n2-2

          x1=x(i)
          x2=x(i+1)
          x3=x(i+2)

          y1=y(i)
          y2=y(i+1)
          y3=y(i+2)

          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
          q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
          q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

          ! Takes value for first and last sections but averages over sections where you
          ! have two independent estimates of the area
          IF(n==3) THEN
             sum=sum+q3-q1
          ELSE IF(i==1) THEN
             sum=sum+(q2-q1)+(q3-q2)/2.d0
          ELSE IF(i==n-2) THEN
             sum=sum+(q2-q1)/2.d0+(q3-q2)
          ELSE
             sum=sum+(q3-q1)/2.
          END IF

       END DO

    ELSE IF(iorder==3) THEN

       DO i=n1,n2-1

          ! First choose the integers used for defining cubics for each section
          ! First and last are different because the section does not lie in the *middle* of a cubic

          IF(i==1) THEN

             i1=1
             i2=2
             i3=3
             i4=4

          ELSE IF(i==n-1) THEN

             i1=n-3
             i2=n-2
             i3=n-1
             i4=n

          ELSE

             i1=i-1
             i2=i
             i3=i+1
             i4=i+2

          END IF

          x1=x(i1)
          x2=x(i2)
          x3=x(i3)
          x4=x(i4)

          y1=y(i1)
          y2=y(i2)
          y3=y(i3)
          y4=y(i4)

          CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

          ! These are the limits of the particular section of integral
          xi=x(i)
          xf=x(i+1)

          qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
          qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

          sum=sum+qf-qi

       END DO

    ELSE

       STOP 'INTEGRATE_TABLE: Error, order not specified correctly'

    END IF

    integrate_table=REAL(sum)

  END FUNCTION integrate_table

END MODULE calculus_table
