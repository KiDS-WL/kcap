MODULE fix_polynomial
  
CONTAINS

  FUNCTION Lagrange_polynomial(x,n,xv,yv)

    !Computes the result of the nth order Lagrange polynomial at point x, L(x)
    IMPLICIT NONE
    REAL :: Lagrange_polynomial
    REAL, INTENT(IN) :: x, xv(n+1), yv(n+1)
    REAL :: l(n+1)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j

    !Initialise variables, one for sum and one for multiplication
    Lagrange_polynomial=0.
    l=1.

    !Loops to find the polynomials, one is a sum and one is a multiple
    DO i=0,n
       DO j=0,n
          IF(i .NE. j) l(i+1)=l(i+1)*(x-xv(j+1))/(xv(i+1)-xv(j+1))
       END DO
       Lagrange_polynomial=Lagrange_polynomial+l(i+1)*yv(i+1)
    END DO
    
  END FUNCTION Lagrange_polynomial
  
  SUBROUTINE fix_line(a1,a0,x1,y1,x2,y2)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2

    !Given xi, yi i=1,2 fits a line between these points

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

  END SUBROUTINE fix_line

  SUBROUTINE fix_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3

    !Given xi, yi i=1,2,3 fits a quadratic between these points

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.)-a1*x1

  END SUBROUTINE fix_quadratic

  SUBROUTINE fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3

    !Given xi, yi i=1,2,3,4 fits a cubic between these points

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2.+x4*x1+x1**2.)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3.-b*x1**2.-c*x1

  END SUBROUTINE fix_cubic

END MODULE fix_polynomial
