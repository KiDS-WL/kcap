MODULE interpolate

  USE fix_polynomial
  USE table_integer
  USE array_operations

CONTAINS

  FUNCTION find(x,xin,yin,n,iorder,ifind,imeth)

    !Given two arrays x and y this routine interpolates to find the y_i value at position x_i
    IMPLICIT NONE
    REAL :: find
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xin(n), yin(n)
    REAL ::  xtab(n), ytab(n)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i
    INTEGER, INTENT(IN) :: iorder, ifind, imeth

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    !ifind = 1 => find x in xtab quickly assuming the table is linearly spaced
    !ifind = 2 => find x in xtab by crudely searching from x(1) to x(n)
    !ifind = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !imeth = 1 => Uses standard polynomials for interpolation
    !imeth = 2 => Uses Lagrange polynomials for interpolation

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab,n)
       CALL reverse(ytab,n)
    END IF

    IF(x<xtab(1)) THEN

       !Do a linear interpolation beyond the table boundary

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       IF(imeth==1) THEN
          CALL fix_line(a,b,x1,y1,x2,y2)
          find=a*x+b
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

    ELSE IF(x>xtab(n)) THEN

       !Do a linear interpolation beyond the table boundary

       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       IF(imeth==1) THEN
          CALL fix_line(a,b,x1,y1,x2,y2)
          find=a*x+b
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

    ELSE IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x1=xtab(1)
          x2=xtab(2)

          y1=ytab(1)
          y2=ytab(2)

       ELSE IF (x>=xtab(n-1)) THEN

          x1=xtab(n-1)
          x2=xtab(n)

          y1=ytab(n-1)
          y2=ytab(n)

       ELSE

          i=select_table_integer(x,xtab,n,ifind)

          x1=xtab(i)
          x2=xtab(i+1)

          y1=ytab(i)
          y2=ytab(i+1)

       END IF

       IF(imeth==1) THEN
          CALL fix_line(a,b,x1,y1,x2,y2)
          find=a*x+b
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x1=xtab(1)
             x2=xtab(2)
             x3=xtab(3)

             y1=ytab(1)
             y2=ytab(2)
             y3=ytab(3)

          ELSE IF (x>=xtab(n-1)) THEN

             x1=xtab(n-2)
             x2=xtab(n-1)
             x3=xtab(n)

             y1=ytab(n-2)
             y2=ytab(n-1)
             y3=ytab(n)

          END IF

          IF(imeth==1) THEN
             CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
             find=a*(x**2)+b*x+c
          ELSE IF(imeth==2) THEN
             find=Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))
          ELSE
             STOP 'FIND: Error, method not specified correctly'
          END IF

       ELSE

          i=select_table_integer(x,xtab,n,ifind)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          IF(imeth==1) THEN
             !In this case take the average of two separate quadratic spline values
             CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
             find=(a*x**2+b*x+c)/2.
             CALL fix_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
             find=find+(a*x**2+b*x+c)/2.
          ELSE IF(imeth==2) THEN
             !In this case take the average of two quadratic Lagrange polynomials
             find=(Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))+Lagrange_polynomial(x,2,(/x2,x3,x4/),(/y2,y3,y4/)))/2.
          ELSE
             STOP 'FIND: Error, method not specified correctly'
          END IF

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x1=xtab(1)
          x2=xtab(2)
          x3=xtab(3)
          x4=xtab(4)        

          y1=ytab(1)
          y2=ytab(2)
          y3=ytab(3)
          y4=ytab(4)

       ELSE IF (x>=xtab(n-2)) THEN

          x1=xtab(n-3)
          x2=xtab(n-2)
          x3=xtab(n-1)
          x4=xtab(n)

          y1=ytab(n-3)
          y2=ytab(n-2)
          y3=ytab(n-1)
          y4=ytab(n)

       ELSE

          i=select_table_integer(x,xtab,n,ifind)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       IF(imeth==1) THEN
          CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
          find=a*x**3+b*x**2+c*x+d
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

    ELSE

       STOP 'FIND: Error, interpolation order specified incorrectly'

    END IF

  END FUNCTION find

  FUNCTION find2d(x,xin,y,yin,fin,nx,ny,iorder,ifind,imeth)

    !A 2D interpolation routine to find value f(x,y) at position x, y
    IMPLICIT NONE
    REAL :: find2d
    INTEGER, INTENT(IN) :: nx, ny
    REAL, INTENT(IN) :: x, xin(nx), y, yin(ny), fin(nx,ny)
    REAL ::  xtab(nx), ytab(ny), ftab(nx,ny)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    REAL :: f11, f12, f13, f14
    REAL :: f21, f22, f23, f24
    REAL :: f31, f32, f33, f34
    REAL :: f41, f42, f43, f44
    REAL :: f10, f20, f30, f40
    REAL :: f01, f02, f03, f04
    INTEGER :: i1, i2, i3, i4
    INTEGER :: j1, j2, j3, j4
    REAL :: findx, findy
    INTEGER :: i, j
    INTEGER, INTENT(IN) :: iorder, ifind, imeth

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    !ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
    !ifind = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !imeth = 1 => Uses cubic polynomials for interpolation
    !imeth = 2 => Uses Lagrange polynomials for interpolation

    IF(imeth==2) STOP 'No Lagrange polynomials for you'

    xtab=xin
    ytab=yin
    ftab=fin

    IF(xtab(1)>xtab(nx)) STOP 'FIND2D: x table in wrong order'
    IF(ytab(1)>ytab(ny)) STOP 'FIND2D: y table in wrong order'

    IF((x<xtab(1) .OR. x>xtab(nx)) .AND. (y>ytab(ny) .OR. y<ytab(1))) THEN
       WRITE(*,*) 'FIND2D: point xmin:', xtab(1)
       WRITE(*,*) 'FIND2D: point xmax:', xtab(nx)
       WRITE(*,*) 'FIND2D: point x:', x
       WRITE(*,*) 'FIND2D: point ymin:', ytab(1)
       WRITE(*,*) 'FIND2D: point ymax:', ytab(ny)
       WRITE(*,*) 'FIND2D: point y:', y
       STOP 'FIND2D: Desired point is outside x AND y table range'
    END IF

    IF(iorder==1) THEN

       IF(nx<2) STOP 'FIND2D: Not enough x points in your table for linear interpolation'
       IF(ny<2) STOP 'FIND2D: Not enough y points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          i=1

       ELSE IF (x>=xtab(nx-1)) THEN

          i=nx-1

       ELSE

          i=select_table_integer(x,xtab,nx,ifind)

       END IF

       i1=i
       i2=i+1

       x1=xtab(i1)
       x2=xtab(i2)      

       IF(y<=ytab(2)) THEN

          j=1

       ELSE IF (y>=ytab(ny-1)) THEN

          j=ny-1

       ELSE

          j=select_table_integer(y,ytab,ny,ifind)

       END IF

       j1=j
       j2=j+1

       y1=ytab(j1)
       y2=ytab(j2)

       !

       f11=ftab(i1,j1)
       f12=ftab(i1,j2)

       f21=ftab(i2,j1)
       f22=ftab(i2,j2)

       !y direction interpolation

       CALL fix_line(a,b,x1,f11,x2,f21)
       f01=a*x+b

       CALL fix_line(a,b,x1,f12,x2,f22)
       f02=a*x+b

       CALL fix_line(a,b,y1,f01,y2,f02)
       findy=a*y+b

       !x direction interpolation

       CALL fix_line(a,b,y1,f11,y2,f12)
       f10=a*y+b

       CALL fix_line(a,b,y1,f21,y2,f22)
       f20=a*y+b

       CALL fix_line(a,b,x1,f10,x2,f20)
       findx=a*x+b

       !

       !Final result is an average over each direction
       find2d=(findx+findy)/2.

    ELSE IF(iorder==2) THEN

       STOP 'FIND2D: Quadratic 2D interpolation not implemented - also probably pointless'

    ELSE IF(iorder==3) THEN

       IF(x<xtab(1) .OR. x>xtab(nx)) THEN

          IF(nx<2) STOP 'FIND2D: Not enough x points in your table for linear interpolation'
          IF(ny<4) STOP 'FIND2D: Not enough y points in your table for cubic interpolation'

          !x is off the table edge

          IF(x<xtab(1)) THEN

             i1=1
             i2=2

          ELSE

             i1=nx-1
             i2=nx

          END IF

          x1=xtab(i1)
          x2=xtab(i2)

          IF(y<=ytab(4)) THEN

             j=2

          ELSE IF (y>=ytab(ny-3)) THEN

             j=ny-2

          ELSE

             j=select_table_integer(y,ytab,ny,ifind)

          END IF

          j1=j-1
          j2=j
          j3=j+1
          j4=j+2

          y1=ytab(j1)
          y2=ytab(j2)
          y3=ytab(j3)
          y4=ytab(j4)

          f11=ftab(i1,j1)
          f12=ftab(i1,j2)
          f13=ftab(i1,j3)
          f14=ftab(i1,j4)

          f21=ftab(i2,j1)
          f22=ftab(i2,j2)
          f23=ftab(i2,j3)
          f24=ftab(i2,j4)

          !y interpolation
          CALL fix_cubic(a,b,c,d,y1,f11,y2,f12,y3,f13,y4,f14)
          f10=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f21,y2,f22,y3,f23,y4,f24)
          f20=a*y**3+b*y**2+c*y+d

          !x interpolation
          CALL fix_line(a,b,x1,f10,x2,f20)
          find2d=a*x+b

       ELSE IF(y<ytab(1) .OR. y>ytab(ny)) THEN

          !y is off the table edge

          IF(nx<4) STOP 'FIND2D: Not enough x points in your table for cubic interpolation'
          IF(ny<2) STOP 'FIND2D: Not enough y points in your table for linear interpolation'

          IF(x<=xtab(4)) THEN

             i=2

          ELSE IF (x>=xtab(nx-3)) THEN

             i=nx-2

          ELSE

             i=select_table_integer(x,xtab,nx,ifind)

          END IF

          i1=i-1
          i2=i
          i3=i+1
          i4=i+2

          x1=xtab(i1)
          x2=xtab(i2)
          x3=xtab(i3)
          x4=xtab(i4)

          IF(y<ytab(1)) THEN

             j1=1
             j2=2

          ELSE

             j1=ny-1
             j2=ny

          END IF

          y1=ytab(j1)
          y2=ytab(j2)

          f11=ftab(i1,j1)
          f21=ftab(i2,j1)
          f31=ftab(i3,j1)
          f41=ftab(i4,j1)

          f12=ftab(i1,j2)
          f22=ftab(i2,j2)
          f32=ftab(i3,j2)
          f42=ftab(i4,j2)

          !x interpolation

          CALL fix_cubic(a,b,c,d,x1,f11,x2,f21,x3,f31,x4,f41)
          f01=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f12,x2,f22,x3,f32,x4,f42)
          f02=a*x**3+b*x**2+c*x+d

          !y interpolation

          CALL fix_line(a,b,y1,f01,y2,f02)
          find2d=a*y+b

       ELSE

          !Points exists within table boundardies (normal)

          IF(nx<4) STOP 'FIND2D: Not enough x points in your table for cubic interpolation'
          IF(ny<4) STOP 'FIND2D: Not enough y points in your table for cubic interpolation'

          IF(x<=xtab(4)) THEN

             i=2

          ELSE IF (x>=xtab(nx-3)) THEN

             i=nx-2

          ELSE

             i=select_table_integer(x,xtab,nx,ifind)

          END IF

          i1=i-1
          i2=i
          i3=i+1
          i4=i+2

          x1=xtab(i1)
          x2=xtab(i2)
          x3=xtab(i3)
          x4=xtab(i4)

          IF(y<=ytab(4)) THEN

             j=2

          ELSE IF (y>=ytab(ny-3)) THEN

             j=ny-2

          ELSE

             j=select_table_integer(y,ytab,ny,ifind)

          END IF

          j1=j-1
          j2=j
          j3=j+1
          j4=j+2

          y1=ytab(j1)
          y2=ytab(j2)
          y3=ytab(j3)
          y4=ytab(j4)

          !

          f11=ftab(i1,j1)
          f12=ftab(i1,j2)
          f13=ftab(i1,j3)
          f14=ftab(i1,j4)

          f21=ftab(i2,j1)
          f22=ftab(i2,j2)
          f23=ftab(i2,j3)
          f24=ftab(i2,j4)

          f31=ftab(i3,j1)
          f32=ftab(i3,j2)
          f33=ftab(i3,j3)
          f34=ftab(i3,j4)

          f41=ftab(i4,j1)
          f42=ftab(i4,j2)
          f43=ftab(i4,j3)
          f44=ftab(i4,j4)

          !x interpolation

          CALL fix_cubic(a,b,c,d,x1,f11,x2,f21,x3,f31,x4,f41)
          f01=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f12,x2,f22,x3,f32,x4,f42)
          f02=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f13,x2,f23,x3,f33,x4,f43)
          f03=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f14,x2,f24,x3,f34,x4,f44)
          f04=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,y1,f01,y2,f02,y3,f03,y4,f04)
          findy=a*y**3+b*y**2+c*y+d

          !y interpolation

          CALL fix_cubic(a,b,c,d,y1,f11,y2,f12,y3,f13,y4,f14)
          f10=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f21,y2,f22,y3,f23,y4,f24)
          f20=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f31,y2,f32,y3,f33,y4,f34)
          f30=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f41,y2,f42,y3,f43,y4,f44)
          f40=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,x1,f10,x2,f20,x3,f30,x4,f40)
          findx=a*x**3+b*x**2+c*x+d

          !Final result is an average over each direction
          find2d=(findx+findy)/2.

       END IF

    ELSE

       STOP 'FIND2D: order for interpolation not specified correctly'

    END IF

  END FUNCTION find2d

  SUBROUTINE interpolate_array(x1,y1,n1,x2,y2,n2,iorder,ifind,imeth)

    !Interpolates array 'x1-y1' onto new 'x' values x2 and output y2
    IMPLICIT NONE
    REAL, INTENT(IN) :: x1(n1), y1(n1), x2(n2)
    REAL, INTENT(OUT) :: y2(n2)
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER, INTENT(IN) :: iorder, ifind, imeth
    INTEGER :: i

    !Could be more efficient, but probably not worth the hassle
    !It does 'find integer' every time

    DO i=1,n2
       y2(i)=find(x2(i),x1,y1,n1,iorder,ifind,imeth)
    END DO

  END SUBROUTINE interpolate_array

END MODULE interpolate
