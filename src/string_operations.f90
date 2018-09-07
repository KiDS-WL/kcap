MODULE string_operations

CONTAINS

  !I cannot remember where (or why) I stole this from
  elemental subroutine str2int(str,int,stat)
    
    implicit none
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat

    read(str,*,iostat=stat) int
    
  end subroutine str2int

  FUNCTION number_file(fbase,i,fext)

    IMPLICIT NONE
    CHARACTER(len=256) :: number_file
    CHARACTER(len=*), INTENT(IN) :: fbase, fext
    INTEGER, INTENT(IN) :: i
    CHARACTER(len=8) num

    IF(i<0) STOP 'NUMBER_FILE: Error: cannot write negative number file names'

    IF(i<10) THEN       
       WRITE(num,fmt='(I1)') i
    ELSE IF(i<100) THEN
       WRITE(num,fmt='(I2)') i
    ELSE IF(i<1000) THEN
       WRITE(num,fmt='(I3)') i
    END IF

    number_file=TRIM(fbase)//TRIM(num)//TRIM(fext)

  END FUNCTION number_file

  FUNCTION number_file2(fbase,i1,mid,i2,fext)

    IMPLICIT NONE
    CHARACTER(len=256) ::number_file2
    CHARACTER(len=*), INTENT(IN) :: fbase, fext, mid
    INTEGER, INTENT(IN) :: i1, i2
    CHARACTER(len=8) :: num1, num2

    IF(i1<10) THEN       
       WRITE(num1,fmt='(I1)') i1
    ELSE IF(i1<100) THEN
       WRITE(num1,fmt='(I2)') i1
    ELSE IF(i1<1000) THEN
       WRITE(num1,fmt='(I3)') i1
    END IF

    IF(i2<10) THEN       
       WRITE(num2,fmt='(I1)') i2
    ELSE IF(i2<100) THEN
       WRITE(num2,fmt='(I2)') i2
    ELSE IF(i2<1000) THEN
       WRITE(num2,fmt='(I3)') i2
    END IF

    number_file2=TRIM(fbase)//TRIM(num1)//TRIM(mid)//TRIM(num2)//TRIM(fext)

  END FUNCTION number_file2

  FUNCTION number_file_zeroes(fbase,i,num,fext)

    !Number a file with zero padding
    !Num specifies the number of digits
    IMPLICIT NONE
    CHARACTER(len=256) :: number_file_zeroes
    CHARACTER(len=*), INTENT(IN) :: fbase, fext
    CHARACTER(len=4) :: num4
    CHARACTER(len=3) :: num3
    CHARACTER(len=2) :: num2
    CHARACTER(len=1) :: num1
    INTEGER, INTENT(IN) :: i
    INTEGER :: maxnum
    INTEGER, INTENT(IN) :: num

    maxnum=4

    IF(i<0) STOP 'NUMBER_FILE_ZEROES: Error: cannot write negative number file names'

    IF(num>maxnum) STOP 'NUMBER_FILE_ZEROES: Error: need to add extra number capacity'

    IF(num==1) THEN

       IF(i>=10) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
       WRITE(num1,fmt='(I0.1)') i
       number_file_zeroes=TRIM(fbase)//num1//TRIM(fext)

    ELSE IF(num==2) THEN

       IF(i>=100) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
       WRITE(num2,fmt='(I0.2)') i
       number_file_zeroes=TRIM(fbase)//num2//TRIM(fext)

    ELSE IF(num==3) THEN

       IF(i>=1000) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
       WRITE(num3,fmt='(I0.3)') i
       number_file_zeroes=TRIM(fbase)//num3//TRIM(fext)

    ELSE IF(num==4) THEN

       IF(i>=10000) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
       WRITE(num4,fmt='(I0.4)') i
       number_file_zeroes=TRIM(fbase)//TRIM(num4)//TRIM(fext)

    END IF

  END FUNCTION number_file_zeroes

END MODULE string_operations
