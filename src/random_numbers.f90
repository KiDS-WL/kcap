MODULE random_numbers

  !TODO: Think about using intrinsic 'random_number' instead of 'rand()'
  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE RNG_set(seed)

    ! Seeds the RNG
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: seed
    INTEGER :: int, timearray(3)
    REAL*4 :: rand ! Necessary to define for ifort, also the *4 is necessary

    WRITE(*,*) 'RNG_SET: Initialising random number generator'
    WRITE(*,*) 'RNG_SET: Seed:', seed

    IF(seed==0) THEN
       
       ! This fills the time array using the system clock!
       ! If called within the same second the numbers will be identical!
       CALL itime(timeArray)
       
       ! This then initialises the generator!
       int=FLOOR(rand(timeArray(1)+timeArray(2)+timeArray(3)))
       
    ELSE
       
       ! In this case you can keep track of the seed
       int=FLOOR(rand(seed))
       
    END IF
    
    WRITE(*,*) 'RNG_SET: Done'
    WRITE(*,*)

  END SUBROUTINE RNG_set

  FUNCTION random_integer(i1,i2)

    ! Picks an integer with uniform random probability between i1 and i2 spaced with 1
    IMPLICIT NONE
    INTEGER :: random_integer
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER :: n
    REAL*4 :: rand ! Necessary to define for ifort

    ! Range for the number
    n=1+i2-i1

    random_integer=i1-1+CEILING(rand(0)*REAL(n))

  END FUNCTION random_integer

  FUNCTION random_uniform(x1,x2)

    ! Produces a uniform random number between x1 and x2
    IMPLICIT NONE
    REAL :: random_uniform
    REAL, INTENT(IN) :: x1,x2
    REAL*4 :: rand !I think this needs to be defined for ifort

    ! rand is some inbuilt function
    random_uniform=x1+(x2-x1)*(rand(0))

  END FUNCTION random_uniform

  FUNCTION random_Rayleigh(sigma)

    ! Produces a Rayleigh-distributed random number
    IMPLICIT NONE
    REAL :: random_Rayleigh
    REAL, INTENT(IN) :: sigma
    REAL, PARAMETER :: small=1e-10

    ! Problems if small=0. because log(0.) gets called sometimes
    random_Rayleigh=sigma*sqrt(-2.*log(random_uniform(small,1.)))

  END FUNCTION random_Rayleigh

  FUNCTION random_Lorentzian()

    ! Produces a Lorentzian-distributed random number
    USE constants
    IMPLICIT NONE
    REAL :: random_Lorentzian

    random_Lorentzian=tan(random_uniform(0.,pi/2.))

  END FUNCTION random_Lorentzian

  FUNCTION random_Gaussian_pair(mean,sigma)
    
    ! Gets a pair of Gaussian random numbers
    USE constants
    IMPLICIT NONE
    REAL :: random_Gaussian_pair(2)
    REAL, INTENT(IN) :: mean, sigma
    REAL :: r, theta

    r=random_Rayleigh(sigma)
    theta=random_uniform(0.,twopi)

    ! Both of these numbers are Gaussian
    random_Gaussian_pair(1)=r*sin(theta)+mean
    random_Gaussian_pair(2)=r*cos(theta)+mean

  END FUNCTION random_Gaussian_pair

  FUNCTION random_Gaussian(mean,sigma)
    
    ! Gets a single Gaussian random number
    IMPLICIT NONE
    REAL :: random_Gaussian
    REAL, INTENT(IN) :: mean, sigma
    REAL :: G(2)

    ! This is wasteful as there is a second, independent Gaussian random number
    G=random_Gaussian_pair(mean,sigma)
    random_Gaussian=G(1)

  END FUNCTION random_Gaussian

  FUNCTION random_lognormal(mean_x,sigma_lnx)
    
    ! Gets a single Gaussian random number
    ! mean_x: <x>
    ! sigma_lnx: rms of the logarithm of x
    IMPLICIT NONE
    REAL :: random_lognormal
    REAL, INTENT(IN) :: mean_x, sigma_lnx
    REAL :: mu, sigma, G(2)

    sigma=sigma_lnx
    mu=log(mean_x)-0.5*sigma**2

    ! This is wasteful as there is a second, independent Gaussian random number
    G=random_Gaussian_pair(mu,sigma)
    random_lognormal=exp(G(1))

  END FUNCTION random_lognormal

  REAL FUNCTION random_exponential(mean)

    ! Produces a exponentially-distributed random number
    IMPLICIT NONE
    REAL, INTENT(IN) :: mean

    ! small is introducted because there will be problems here if log(0) is ever called
    REAL, PARAMETER :: small=1e-10
  
    random_exponential=-mean*log(random_uniform(small,1.))

  END FUNCTION random_exponential

  REAL FUNCTION random_polynomial(n)

    ! Generate a polynomailly distributed number [x:0->1]
    IMPLICIT NONE
    REAL, INTENT(IN) :: n

    random_polynomial=(random_uniform(0.,1.))**(1./(n+1))

  END FUNCTION random_polynomial

  REAL FUNCTION random_theta()

    ! A random spherical angle such that the space is equally populated
    IMPLICIT NONE

    random_theta=acos(random_uniform(-1.,1.))

  END FUNCTION random_theta

END MODULE random_numbers
