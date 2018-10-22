module HMx_setup
  use cosmology_functions, only: cosmology
  use HMx, only: halomod
  implicit none

  type :: field
    character(len=256) :: name
    integer :: id
  end type field

  type :: HMx_setup_config    
    real(8) :: kmin, kmax
    integer :: nk

    real(8) :: zmin, zmax, amin, amax, mmin, mmax
    integer :: nz

    character(len=256) :: p_lin_source, hm_mode, matter_matter_section_name

    integer :: ihm, iw, icosmo
    logical :: response, dimensionless_power_spectrum
    
    type(field), dimension(:), allocatable :: fields
    integer, dimension(:), allocatable :: itype
    integer :: nt

    type(cosmology) :: cosm
    type(halomod) :: hm

    real(8), dimension(:), allocatable :: k, a
    logical :: verbose
  end type HMx_setup_config

  contains
    function read_fields(string) result(fields)
      character(len=512), intent(in) :: string
      type(field), dimension(:), allocatable:: fields

      integer, parameter :: n_field_max = 10
      type(field), dimension(n_field_max) :: fields_tmp
      character(len=512) :: rest, field_str
      integer :: i, id, idx

      rest = trim(string)
      i = 1
      do while (trim(rest) /= "" .and. i < n_field_max)
        idx = scan(trim(rest), " ", back=.true.)
        if(idx <= 1) then
          field_str = rest
          rest = ""
        else
          field_str = rest(idx+1:)
          rest = rest(:idx)
        end if

        if(trim(field_str) == "matter") then
          fields_tmp(i)%name = "matter"
          fields_tmp(i)%id = 0
        else if(trim(field_str) == "dm") then
          fields_tmp(i)%name = "dm"
          fields_tmp(i)%id = 1
        else if(trim(field_str) == "gas") then
          fields_tmp(i)%name = "gas"
          fields_tmp(i)%id = 2
        else if(trim(field_str) == "stars") then
          fields_tmp(i)%name = "stars"
          fields_tmp(i)%id = 3
        else if(trim(field_str) == "pressure") then
          fields_tmp(i)%name = "pressure"
          fields_tmp(i)%id = 6
        else
          write(*,*) "Field ", trim(field_str), " is not supported."
        end if

        i = i + 1
      end do

      allocate(fields, source=fields_tmp(:i-1))
      fields = [(fields(i), i=size(fields),1,-1)]
    end function read_fields

end module HMx_setup

function setup(options) result(result)
  use HMx_setup
  use HMx
  use array_operations
  use cosmosis_modules
  implicit none

  ! Arguments
  integer(cosmosis_block), value :: options
  ! Return value
  type(c_ptr) :: result
  ! Variables
  integer(cosmosis_status) :: status
  type(HMx_setup_config), pointer :: HMx_config
  integer :: i, verbose, response, icosmo, dimensionless_power_spectrum
  character(len=512) :: field_list_string

  allocate(HMx_config)

  if(datablock_get(options, option_section, "nz", HMx_config%nz) /= 0) then
     write(*,*) "Could not load nz."
     stop
  end if
  if(datablock_get(options, option_section, "zmin", HMx_config%zmin) /= 0) then
     write(*,*) "Could not load zmin."
     stop
  end if
  if(datablock_get(options, option_section, "zmax", HMx_config%zmax) /= 0) then
     write(*,*) "Could not load zmax."
     stop
  end if
  HMx_config%amin = 1.0/(1+HMx_config%zmax)
  HMx_config%amax = 1.0/(1+HMx_config%zmin)

  if(datablock_get(options, option_section, "nk", HMx_config%nk) /= 0) then
     write(*,*) "Could not load nk."
     stop
  end if
  if(datablock_get(options, option_section, "kmin", HMx_config%kmin) /= 0) then
     write(*,*) "Could not load kmin."
     stop
  end if
  if(datablock_get(options, option_section, "kmax", HMx_config%kmax) /= 0) then
     write(*,*) "Could not load kmax."
     stop
  end if

  status = datablock_get_string(options, option_section, "fields", field_list_string)
  HMx_config%fields = read_fields(field_list_string)
  HMx_config%nt = size(HMx_config%fields)
  allocate(HMx_config%itype(HMx_config%nt))
  HMx_config%itype = [(HMx_config%fields(i)%id, i=1,HMx_config%nt)]

  write(*,*) "Read ", HMx_config%nt, " fields."
  do i=1,HMx_config%nt
    write(*,*) HMx_config%fields(i)%name
  end do
  status = datablock_get_string_default(options, option_section, "matter_matter_section_name", "", HMx_config%matter_matter_section_name)

  status = datablock_get_double_default(options, option_section, "mmin", 1e7, HMx_config%mmin)
  status = datablock_get_double_default(options, option_section, "mmax", 1e17, HMx_config%mmax)

  status = datablock_get_int_default(options, option_section, "verbose", 1, verbose)

  ! Get halo model mode
  status = datablock_get_string_default(options, option_section, "hm_mode", "hmx", HMx_config%hm_mode)
  if(trim(HMx_config%hm_mode) == "hmx") then
    HMx_config%ihm = 6
  else if(trim(HMx_config%hm_mode) == "hmcode") then
    HMx_config%ihm = 1
  else if(trim(HMx_config%hm_mode) == "vanilla_halo_model") then
    HMx_config%ihm = 3
  else
    write(*,*) "hm_mode = ", trim(HMx_config%hm_mode), " not supported."
    stop
  end if
  ! Get ihm value directly if supplied
  status = datablock_get_int_default(options, option_section, "ihm", HMx_config%ihm, HMx_config%ihm)

  HMx_config%verbose = verbose > 0

  HMx_config%icosmo = 1
  
  status = datablock_get_string_default(options, option_section, "p_lin_source", "eh", HMx_config%p_lin_source)
  if(trim(HMx_config%p_lin_source) == "external") then
    ! Use linear power spectrum provided by CosmoSIS
    HMx_config%cosm%itk = 4
  else if(trim(HMx_config%p_lin_source) == "eh") then
    ! Use Eisenstein & Hu transfer function
    HMx_config%cosm%itk = 1
  else
    write(*,*) "p_lin_source = ", trim(HMx_config%p_lin_source), " not supported."
    stop
  end if

  status = datablock_get_int_default(options, option_section, "de_model", 1, HMx_config%cosm%iw)

  status = datablock_get_int_default(options, option_section, "response", 0, response)
  HMx_config%response = response == 1

  !Set power spectrum units
  status = datablock_get_int_default(options, option_section, "dimensionless_power_spectrum", 1, dimensionless_power_spectrum)
  HMx_config%dimensionless_power_spectrum = dimensionless_power_spectrum == 1

  ! Create k array (log spacing)
  call fill_array(log(HMx_config%kmin), log(HMx_config%kmax), HMx_config%k, HMx_config%nk)
  HMx_config%k = exp(HMx_config%k)
  
  ! Create a arrays. The projection module wants increasing z, so a is decreasing
  call fill_array(HMx_config%amax, HMx_config%amin, HMx_config%a, HMx_config%nz)

  result = c_loc(HMx_config)

end function setup

function execute(block, config) result(status)
  use cosmosis_modules
  use HMx_setup
  use HMx, only : calculate_HMx, assign_halomod
  use cosmology_functions, only : init_cosmology, print_cosmology, assign_cosmology
  use constants

  implicit none
  !Arguments
  integer(cosmosis_block), value :: block
  type(c_ptr), value :: config
  !Return value
  integer(cosmosis_status) :: status
  !Variables
  type(HMx_setup_config), pointer :: HMx_config
  integer :: i, j
  real(8) :: log10_eps, log10_M0, log10_whim, log10_Theat
  real(8), dimension(:), allocatable :: k_plin, z_plin
  real(8), dimension(:,:), allocatable :: pk_lin
  real(8), dimension(:,:,:,:), allocatable :: pk_1h, pk_2h, pk_full
  character(len=256) :: pk_section

  call c_f_pointer(config, HMx_config)

  if(allocated(k_plin)) deallocate(k_plin)
  if(allocated(z_plin)) deallocate(z_plin)
  if(allocated(pk_lin)) deallocate(pk_lin)
  if(allocated(pk_1h)) deallocate(pk_1h)
  if(allocated(pk_2h)) deallocate(pk_2h)
  if(allocated(pk_full)) deallocate(pk_full)

  if(allocated(HMx_config%cosm%log_k_plin)) deallocate(HMx_config%cosm%log_k_plin)
  if(allocated(HMx_config%cosm%log_plin)) deallocate(HMx_config%cosm%log_plin)

! Assign default values.
  call assign_cosmology(HMx_config%icosmo, HMx_config%cosm, HMx_config%verbose)
  call assign_halomod(HMx_config%ihm, HMx_config%hm, HMx_config%verbose)

  ! Cosmology parameters
  status = datablock_get_double(block, cosmological_parameters_section, "omega_m", HMx_config%cosm%om_m)
  status = datablock_get_double(block, cosmological_parameters_section, "omega_lambda", HMx_config%cosm%om_v)
  status = datablock_get_double(block, cosmological_parameters_section, "omega_b", HMx_config%cosm%om_b)
  status = datablock_get_double(block, cosmological_parameters_section, "omega_nu", HMx_config%cosm%om_nu)
  status = datablock_get_double(block, cosmological_parameters_section, "h0", HMx_config%cosm%h)
  status = datablock_get_double(block, cosmological_parameters_section, "sigma_8", HMx_config%cosm%sig8)
  status = datablock_get_double(block, cosmological_parameters_section, "n_s", HMx_config%cosm%n)
  status = datablock_get_double(block, cosmological_parameters_section, "w", HMx_config%cosm%w)
  status = datablock_get_double(block, cosmological_parameters_section, "wa", HMx_config%cosm%wa)
  status = datablock_get_double(block, cosmological_parameters_section, "neff", HMx_config%cosm%neff)
  status = datablock_get_double(block, cosmological_parameters_section, "T_cmb", HMx_config%cosm%T_cmb)
  status = datablock_get_double(block, cosmological_parameters_section, "z_cmb", HMx_config%cosm%z_cmb)
  status = datablock_get_double(block, cosmological_parameters_section, "Y_H", HMx_config%cosm%YH)
  status = datablock_get_double(block, cosmological_parameters_section, "omega_w", HMx_config%cosm%om_w)
  status = datablock_get_double(block, cosmological_parameters_section, "inv_m_wdm", HMx_config%cosm%inv_m_wdm)

  ! HMCode parameters
  status = datablock_get_double(block, halo_model_parameters_section, "Dv0", HMx_config%hm%Dv0)
  status = datablock_get_double(block, halo_model_parameters_section, "Dv1", HMx_config%hm%Dv1)
  status = datablock_get_double(block, halo_model_parameters_section, "dc0", HMx_config%hm%dc0)
  status = datablock_get_double(block, halo_model_parameters_section, "dc1", HMx_config%hm%dc1)
  status = datablock_get_double(block, halo_model_parameters_section, "eta0", HMx_config%hm%eta0)
  status = datablock_get_double(block, halo_model_parameters_section, "eta1", HMx_config%hm%eta1)
  status = datablock_get_double(block, halo_model_parameters_section, "f0", HMx_config%hm%f0)
  status = datablock_get_double(block, halo_model_parameters_section, "f1", HMx_config%hm%f1)
  status = datablock_get_double(block, halo_model_parameters_section, "ks", HMx_config%hm%ks)
  status = datablock_get_double(block, halo_model_parameters_section, "As", HMx_config%hm%As)
  status = datablock_get_double(block, halo_model_parameters_section, "alp0", HMx_config%hm%alp0)
  status = datablock_get_double(block, halo_model_parameters_section, "alp1", HMx_config%hm%alp1)

  ! HOD parameters
  status = datablock_get_double(block, halo_model_parameters_section, "m_gal", HMx_config%hm%mgal)
  status = datablock_get_double(block, halo_model_parameters_section, "HImin", HMx_config%hm%HImin)
  status = datablock_get_double(block, halo_model_parameters_section, "HImax", HMx_config%hm%HImax)
  
  ! HMx baryon parameters
  status = datablock_get_double(block, halo_model_parameters_section, "alpha", HMx_config%hm%alpha)
  status = datablock_get_double_default(block, halo_model_parameters_section, "log10_eps", log10(HMx_config%hm%eps), log10_eps)
  status = datablock_get_double(block, halo_model_parameters_section, "Gamma", HMx_config%hm%Gamma)
  status = datablock_get_double_default(block, halo_model_parameters_section, "log10_M0", log10(HMx_config%hm%M0), log10_M0)
  status = datablock_get_double(block, halo_model_parameters_section, "Astar", HMx_config%hm%Astar)
  status = datablock_get_double_default(block, halo_model_parameters_section, "log10_whim", log10(HMx_config%hm%Twhim), log10_whim)
  status = datablock_get_double(block, halo_model_parameters_section, "cstar", HMx_config%hm%cstar)
  status = datablock_get_double(block, halo_model_parameters_section, "sstar", HMx_config%hm%sstar)
  status = datablock_get_double(block, halo_model_parameters_section, "mstar", HMx_config%hm%mstar)
  status = datablock_get_double_default(block, halo_model_parameters_section, "log10_Theat", log10(HMx_config%hm%Theat), log10_Theat)
  status = datablock_get_double(block, halo_model_parameters_section, "fcold", HMx_config%hm%fcold)

  ! Exponentiate those parameters that will be explored in log space
  HMx_config%hm%eps = 10**log10_eps
  HMx_config%hm%M0 = 10**log10_M0
  HMx_config%hm%Twhim = 10**log10_whim
  HMx_config%hm%Theat = 10**log10_Theat

  if(trim(HMx_config%p_lin_source) == "external") then
     status = datablock_get_double_grid(block, matter_power_lin_section, &
          "k_h", k_plin, &
          "z", z_plin, &
          "p_k", pk_lin)
     if(status /= 0) then
        write(*,*) "Could not load load linear power spectrum."
        stop
     end if
     HMx_config%cosm%has_power = .true.
     HMx_config%cosm%n_plin = size(k_plin)
     if(.not. allocated(HMx_config%cosm%log_k_plin)) allocate(HMx_config%cosm%log_k_plin(size(k_plin)))
     if(.not. allocated(HMx_config%cosm%log_plin)) allocate(HMx_config%cosm%log_plin(size(k_plin)))
     HMx_config%cosm%log_k_plin = log(k_plin)
     HMx_config%cosm%log_plin = log(pk_lin(:,1)*k_plin**3/(2*pi**2))
  end if

  call init_cosmology(HMx_config%cosm)
  if(HMx_config%verbose) then
     call print_cosmology(HMx_config%cosm)
     WRITE(*,*) "HALOMOD parameters:"
     WRITE(*,*) "alpha    :", HMx_config%hm%alpha
     WRITE(*,*) "eps      :", HMx_config%hm%eps
     WRITE(*,*) "Gamma    :", HMx_config%hm%Gamma
     WRITE(*,*) "M0       :", HMx_config%hm%M0
     WRITE(*,*) "Astar    :", HMx_config%hm%Astar
     WRITE(*,*) "Twhim    :", HMx_config%hm%Twhim
     WRITE(*,*) "cstar    :", HMx_config%hm%cstar
     WRITE(*,*) "sstar    :", HMx_config%hm%sstar
     WRITE(*,*) "mstar    :", HMx_config%hm%mstar
     WRITE(*,*) "Theat    :", HMx_config%hm%Theat
     WRITE(*,*) "fcold    :", HMx_config%hm%fcold
  end if

  call calculate_HMx(HMx_config%itype, HMx_config%nt, &
                     HMx_config%mmin, HMx_config%mmax, &
                     HMx_config%k, HMx_config%nk, &
                     HMx_config%a, HMx_config%nz, &
                     pk_lin, pk_2h, pk_1h, pk_full, &
                     HMx_config%hm, HMx_config%cosm, &
                     HMx_config%verbose, HMx_config%response)

  if(.not. HMx_config%dimensionless_power_spectrum) then
    ! Remove the k^3/2pi^2 factor
    forall (i=1:HMx_config%nk) pk_full(:,:,i,:) = pk_full(:,:,i,:)*2*pi**2/HMx_config%k(i)**3
  endif

  do i=1,HMx_config%nt
    do j=1,i
      pk_section = trim(HMx_config%fields(i)%name)//"_"//trim(HMx_config%fields(j)%name)//"_power_spectrum"
      if(i==1 .and. j==1 .and. trim(HMx_config%matter_matter_section_name) /= "") then
        pk_section = HMx_config%matter_matter_section_name
      end if

      status = datablock_put_double_grid(block, pk_section, &
                                         "k_h", HMx_config%k, &
                                         "z", 1.0/HMx_config%a-1.0, &
                                         "p_k", pk_full(i,j,:,:))
    end do
  end do

  if(status /= 0) then
     write(*,*) "Failed to write NL power spectrum to datablock."
  end if

end function execute

function cleanup(config) result(status)
  use cosmosis_modules
  use HMx_setup

  ! Arguments
  type(c_ptr), value :: config
  ! Return value
  integer(cosmosis_status) :: status
  ! Variables
  type(HMx_setup_config), pointer :: HMx_config  

  ! Free memory allocated in the setup function
  call c_f_pointer(config, HMx_config)
  !deallocate(HMx_config)

  status = 0
end function cleanup
