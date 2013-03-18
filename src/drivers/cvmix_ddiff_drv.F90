!BOP
!\newpage
! !ROUTINE: cvmix_ddiff_driver

! !DESCRIPTION: A routine to test the double diffusion mixing module.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_ddiff_driver(nlev)

! !USES:

  use cvmix_kinds_and_types, only : one,                      &
                                    cvmix_r8,                 &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_ddiff_params_type
  use cvmix_ddiff,           only : cvmix_init_ddiff,         &
                                    cvmix_coeffs_ddiff
  use cvmix_put_get,         only : cvmix_put
  use cvmix_output,          only : cvmix_output_open,        &
                                    cvmix_output_write,       &
                                    cvmix_output_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)          :: CVmix_vars
  type(cvmix_global_params_type) :: CVmix_params
  type(cvmix_ddiff_params_type)  :: CVmix_ddiff_params

  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity
  real(cvmix_r8), dimension(:),   allocatable, target :: Rrho_num, Rrho_denom

  ! column / file indices
  integer :: k, fid

  ! Namelist variables
  real(cvmix_r8) :: ddiff_exp1, strat_param_max, kappa_ddiff_t

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/ddiff_nml/ddiff_exp1, strat_param_max, kappa_ddiff_t

  ! Allocate memory to store diffusivity values
  ! Also store numerator / denominator for stratification parameter
  allocate(diffusivity(nlev+1,2))
  allocate(Rrho_num(nlev+1), Rrho_denom(nlev+1))
  do k=1,nlev+1
    Rrho_num(k) = dble(k-1)/dble(nlev)+one
    Rrho_denom(k) = one
  end do

  ! Initialization for CVMix data types
  call cvmix_put(CVmix_params,  'max_nlev', nlev)
  call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
  call cvmix_put(CVmix_vars,    'nlev',     nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars%diff_iface => diffusivity
  CVmix_vars%strat_param_num => Rrho_num
  CVmix_vars%strat_param_denom => Rrho_denom

  ! Read / set double diffusion parameters
  read(*, nml=ddiff_nml)
  call cvmix_init_ddiff(CVmix_ddiff_params,'mks', ddiff_exp1=ddiff_exp1, &
                        strat_param_max=strat_param_max,                 &
                        kappa_ddiff_t=kappa_ddiff_t)
  ! Debug option
  if (.false.) then
    print*, "Parameters are as follows:"
    print*, "strat_param_max = ", CVmix_ddiff_params%strat_param_max
    print*, "ddiff_exp1 = ", CVmix_ddiff_params%ddiff_exp1
    print*, "ddiff_exp2 = ", CVmix_ddiff_params%ddiff_exp2
    print*, "kappa_ddiff_param1 = ", CVmix_ddiff_params%kappa_ddiff_param1
    print*, "kappa_ddiff_param2 = ", CVmix_ddiff_params%kappa_ddiff_param2
    print*, "kappa_ddiff_param3 = ", CVmix_ddiff_params%kappa_ddiff_param3
    print*, "kappa_ddiff_t = ", CVmix_ddiff_params%kappa_ddiff_t
    print*, "kappa_ddiff_s = ", CVmix_ddiff_params%kappa_ddiff_s
    print*, "mol_diff = ", CVmix_ddiff_params%mol_diff
  end if

  call cvmix_coeffs_ddiff(CVmix_vars, CVmix_ddiff_params)
  ! For continuity of plot, set diffusivity when Rrho = 1
  diffusivity(1,1) = CVmix_ddiff_params%kappa_ddiff_t

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call cvmix_output_open(fid, "data.nc", "nc")
#else
  call cvmix_output_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"Rrho", "diff"/))

  call cvmix_output_close(fid)

!EOC

End Subroutine cvmix_ddiff_driver
