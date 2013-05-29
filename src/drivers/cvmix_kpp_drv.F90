!BOP
!\newpage
! !ROUTINE: cvmix_kpp_driver

! !DESCRIPTION: A routine to test the KPP module.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_kpp_driver(nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_data_type,          &
                                    cvmix_kpp_params_type
  use cvmix_kpp,             only : cvmix_init_kpp,           &
                                    cvmix_coeffs_kpp  
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)       :: CVmix_vars
  type(cvmix_kpp_params_type) :: CVmix_kpp_params

  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity
  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  integer :: fid

  allocate(diffusivity(nlev,2))
  allocate(viscosity(nlev))

  call cvmix_put(CVmix_vars, 'nlev', nlev)
  CVmix_vars%diff_iface => diffusivity(:,:)
  CVmix_vars%visc_iface => viscosity(:)

  call cvmix_init_kpp(CVmix_kpp_params)
  call cvmix_coeffs_kpp(CVMix_vars, CVmix_kpp_params)

#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"visc", "diff"/))

  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_kpp_driver
