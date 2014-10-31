!BOP
!\newpage
! !ROUTINE: cvmix_ddiff_driver

! !DESCRIPTION: A routine to test the double diffusion mixing module.
!\\
!\\

! !INTERFACE:

Subroutine cvmix_ddiff_driver(nlev, max_nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                &
                                    cvmix_one,               &
                                    cvmix_data_type
  use cvmix_ddiff,           only : cvmix_init_ddiff,        &
                                    cvmix_coeffs_ddiff,      &
                                    cvmix_get_ddiff_real
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,           &
                                    cvmix_output_write,      &
#ifdef _NETCDF
                                    cvmix_output_write_att,  &
#endif
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev,               &! number of levels for column
                         max_nlev             ! number of columns in memory

!EOP
!BOC

  integer, parameter :: ncol = 2

  ! CVMix datatypes
  type(cvmix_data_type), dimension(ncol) :: CVmix_vars

  real(cvmix_r8), dimension(:,:), allocatable, target :: Tdiff, Sdiff
  real(cvmix_r8), dimension(:),   allocatable, target :: Rrho_num, Rrho_denom

  ! column / file indices
  integer :: k, fid, ic

  ! Namelist variables
  real(cvmix_r8) :: ddiff_exp1, strat_param_max

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/ddiff_nml/ddiff_exp1, strat_param_max

  print*, "Active levels: ", nlev
  print*, "Levels allocated in memory: ", max_nlev

  ! Allocate memory to store diffusivity values
  ! Also store numerator / denominator for stratification parameter
  allocate(Tdiff(max_nlev+1,ncol), Sdiff(max_nlev+1,ncol))
  allocate(Rrho_num(max_nlev), Rrho_denom(max_nlev))
  do k=1,nlev/2
    ! For first column, Rrho varies from 1 to 2
    Rrho_num(k) = real(k-1,cvmix_r8)/real(nlev/2-1,cvmix_r8)+cvmix_one
    Rrho_denom(k) = cvmix_one
    ! For second column, 1/Rrho varies from 1 to 10
    ! (Note: last column has diff=0, hence only using nlev instead of nlev+1)
    Rrho_num(k+nlev/2)   = -cvmix_one
    Rrho_denom(k+nlev/2) = -real(9*(k-1),cvmix_r8)/real(nlev/2-1,cvmix_r8) -  &
                            cvmix_one
  end do

  ! Point CVmix_vars values to memory allocated above
  do ic=1,ncol
    call cvmix_put(CVmix_vars(ic), 'nlev', nlev)
    call cvmix_put(CVmix_vars(ic), 'max_nlev', nlev)
    CVmix_vars(ic)%Tdiff_iface => Tdiff(:,ic)
    CVmix_vars(ic)%Sdiff_iface => Sdiff(:,ic)
    CVmix_vars(ic)%strat_param_num => Rrho_num(:)
    CVmix_vars(ic)%strat_param_denom => Rrho_denom(:)
  end do

  ! Read / set double diffusion parameters
  read(*, nml=ddiff_nml)
  call cvmix_init_ddiff(ddiff_exp1=ddiff_exp1, diff_conv_type="MC76",         &
                        strat_param_max=strat_param_max)
  call cvmix_coeffs_ddiff(CVmix_vars(1))

  call cvmix_init_ddiff(ddiff_exp1=ddiff_exp1, diff_conv_type="K88",          &
                        strat_param_max=strat_param_max)
  call cvmix_coeffs_ddiff(CVmix_vars(2))

  ! Output
#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"Rrho ", "Tdiff", "Sdiff"/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "long_name", "double diffusion " //        &
                              "stratification parameter", var_name="Rrho")
  call cvmix_output_write_att(fid, "units", "unitless", var_name="Rrho")
  call cvmix_output_write_att(fid, "long_name", "temperature diffusivity",    &
                              var_name="Tdiff")
  call cvmix_output_write_att(fid, "long_name", "salinity diffusivity",       &
                              var_name="Sdiff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="Tdiff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="Sdiff")
#endif
  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_ddiff_driver
