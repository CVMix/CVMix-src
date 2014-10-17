!BOP
!\newpage
! !ROUTINE: cvmix_shear_driver

! !DESCRIPTION: A routine to test the Large, et al., implementation of shear
!  mixing. Inputs are the coefficients used in Equation (28) of the paper.
!  The diffusivity coefficient is output from a single column to allow
!  recreation of the paper's Figure 3. Note that here each "level" of the
!  column denotes a different local gradient Richardson number rather than a
!  physical ocean level. All memory is declared in the driver, and the CVMix
!  data type points to the local variables.
!\\
!\\

! !INTERFACE:

Subroutine cvmix_shear_driver(nlev, max_nlev)

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_zero,               &
                                    cvmix_one,                &
                                    cvmix_data_type
  use cvmix_shear,           only : cvmix_init_shear,         &
                                    cvmix_coeffs_shear
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
#ifdef _NETCDF
                                    cvmix_output_write_att,   &
#endif
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev,               &! number of levels for column
                         max_nlev             ! number of columns in memory

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)          :: CVmix_vars

  real(cvmix_r8), dimension(:), allocatable, target :: Ri_g
  real(cvmix_r8), dimension(:), allocatable, target :: Mdiff, Tdiff

  integer :: kw, fid

  ! Namelist variables
  ! KPP mixing parameters for column
  real(cvmix_r8) :: KPP_nu_zero, KPP_Ri_zero, KPP_exp

  ! Namelist with shear mixing parameters
  namelist/KPP_nml/KPP_nu_zero, KPP_Ri_zero, KPP_exp

  print*, "Active levels: ", nlev
  print*, "Levels allocated in memory: ", max_nlev

  ! Ri_g should increase from 0 to 1 in active portion of level
  allocate(Ri_g(max_nlev+1))
  Ri_g(1) = cvmix_zero
  do kw = 2,max_nlev+1
    Ri_g(kw) = Ri_g(kw-1) + cvmix_one/real(nlev,cvmix_r8)
  end do

  ! Allocate memory to store viscosity and diffusivity values
  allocate(Mdiff(max_nlev+1), Tdiff(max_nlev+1))

  ! Initialization for CVMix data type
  call cvmix_put(CVmix_vars,    'nlev',     nlev)
  call cvmix_put(CVmix_vars,    'max_nlev', max_nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars%Mdiff_iface => Mdiff
  CVmix_vars%Tdiff_iface => Tdiff
  CVmix_vars%ShearRichardson_iface => Ri_g

  ! Read / set KPP parameters
  read(*, nml=KPP_nml)
  call cvmix_init_shear(mix_scheme='KPP', KPP_nu_zero=KPP_nu_zero,            &
                        KPP_Ri_zero=KPP_Ri_zero, KPP_exp=KPP_exp)
  call cvmix_coeffs_shear(CVmix_vars)

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"Ri   ", "Tdiff"/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "long_name", "Richardson number",          &
                              var_name="ShearRichardson")
  call cvmix_output_write_att(fid, "units", "unitless",                       &
                              var_name="ShearRichardson")
  call cvmix_output_write_att(fid, "long_name", "temperature diffusivity",    &
                              var_name="Tdiff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="Tdiff")
#endif
  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_shear_driver
