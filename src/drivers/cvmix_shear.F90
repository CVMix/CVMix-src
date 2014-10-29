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

  ! Global parameter
  integer, parameter :: ncol = 6

  ! CVMix datatypes
  type(cvmix_data_type), dimension(ncol) :: CVmix_vars_PP
  type(cvmix_data_type)                  :: CVmix_vars_LMD

  ! Both regression tests look at Richardson numbers in [0,1]
  real(cvmix_r8), dimension(:), allocatable, target :: Ri_g

  ! LMD test is single column, PP test is multi-column
  real(cvmix_r8), dimension(:),   allocatable, target :: LMD_Mdiff, LMD_Tdiff
  real(cvmix_r8), dimension(:,:), allocatable, target :: PP_Mdiff,  PP_Tdiff

  ! Hard-code in parameters for each case
  real(cvmix_r8), dimension(ncol) :: PP_nu_zero, PP_n, PP_alpha

  integer :: icol, kw, fid

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
  allocate(LMD_Mdiff(max_nlev+1), LMD_Tdiff(max_nlev+1))
  allocate(PP_Mdiff(ncol,max_nlev+1), PP_Tdiff(ncol,max_nlev+1))

  ! Set Pacanowski-Philander coefficients
  ! (See Table 1 from paper, converted from cm^2/s to m^2/s)
  PP_nu_zero  = real((/0.02, 0.05, 0.1, 0.1, 0.15, 0.15/), cvmix_r8)
  PP_n(:)     = real(2.0,cvmix_r8)
  PP_alpha(:) = real(5.0,cvmix_r8)
  PP_n(4)     = real(1.0,cvmix_r8)
  PP_alpha(5) = real(10.0,cvmix_r8)

  ! Initialization for CVMix data type
  call cvmix_put(CVmix_vars_LMD, 'nlev',     nlev)
  call cvmix_put(CVmix_vars_LMD, 'max_nlev', max_nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars_LMD%Mdiff_iface           => LMD_Mdiff
  CVmix_vars_LMD%Tdiff_iface           => LMD_Tdiff
  CVmix_vars_LMD%ShearRichardson_iface => Ri_g

  ! Read / set KPP parameters
  read(*, nml=KPP_nml)
  call cvmix_init_shear(mix_scheme='KPP', KPP_nu_zero=KPP_nu_zero,            &
                        KPP_Ri_zero=KPP_Ri_zero, KPP_exp=KPP_exp)
  call cvmix_coeffs_shear(CVmix_vars_LMD)
  do icol=1,ncol
      call cvmix_put(CVmix_vars_PP(icol), 'nlev',     nlev)
      call cvmix_put(CVmix_vars_PP(icol), 'max_nlev', max_nlev)
      CVmix_vars_PP(icol)%Mdiff_iface           => PP_Mdiff(icol,:)
      CVmix_vars_PP(icol)%Tdiff_iface           => PP_Tdiff(icol,:)
      CVmix_vars_PP(icol)%ShearRichardson_iface => Ri_g
      call cvmix_init_shear(mix_scheme='PP', PP_nu_zero=PP_nu_zero(icol),     &
                            PP_alpha=PP_alpha(icol), PP_exp=PP_n(icol))
      call cvmix_coeffs_shear(CVmix_vars_PP(icol))
  end do

  ! Output
  ! (1) LMD column
#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars_LMD, (/"Ri   ", "Tdiff"/))
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

  ! (2) PP columns
#ifdef _NETCDF
  call cvmix_io_open(fid, "data_PP.nc", "nc")
#else
  call cvmix_io_open(fid, "data_PP.out", "ascii")
#endif
  call cvmix_output_write(fid, CVmix_vars_PP, (/"Ri   ", "Tdiff"/))
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
