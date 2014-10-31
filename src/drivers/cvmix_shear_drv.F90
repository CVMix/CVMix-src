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
  type(cvmix_data_type)                  :: CVmix_vars_LMD_1D, CVmix_vars_PP_1D
  type(cvmix_data_type), dimension(ncol) :: CVmix_vars_PP_2D

  ! All regression tests look at Richardson numbers in [0,1]
  real(cvmix_r8), dimension(:), allocatable, target :: Ri_g

  ! "1D" variables will be 2 x nlev+1 (1 for LMD, 1 for PP)
  ! "2D" variables will be ncol x nlev+1 (for PP)
  real(cvmix_r8), dimension(:,:), allocatable, target :: Mdiff_1D, Tdiff_1D
  real(cvmix_r8), dimension(:,:), allocatable, target :: Mdiff_2D, Tdiff_2D

  ! Hard-code in parameters for each c in Table 1
  real(cvmix_r8), dimension(ncol) :: PP_nu_zero_2D, PP_exp_2D, PP_alpha_2D
  real(cvmix_r8)                  :: PP_nu_zero, PP_exp, PP_alpha

  integer :: icol, kw, fid

  ! Namelist variables
  ! KPP mixing parameters for column
  real(cvmix_r8) :: LMD_nu_zero, LMD_Ri_zero, LMD_exp

  namelist/LMD_nml/LMD_nu_zero, LMD_Ri_zero, LMD_exp
  namelist/PP_nml/PP_nu_zero, PP_alpha, PP_exp

  print*, "Active levels: ", nlev
  print*, "Levels allocated in memory: ", max_nlev

  ! Namelist (set defaults then read from file)
  LMD_nu_zero = 5e-3_cvmix_r8
  LMD_Ri_zero = 0.7_cvmix_r8
  LMD_exp     = real(3,  cvmix_r8)

  PP_nu_zero = 5e-3_cvmix_r8
  PP_alpha   = real(5, cvmix_r8)
  PP_exp     = real(2, cvmix_r8)

  read(*, nml=LMD_nml)
  read(*, nml=PP_nml)

  print*, ""
  print*, "Parameters Used in LMD test"
  print*, "----"
  print*, "KPP_nu_zero = ", LMD_nu_zero
  print*, "KPP_Ri_zero = ", LMD_Ri_zero
  print*, "KPP_exp = ", LMD_exp

  print*, ""
  print*, "Parameters Used in PP test"
  print*, "----"
  print*, "PP_nu_zero = ", PP_nu_zero
  print*, "PP_alpha = ", PP_alpha
  print*, "PP_exp = ", PP_exp


  ! Ri_g should increase from 0 to 1 in active portion of level
  allocate(Ri_g(max_nlev+1))
  Ri_g(1) = cvmix_zero
  do kw = 2,max_nlev+1
    Ri_g(kw) = Ri_g(kw-1) + cvmix_one/real(nlev,cvmix_r8)
  end do

  ! Allocate memory to store viscosity and diffusivity values
  allocate(Mdiff_1D(2,max_nlev+1), Tdiff_1D(2,max_nlev+1))
  allocate(Mdiff_2D(ncol,max_nlev+1), Tdiff_2D(ncol,max_nlev+1))

  ! Set Pacanowski-Philander coefficients
  ! (See Table 1 from paper, converted from cm^2/s to m^2/s)
  PP_nu_zero_2D  = 0.001_cvmix_r8 * real((/2, 5, 10, 10, 15, 15/), cvmix_r8)
  PP_exp_2D(:)   = real(2,cvmix_r8)
  PP_alpha_2D(:) = real(5,cvmix_r8)
  PP_exp_2D(4)   = real(1,cvmix_r8)
  PP_alpha_2D(6) = real(10,cvmix_r8)

  ! Initialization for LMD test
  call cvmix_put(CVmix_vars_LMD_1D, 'nlev',     nlev)
  call cvmix_put(CVmix_vars_LMD_1D, 'max_nlev', max_nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars_LMD_1D%Mdiff_iface           => Mdiff_1D(1,:)
  CVmix_vars_LMD_1D%Tdiff_iface           => Tdiff_1D(1,:)
  CVmix_vars_LMD_1D%ShearRichardson_iface => Ri_g

  ! Initialization for 1D PP test
  call cvmix_put(CVmix_vars_PP_1D, 'nlev',     nlev)
  call cvmix_put(CVmix_vars_PP_1D, 'max_nlev', max_nlev)
  ! Point CVmix_vars values to memory allocated above
  CVmix_vars_PP_1D%Mdiff_iface           => Mdiff_1D(2,:)
  CVmix_vars_PP_1D%Tdiff_iface           => Tdiff_1D(2,:)
  CVmix_vars_PP_1D%ShearRichardson_iface => Ri_g

  ! Initialization for 2D PP test
  do icol=1,ncol
    call cvmix_put(CVmix_vars_PP_2D(icol), 'nlev',     nlev)
    call cvmix_put(CVmix_vars_PP_2D(icol), 'max_nlev', max_nlev)
    CVmix_vars_PP_2D(icol)%Mdiff_iface           => Mdiff_2D(icol,:)
    CVmix_vars_PP_2D(icol)%Tdiff_iface           => Tdiff_2D(icol,:)
    CVmix_vars_PP_2D(icol)%ShearRichardson_iface => Ri_g
  end do

  ! Set LMD94 parameters
  call cvmix_init_shear(mix_scheme='KPP', KPP_nu_zero=LMD_nu_zero,            &
                        KPP_Ri_zero=LMD_Ri_zero, KPP_exp=LMD_exp)
  call cvmix_coeffs_shear(CVmix_vars_LMD_1D)

  ! Set PP81 single column parameters
  call cvmix_init_shear(mix_scheme='PP', PP_nu_zero=PP_nu_zero,               &
                        PP_alpha=PP_alpha, PP_exp=PP_exp)
  call cvmix_coeffs_shear(CVmix_vars_PP_1D)

  ! Set PP81 multiple column
  do icol=1,ncol
    call cvmix_init_shear(mix_scheme='PP', PP_nu_zero=PP_nu_zero_2D(icol),    &
                          PP_alpha=PP_alpha_2D(icol), PP_exp=PP_exp_2D(icol))
    call cvmix_coeffs_shear(CVmix_vars_PP_2D(icol))
  end do

  ! Output
  ! (1) LMD column
#ifdef _NETCDF
  call cvmix_io_open(fid, "data_LMD.nc", "nc")
#else
  call cvmix_io_open(fid, "data_LMD.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars_LMD_1D, (/"Ri   ", "Tdiff"/))
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

  ! (2) PP single column
#ifdef _NETCDF
  call cvmix_io_open(fid, "data_PP1d.nc", "nc")
#else
  call cvmix_io_open(fid, "data_PP1d.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars_PP_1D, (/"Ri   ", "Mdiff"/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "long_name", "Richardson number",          &
                              var_name="ShearRichardson")
  call cvmix_output_write_att(fid, "units", "unitless",                       &
                              var_name="ShearRichardson")
  call cvmix_output_write_att(fid, "long_name", "momentum diffusivity",       &
                              var_name="Mdiff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="Mdiff")
#endif
  call cvmix_io_close(fid)

  ! (3) PP multiple columns
#ifdef _NETCDF
  call cvmix_io_open(fid, "data_PP2d.nc", "nc")
#else
  call cvmix_io_open(fid, "data_PP2d.out", "ascii")
#endif
  call cvmix_output_write(fid, CVmix_vars_PP_2D, (/"Ri   ", "Mdiff"/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "long_name", "Richardson number",          &
                              var_name="ShearRichardson")
  call cvmix_output_write_att(fid, "units", "unitless",                       &
                              var_name="ShearRichardson")
  call cvmix_output_write_att(fid, "long_name", "momentum diffusivity",       &
                              var_name="Mdiff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="Mdiff")
#endif
  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_shear_driver
