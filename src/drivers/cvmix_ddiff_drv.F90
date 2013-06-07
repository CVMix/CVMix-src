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

  use cvmix_kinds_and_types, only : one,                     &
                                    cvmix_r8,                &
                                    cvmix_data_type
  use cvmix_ddiff,           only : cvmix_init_ddiff,        &
                                    cvmix_coeffs_ddiff,      &
                                    cvmix_get_ddiff_real,    &
                                    cvmix_ddiff_params_type
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,           &
                                    cvmix_output_write,      &
#ifdef _NETCDF
                                    cvmix_output_write_att,  &
#endif
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type),         dimension(2) :: CVmix_vars
  type(cvmix_ddiff_params_type)               :: CVmix_ddiff_params

  real(cvmix_r8), dimension(:,:,:), allocatable, target :: diffusivity
  real(cvmix_r8), dimension(:,:),   allocatable, target :: Rrho_num, Rrho_denom

  ! column / file indices
  integer :: k, fid, ncol, ic

  ! Namelist variables
  real(cvmix_r8) :: ddiff_exp1, strat_param_max, kappa_ddiff_t

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/ddiff_nml/ddiff_exp1, strat_param_max, kappa_ddiff_t

  ! Allocate memory to store diffusivity values
  ! Also store numerator / denominator for stratification parameter
  ncol = 2
  allocate(diffusivity(nlev+1,2,ncol))
  allocate(Rrho_num(nlev,ncol), Rrho_denom(nlev,ncol))
  do k=1,nlev
    ! For first column, Rrho varies from 1 to 2
    Rrho_num(k,1) = real(k-1,cvmix_r8)/real(nlev-1,cvmix_r8)+one
    Rrho_denom(k,1) = one
    ! For second column, 1/Rrho varies from 1 to 10
    ! (Note: last column has diff=0, hence only using nlev instead of nlev+1)
    Rrho_num(k,2) = -one/(real(9.0*(k-1),cvmix_r8)/real(nlev-1,cvmix_r8)+one)
    Rrho_denom(k,2) = -one
  end do

  ! Point CVmix_vars values to memory allocated above
  do ic=1,ncol
    call cvmix_put(CVmix_vars(ic), 'nlev', nlev)
    CVmix_vars(ic)%diff_iface => diffusivity(:,:,ic)
    CVmix_vars(ic)%strat_param_num => Rrho_num(:,ic)
    CVmix_vars(ic)%strat_param_denom => Rrho_denom(:,ic)
  end do

  ! Read / set double diffusion parameters
  read(*, nml=ddiff_nml)
  call cvmix_init_ddiff(CVmix_ddiff_params,'mks', ddiff_exp1=ddiff_exp1, &
                        strat_param_max=strat_param_max,                 &
                        kappa_ddiff_t=kappa_ddiff_t)

  call cvmix_coeffs_ddiff(CVmix_vars(1), CVmix_ddiff_params)
  call cvmix_coeffs_ddiff(CVmix_vars(2), CVmix_ddiff_params)
  ! For continuity of plot, set diffusivity when Rrho = 1
 !diffusivity(1,1,1) = CVmix_ddiff_params%kappa_ddiff_t
 !diffusivity(1,1,2) = CVmix_ddiff_params%mol_diff*&
 !                     CVmix_ddiff_params%kappa_ddiff_param1*&
 !                     exp(CVmix_ddiff_params%kappa_ddiff_param2)
  diffusivity(1,1,1) = kappa_ddiff_t
  diffusivity(1,1,2) = cvmix_get_ddiff_real(CVmix_ddiff_params,'mol_diff')*&
                       cvmix_get_ddiff_real(CVmix_ddiff_params,'kappa_ddiff_param1')*&
                       exp(cvmix_get_ddiff_real(CVmix_ddiff_params,'kappa_ddiff_param2'))

  ! Output
  ! data will have diffusivity from both columns (needed for NCL script)
#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"Rrho", "diff"/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "long_name", "double diffusion " //        &
                              "stratification parameter", var_name="Rrho")
  call cvmix_output_write_att(fid, "units", "unitless", var_name="Rrho")
  call cvmix_output_write_att(fid, "long_name", "tracer diffusivity",         &
                              var_name="diff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="diff")
#endif
  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_ddiff_driver
