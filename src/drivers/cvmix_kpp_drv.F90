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
                                    cvmix_strlen,             &
                                    cvmix_data_type
  use cvmix_kpp,             only : cvmix_init_kpp,           &
                                    cvmix_coeffs_kpp
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_output_write,       &
                                    cvmix_output_write_att,   &
                                    cvmix_io_close

  Implicit None

! !INPUT PARAMETERS:
  integer, intent(in) :: nlev

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type)       :: CVmix_vars

  real(cvmix_r8), dimension(:,:), allocatable, target :: diffusivity
  real(cvmix_r8), dimension(:),   allocatable, target :: viscosity
  real(cvmix_r8), dimension(:),   allocatable, target :: zt, zw_iface, Ri_bulk
  integer :: fid, kt, kw
  real(cvmix_r8) :: hmix, ri_crit
  character(len=cvmix_strlen) :: interp_type

  namelist/kpp_nml/hmix, ri_crit, interp_type

  hmix = -15.0_cvmix_r8
  ri_crit = 0.3_cvmix_r8
  interp_type = 'quadratic'
  read(*, nml=kpp_nml)

  allocate(diffusivity(nlev+1,2))
  allocate(viscosity(nlev+1))
  allocate(zt(nlev), zw_iface(nlev+1), Ri_bulk(nlev))
  do kw=1,nlev+1
    zw_iface(kw) = real(-10*(kw-1), cvmix_r8)
  end do
  do kt=1,nlev
    zt(kt) = 0.5_cvmix_r8*(zw_iface(kt)+zw_iface(kt+1))
    if (zt(kt).gt.hmix) then
      Ri_bulk(kt) = 0.0_cvmix_r8
    else
      if (Ri_bulk(kt-1).eq.0) then
        ! Exact integration for average value over first cell with non-zero
        ! Ri_bulk
        Ri_bulk(kt) = 0.025_cvmix_r8*ri_crit*(zw_iface(kt+1)-hmix)**2
      else
        Ri_bulk(kt) = 0.5_cvmix_r8*ri_crit*(hmix-zt(kt))
      end if
    end if
    print*, zt(kt), Ri_bulk(kt)
  end do

  call cvmix_put(CVmix_vars, 'nlev', nlev)
  call cvmix_put(CVmix_vars, 'ocn_depth', 10*nlev)
  CVmix_vars%diff_iface => diffusivity(:,:)
  CVmix_vars%visc_iface => viscosity(:)
  CVmix_vars%zt         => zt(:)
  CVmix_vars%zw_iface   => zw_iface(:)
  CVmix_vars%Rib        => Ri_bulk(:)

  call cvmix_init_kpp(ri_crit=ri_crit, interp_type=interp_type)
  call cvmix_coeffs_kpp(CVmix_vars)
  print*, "OBL depth = ", CVmix_vars%OBL_depth

#ifdef _NETCDF
  call cvmix_io_open(fid, "data.nc", "nc")
#else
  call cvmix_io_open(fid, "data.out", "ascii")
#endif

  call cvmix_output_write(fid, CVmix_vars, (/"zt     ", "zw     ", "Ri_bulk"/))
#ifdef _NETCDF
  call cvmix_output_write_att(fid, "Interpolation", interp_type)
  call cvmix_output_write_att(fid, "OBL_depth", CVmix_vars%OBL_depth)
#endif

  call cvmix_io_close(fid)

!EOC

End Subroutine cvmix_kpp_driver
