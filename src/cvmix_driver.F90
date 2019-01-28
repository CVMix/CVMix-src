!BOP

! \newpage
! !MODULE: Main Program (Stand-Alone)
! !ROUTINE: cvmix_driver

! !DESCRIPTION: The stand-alone driver for the CVMix package. This driver
!  reads in the cvmix\_nml namelist to determine what type of mixing has
!  been requested, and also reads in mixing-specific parameters from a
!  mixingtype\_nml namelist.
!\\
!\\

! !INTERFACE:

Program cvmix_driver

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_zero,                               &
                                    cvmix_strlen

!EOP
!BOC
  integer :: nlev, max_nlev
  real(kind=cvmix_r8) :: ocn_depth
  character(len=cvmix_strlen) :: mix_type

  namelist/cvmix_nml/mix_type, nlev, max_nlev, ocn_depth

  mix_type = 'unset'
  nlev = 0
  max_nlev = 0
  ocn_depth = cvmix_zero
  read(*, nml=cvmix_nml)
  if (max_nlev.eq.0) then
    max_nlev = nlev
  end if

  select case (trim(mix_type))
    case ('BryanLewis')
      call cvmix_BL_driver(nlev, max_nlev, ocn_depth)
    case ('shear')
      call cvmix_shear_driver(nlev, max_nlev)
    case ('tidal')
      call cvmix_tidal_driver()
    case ('ddiff')
      call cvmix_ddiff_driver(2*nlev, 2*max_nlev)
    case ('kpp')
      call cvmix_kpp_driver()
    case DEFAULT
      print*, "WARNING: mix_type = '", trim(mix_type), "' is not supported by this driver."
  end select

End Program cvmix_driver
