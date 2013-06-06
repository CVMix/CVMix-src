!BOI

! !TITLE: In-code documentation for CVMix
! !AUTHORS: Many contributors from GFDL, LANL, and NCAR
! !AFFILIATION: GFDL, LANL, and NCAR
! !DATE: \today

!EOI
!BOP

! !ROUTINE: cvmix_driver

! !DESCRIPTION: The stand-alone driver for the CVMix package. This driver
!  reads in the cvmix\_nml namelist to determine what type of mixing has
!  been requested, and also reads in mixing-specific parameters from a
!  mixingtype\_nml namelist.
!\\
!\\

! !REVISION HISTORY:
!  SVN $Id$
!  SVN $URL$

! !INTERFACE:

Program cvmix_driver

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8, &
                                    cvmix_strlen

!EOP
!BOC
  integer :: nlev
  real(kind=cvmix_r8) :: ocn_depth 
  character(len=cvmix_strlen) :: mix_type

  namelist/cvmix_nml/mix_type, nlev, ocn_depth

  mix_type = 'unset'
  nlev = 0
  ocn_depth = 0.0_cvmix_r8
  read(*, nml=cvmix_nml)

  select case (trim(mix_type))
    case ('BryanLewis_pointer')
      call cvmix_BL_pointer_driver(nlev, ocn_depth)
    case ('BryanLewis_memcopy')
      call cvmix_BL_memcopy_driver(nlev, ocn_depth)
    case ('shear')
      call cvmix_shear_driver(nlev)
    case ('tidal')
      call cvmix_tidal_driver()
    case ('ddiff')
      call cvmix_ddiff_driver(nlev)
    case ('kpp')
      call cvmix_kpp_driver(nlev)
    case DEFAULT
      print*, "WARNING: mix_type = '", trim(mix_type), "' is not supported by this driver."
  end select

End Program cvmix_driver
