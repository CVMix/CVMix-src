module vmix_output

!BOP
! !MODULE: vmix_output
!
! !DESCRIPTION:
!  This module contains routines to output CVmix variables, either to a netCDF
!  file or an ascii file.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id: vmix_background.F90 39453 2012-08-14 20:01:16Z mlevy@ucar.edu $
!  SVN:$URL: https://svn-ccsm-models.cgd.ucar.edu/pop2/branches/vmix_project/source/vmix/vmix_background.F90 $

! !USES:

#ifdef _NETCDF
   use netcdf
#endif
   use vmix_kinds_and_types
!EOP

  implicit none
  private
  save

!BOP
  public :: vmix_output_open
!EOP

contains

!BOP

! !IROUTINE: vmix_output_open
! !INTERFACE:

  subroutine vmix_output_open(file_format, file_id)

! !DESCRIPTION:
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
  character(len=*), intent(in) :: file_format
  integer,          intent(in) :: file_id

! !OUTPUT PARAMETERS:

!EOP
!BOC

  select case (trim(file_format))
    case ('nc')
#ifndef _NETCDF
        print*, "ERROR: you must compile -D_NETCDF to open a netCDF file"
        stop
#endif

    case default
      print*, "ERROR: ", trim(file_format)," is not a valid file type"

  end select
!EOC
  end subroutine vmix_output_open

end module vmix_output

