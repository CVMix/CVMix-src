module vmix_output

!BOP
!\newpage
! !MODULE: vmix_output
!
! !DESCRIPTION:
!  This module contains routines to output CVmix variables to data files.
!  Currently only ascii output is supported, but the plan is to also include
!  plain binary and netCDF output as well.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

#ifdef _NETCDF
   use netcdf
#endif
   use vmix_kinds_and_types, only : vmix_data_type
#ifdef _NETCDF
   use vmix_kinds_and_types, only : vmix_r8
#endif

!EOP

  implicit none
  private
  save

!BOP
! !PUBLIC MEMBER FUNCTIONS:
  public :: vmix_output_open
  public :: vmix_output_write
  public :: vmix_output_close
  public :: print_open_files

  interface vmix_output_write
    module procedure vmix_output_write_single_col
    module procedure vmix_output_write_multi_col
  end interface

! !DEFINED PARAMETERS:
  integer, parameter :: ASCII_FILE_TYPE  = 1
  integer, parameter :: BIN_FILE_TYPE    = 2
  integer, parameter :: NETCDF_FILE_TYPE = 3
  integer, parameter :: FILE_NOT_FOUND   = 404

  ! Probably not the best technique, but going to use a linked list to keep
  ! track of what files are open / what format they are (ascii, bin, or nc)
  type :: vmix_file_entry
    integer :: file_id
    integer :: file_type
    type(vmix_file_entry), pointer :: prev
    type(vmix_file_entry), pointer :: next
  end type

  type(vmix_file_entry), allocatable, target :: file_database(:)
!EOP

contains

!BOP

! !IROUTINE: vmix_output_open
! !INTERFACE:

  subroutine vmix_output_open(file_id, file_name, file_format)

! !DESCRIPTION:
!  Routine to open a file for writing. Goal is to support writing files
!  in plain text (currently working), netCDF, and plain binary. Besides
!  opening the file, this routine also adds an entry to file\_database,
!  a linked list that keeps track of what files are open and what type
!  of file each identifier refers to. So it will be possible to output
!  the same data in ascii and netCDF, for example.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: file_name, file_format

! !OUTPUT PARAMETERS:
    integer, intent(out) :: file_id

! !LOCAL VARIABLES:
    type(vmix_file_entry), pointer :: file_index
!EOP
!BOC

    ! Need routine that will produce unique file_id
    ! Starting with 615 and incrementing by one for now...
    file_id = 615
    if (.not.(allocated(file_database))) then
      allocate(file_database(1))
      file_database(1)%file_id = file_id
      nullify(file_database(1)%prev)
      nullify(file_database(1)%next)
      file_index => file_database(1)
    else
      file_id = file_id+1
      file_index => file_database(1)
      do while(associated(file_index%next))
        file_id = file_id+1
        file_index => file_index%next
      end do
      allocate(file_index%next)
      file_index%next%file_id   = file_id 
      file_index%next%prev     => file_index 
      nullify(file_index%next%next)
      file_index => file_index%next
    end if

    select case (trim(file_format))
      case ('nc')
#ifndef _NETCDF
        print*, "ERROR: you must compile -D_NETCDF to open a netCDF file"
        stop
#else
        file_index%file_type = NETCDF_FILE_TYPE
        call netcdf_check(nf90_create(file_name, NF90_CLOBBER, file_id))
        file_index%file_id = file_id
        ! For outputting params, want vertical dimension to be unlimited?
        ! (Will be looping through the levels)
#endif
      case ('ascii')
        file_index%file_type = ASCII_FILE_TYPE
        open(file_id, file = file_name, status="replace")
      case default
        print*, "ERROR: ", trim(file_format)," is not a valid file type"

    end select
!EOC

  end subroutine vmix_output_open

!BOP

! !IROUTINE: vmix_output_write_single_col
! !INTERFACE:

  subroutine vmix_output_write_single_col(file_id, Vmix_vars, var_names)

! !DESCRIPTION:
!  Routine to write the requested variables from a single column to a file
!  (file must be opened using vmix\_output\_open to ensure it is written
!  correctly). Called with vmix\_output\_write (see interface in PUBLIC
!  MEMBER FUNCTIONS above).
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer,                        intent(in) :: file_id
    type(vmix_data_type),           intent(in) :: Vmix_vars
    character(len=*), dimension(:), intent(in) :: var_names

! !LOCAL VARIABLES:
    integer :: kw, var
#ifdef _NETCDF
    integer                            :: nw, nw_id
    integer, dimension(:), allocatable :: var_id
#endif
!EOP
!BOC

    select case (get_file_type(file_id))
#ifdef _NETCDF
      case (NETCDF_FILE_TYPE)
        nw = Vmix_vars%nlev+1
        call netcdf_check(nf90_def_dim(file_id, "nw", nw, nw_id))
        allocate(var_id(size(var_names)))
        do var=1,size(var_names)
          call netcdf_check(nf90_def_var(file_id, var_names(var), NF90_DOUBLE, &
                                         (/nw_id/), var_id(var)))
        end do
        call netcdf_check(nf90_enddef(file_id))
        do var=1,size(var_names)
          select case(trim(var_names(var)))
            case ("depth")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             Vmix_vars%z_iface(:)))
            case ("Ri")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             Vmix_vars%Ri_iface(:)))
            case ("visc")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             Vmix_vars%visc_iface(:)))
            case ("diff")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             Vmix_vars%diff_iface(:,1)))
            case DEFAULT
              print*, "ERROR: unable to write variable ", var_names(var)
              stop
          end select
        end do
#endif
      case (ASCII_FILE_TYPE)
        do kw=1,Vmix_vars%nlev+1
          do var=1,size(var_names)
            select case(trim(var_names(var)))
              case ("depth")
                write(file_id,"(E24.17E2)",advance='no') &
                      Vmix_vars%z_iface(kw)
              case ("Ri")
                write(file_id,"(E24.17E2)",advance='no') &
                      Vmix_vars%Ri_iface(kw)
              case ("visc")
                write(file_id,"(E24.17E2)",advance='no') &
                      Vmix_vars%visc_iface(kw)
              case ("diff")
                write(file_id,"(E24.17E2)",advance='no') &
                      Vmix_vars%diff_iface(kw,1)
              case DEFAULT
                print*, "ERROR: unable to write variable ", var_names(var)
                stop
            end select
            if (var.ne.size(var_names)) write(file_id, "(1X)", advance='no')
          end do
          write(file_id, *)
        end do
      case DEFAULT
        print*, "ERROR: Invalid file type"
        stop
    end select
!EOC

  end subroutine vmix_output_write_single_col

!BOP

! !IROUTINE: vmix_output_write_multi_col
! !INTERFACE:

  subroutine vmix_output_write_multi_col(file_id, Vmix_vars, var_names)

! !DESCRIPTION:
!  Routine to write the requested variables from multiple columns to a file
!  (file must be opened using vmix\_output\_open to ensure it is written
!  correctly). Called with vmix\_output\_write (see interface in PUBLIC
!  MEMBER FUNCTIONS above).
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer,                            intent(in) :: file_id
    type(vmix_data_type), dimension(:), intent(in) :: Vmix_vars
    character(len=*),     dimension(:), intent(in) :: var_names

! !LOCAL VARIABLES:
    integer :: ncol, nw, icol, kw, var
    logical :: z_err
#ifdef _NETCDF
    integer                                         :: nw_id, ncol_id
    integer, dimension(:), allocatable              :: var_id
    real(kind=vmix_r8), dimension(:,:), allocatable :: lcl_visc, lcl_diff
#endif
!EOP
!BOC

    z_err = .false.
    ncol = size(Vmix_vars)
    nw = Vmix_vars(1)%nlev+1
    ! Make sure all levels are the same
    do icol=2,ncol
      if (Vmix_vars(icol)%nlev+1.ne.nw) then
        z_err = .true.
      else
        if (any(Vmix_vars(icol)%z_iface.ne.Vmix_vars(icol-1)%z_iface)) then
          z_err = .true.
        end if
      end if
    end do

    if (z_err) then
      print*, "ERROR: z-coordinates are not the same in every column!"
      stop
    end if

    select case (get_file_type(file_id))
#ifdef _NETCDF
      case (NETCDF_FILE_TYPE)
        call netcdf_check(nf90_def_dim(file_id, "nw",   nw,   nw_id))
        call netcdf_check(nf90_def_dim(file_id, "ncol", ncol, ncol_id))
        allocate(var_id(size(var_names)))
        do var=1,size(var_names)
          if ((trim(var_names(var)).eq."depth").or.  &
              (trim(var_names(var)).eq."depth")) then
            call netcdf_check(nf90_def_var(file_id, var_names(var),          &
                                NF90_DOUBLE, (/nw_id/), var_id(var)))
          else
            call netcdf_check(nf90_def_var(file_id, var_names(var),          &
                                NF90_DOUBLE, (/ncol_id,nw_id/), var_id(var)))
          end if
          if (trim(var_names(var)).eq."visc") then
            allocate(lcl_visc(ncol,nw))
            do icol=1,ncol
              lcl_visc(icol,:) = Vmix_vars(icol)%visc_iface
            end do
          endif
          if (trim(var_names(var)).eq."diff") then
            allocate(lcl_diff(ncol,nw))
            do icol=1,ncol
              lcl_diff(icol,:) = Vmix_vars(icol)%diff_iface(:,1)
            end do
          endif
        end do
        call netcdf_check(nf90_enddef(file_id))
        do var=1,size(var_names)
          select case(trim(var_names(var)))
            case("depth")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                Vmix_vars(1)%z_iface(:)))
            case("Ri")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                Vmix_vars(1)%Ri_iface(:)))
            case("visc")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                lcl_visc))
              deallocate(lcl_visc)
            case("diff")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                lcl_diff))
              deallocate(lcl_diff)
            case DEFAULT
              print*, "ERROR: unable to write variable ", var_names(var)
              stop
          end select
        end do
#endif
      case (ASCII_FILE_TYPE)
        do kw=1,nw
          do var=1,size(var_names)
            select case(trim(var_names(var)))
              case ("depth")
                write(file_id,"(E24.17E2)",advance='no') &
                      Vmix_vars(1)%z_iface(kw)
              case ("Ri")
                write(file_id,"(E24.17E2)",advance='no') &
                      Vmix_vars(1)%Ri_iface(kw)
              case ("visc")
                do icol=1,ncol
                  write(file_id,"(E24.17E2)",advance='no') &
                        Vmix_vars(icol)%visc_iface(kw)
                  if (icol.ne.ncol) write(file_id, "(1X)", advance='no')
                end do
              case ("diff")
                do icol=1,ncol
                  write(file_id,"(E24.17E2)",advance='no') &
                        Vmix_vars(icol)%diff_iface(kw,1)
                  if (icol.ne.ncol) write(file_id, "(1X)", advance='no')
                end do
              case DEFAULT
                print*, "ERROR: unable to write variable ", var_names(var)
                stop
            end select
            if (var.ne.size(var_names)) write(file_id, "(1X)", advance='no')
          end do
          write(file_id, *)
        end do
      case DEFAULT
        print*, "ERROR: Invalid file type"
        stop
    end select
!EOC

  end subroutine vmix_output_write_multi_col

!BOP

! !IROUTINE: vmix_output_close
! !INTERFACE:

  subroutine vmix_output_close(file_id)

! !DESCRIPTION:
!  Routine to close a file once all writing has been completed. In addition
!  to closing the file, this routine also deletes its entry in file\_database
!  to avoid trying to write to the file in the future.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer, intent(in) :: file_id

! !LOCAL VARIABLES:
    type(vmix_file_entry), pointer :: ifile, file_to_close
    logical :: file_found
    integer :: file_type
!EOP
!BOC

    ! Is fid in the file database?
    nullify(file_to_close)
    if (allocated(file_database)) then
      ifile => file_database(1)
      do while(associated(ifile%next))
        if (ifile%file_id.eq.file_id) then
          file_to_close => ifile
        end if
        ifile => ifile%next
      end do
      if (ifile%file_id.eq.file_id) then
         file_to_close => ifile
      end if
    end if
    file_found = associated(file_to_close)

    if (.not.file_found) then
      write(*,"(A,I0,A)") "Warning: file id ", file_id, " is not an open file!"
      return
    end if
    file_type = file_to_close%file_type

    if (associated(file_to_close%prev)) then
      ifile => file_to_close%prev
      if (associated(file_to_close%next)) then
        ifile%next => file_to_close%next
        ifile%next%prev => ifile
      else
        nullify(ifile%next)
      end if
      deallocate(file_to_close)
    else
      ! file_id is stored in the first entry
      if (associated(file_database(1)%next)) then
        ! Database has more than one entry, so copy last entry into first
        file_to_close => file_database(1)
        do while(associated(file_to_close%next))
          file_to_close => file_to_close%next
        end do
        ifile => file_to_close%prev
        file_database(1)%file_id   = file_to_close%file_id
        file_database(1)%file_type = file_to_close%file_type
        nullify(ifile%next)
        deallocate(file_to_close)
      else
        ! file_id is only entry in database
        deallocate(file_database)
      end if
    end if

    select case (file_type)
#ifdef _NETCDF
      case (NETCDF_FILE_TYPE)
        call netcdf_check(nf90_close(file_id))
#endif
      case (ASCII_FILE_TYPE)
        close(file_id)
      case (BIN_FILE_TYPE)
        close(file_id)
    end select
!EOC

  end subroutine vmix_output_close

!BOP

! !IROUTINE: get_file_type
! !INTERFACE:

  function get_file_type(file_id)

! !DESCRIPTION:
!  Returns the file format (enumerated in DEFINED PARAMETERS section) of a
!  given file. If the file is not in the database, returns FILE\_NOT\_FOUND.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer, intent(in) :: file_id

! !OUTPUT PARAMETERS:
    integer             :: get_file_type

! !LOCAL VARIABLES:
    type(vmix_file_entry), pointer :: ifile
!EOP
!BOC

    ifile => file_database(1)
    if (ifile%file_id.eq.file_id) then
      get_file_type = ifile%file_type
      return
    end if
    do while(associated(ifile%next))
      ifile => ifile%next
      if (ifile%file_id.eq.file_id) then
        get_file_type = ifile%file_type
        return
      end if
    end do
    get_file_type = FILE_NOT_FOUND
!EOC

  end function get_file_type

! DEBUGGING ROUTINE
  subroutine print_open_files()

    type(vmix_file_entry), pointer :: ifile

    if (.not.allocated(file_database)) then
      print*, "No Open files"
    else
      ifile => file_database(1)
      do while (associated(ifile%next))
        print*, "file id: ", ifile%file_id, ifile%file_type
        ifile => ifile%next
      end do
      print*, "file id: ", ifile%file_id, ifile%file_type
    end if
    print*, "----"

  end subroutine print_open_files

  subroutine netcdf_check(status)

    integer, intent(in) :: status

#ifdef _NETCDF
    if (status.ne.nf90_noerr) then
      print*, "netCDF error: ", trim(nf90_strerror(status))
      stop
    end if
#else
    print*, "ERROR: can not call netcdf_check unless compiling -D_NETCDF"
    print*, "The status you passed in = ", status
    stop
#endif

  end subroutine netcdf_check

end module vmix_output

