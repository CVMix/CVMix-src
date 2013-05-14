module cvmix_io

!BOP
!\newpage
! !MODULE: cvmix_io
!
! !DESCRIPTION:
!  This module contains routines to read CVmix variables from data files or 
!  output CVmix variables to data files. Currently only ascii and netCDF output
!  are supported, as well as netCDF input, but the plan is to also include plain
!  binary input / output as well.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use cvmix_kinds_and_types, only : cvmix_data_type, &
                                     cvmix_r8,        &
                                     cvmix_strlen
#ifdef _NETCDF
   use netcdf
#endif

!EOP

  implicit none
  private
  save

!BOP
! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_io_open
  public :: cvmix_input_read
  public :: cvmix_output_write
  public :: cvmix_io_close
  public :: cvmix_io_close_all
  public :: print_open_files

  interface cvmix_input_read
    module procedure cvmix_input_read_2d_double
  end interface

  interface cvmix_output_write
    module procedure cvmix_output_write_single_col
    module procedure cvmix_output_write_multi_col
  end interface

! !DEFINED PARAMETERS:
  integer, parameter :: ASCII_FILE_TYPE  = 1
  integer, parameter :: BIN_FILE_TYPE    = 2
  integer, parameter :: NETCDF_FILE_TYPE = 3
  integer, parameter :: FILE_NOT_FOUND   = 404

  ! Probably not the best technique, but going to use a linked list to keep
  ! track of what files are open / what format they are (ascii, bin, or nc)
  type :: cvmix_file_entry
    integer :: file_id
    integer :: file_type
    character(len=cvmix_strlen) :: file_name
    type(cvmix_file_entry), pointer :: prev
    type(cvmix_file_entry), pointer :: next
  end type

  type(cvmix_file_entry), allocatable, target :: file_database(:)
!EOP

contains

!BOP

! !IROUTINE: cvmix_io_open
! !INTERFACE:

  subroutine cvmix_io_open(file_id, file_name, file_format, read_only)

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
    character(len=*),  intent(in) :: file_name, file_format
    logical, optional, intent(in) :: read_only

! !OUTPUT PARAMETERS:
    integer, intent(out) :: file_id

! !LOCAL VARIABLES:
    type(cvmix_file_entry), pointer :: file_index
    logical                         :: readonly
!EOP
!BOC

    if (present(read_only)) then
      readonly = read_only
    else
      readonly = .false.
    end if
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
    file_index%file_name = trim(file_name)

    select case (trim(file_format))
      case ('nc')
#ifndef _NETCDF
        print*, "ERROR: you must compile -D_NETCDF to open a netCDF file"
        stop
#else
        file_index%file_type = NETCDF_FILE_TYPE
        ! Note: at this point, will either open file with NOWRITE for
        !       read-only, or will clobber file to write new data to it.
        !       Eventually we should add a check to see if the file exists
        !       and open it with NF90_WRITE for non-readonly files, but that
        !       will require checking to see if dims / variables already exist
        !       (and are correct dimension) before trying to define them.
        if (readonly) then
          call netcdf_check(nf90_open(file_name, NF90_NOWRITE, file_id))
        else
          call netcdf_check(nf90_create(file_name, NF90_CLOBBER, file_id))
        end if
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

  end subroutine cvmix_io_open

!BOP

! !IROUTINE: cvmix_input_read_2d_double
! !INTERFACE:

  subroutine cvmix_input_read_2d_double(file_id, var_name, local_copy)

! !DESCRIPTION:
!  Routine to read the requested 2D variable from a netcdf file and save it to
!  a local array (file must be opened using cvmix\_io\_open with the optional
!  argument readonly = .true.). Called with cvmix\_input\_read (see interface
!  in PUBLIC MEMBER FUNCTIONS above). At this time, only works with netcdf
!  files.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer,          intent(in)  :: file_id
    character(len=*), intent(in)  :: var_name
    real(cvmix_r8), dimension(:,:),  intent(out) :: local_copy

! !LOCAL VARIABLES:
    logical :: lerr_in_read
#ifdef _NETCDF
    integer :: varid, nvar, i, ndims, xtype
    integer, dimension(2) :: dims1, dims2
    character(len=cvmix_strlen) :: tmp_name
#endif

!EOP
!BOC

  local_copy = 0.0_cvmix_r8
  lerr_in_read = .false.
    select case (get_file_type(file_id))
#ifdef _NETCDF
      case (NETCDF_FILE_TYPE)
        varid = -1
        ! Find number of variables in file
        call netcdf_check(nf90_inquire(file_id, nVariables=nvar))
        i = 1
        do while((i.le.nvar).and.(varid.eq.-1))
          ! Loop to figure out if var_name is a valid variable in the file
          call netcdf_check(nf90_inquire_variable(file_id, i, name=tmp_name,&
                                                  xtype=xtype, ndims=ndims))
          if (trim(var_name).eq.trim(tmp_name)) then
            varid = i
          else
            i = i+1
          end if
        end do
        lerr_in_read = (varid.eq.-1)

        if (lerr_in_read) then
          write(*,"(A,A,1X,A,A)") "Could not find variable ", trim(var_name), &
                                  "in ", trim(get_file_name(file_id))
        else
          ! A couple more error checks
          if (xtype.ne.NF90_DOUBLE) then
            write(*, "(A,1X,A,1X,A)") "ERROR: variable", trim(var_name), &
                                      "is not a double-precision float!"
            lerr_in_read = .true.
          end if
          if (ndims.ne.2) then
            write(*,"(A,1X,I0,A)") "ERROR: you are trying to read a", ndims, &
                                   "-dimensional array into a 2D array."
            lerr_in_read = .true.
          end if
        end if

        
        if (.not.lerr_in_read) then
          call netcdf_check(nf90_inquire_variable(file_id, varid, dimids=dims1))
          do i=1,2
            call netcdf_check(nf90_inquire_dimension(file_id, dims1(i), &
                              len=dims2(i)))
          end do

          dims1 = shape(local_copy)
          if (all(dims1.eq.dims2)) then
            call netcdf_check(nf90_get_var(file_id, varid, local_copy))
          else
            write(*,"(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
                    "ERROR: you are trying to read a", dims2(1), "by", dims2(2), &
                    "array into a local variable that is", dims1(1), "by", dims1(2)
            lerr_in_read = .true.
          end if
        end if
#endif
      case DEFAULT
        lerr_in_read = .true.
        write(*,"(A,1X,A,1X,A)") "ERROR: no read support for binary files,", &
                                 "use netCDF to read", trim(var_name)
    end select

    if (lerr_in_read) then
      call cvmix_io_close_all
      stop 1
    end if
!EOC

  end subroutine cvmix_input_read_2d_double

!BOP

! !IROUTINE: cvmix_output_write_single_col
! !INTERFACE:

  subroutine cvmix_output_write_single_col(file_id, CVmix_vars, var_names)

! !DESCRIPTION:
!  Routine to write the requested variables from a single column to a file
!  (file must be opened using cvmix\_io\_open to ensure it is written
!  correctly). Called with cvmix\_output\_write (see interface in PUBLIC
!  MEMBER FUNCTIONS above).
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer,                        intent(in) :: file_id
    type(cvmix_data_type),           intent(in) :: CVmix_vars
    character(len=*), dimension(:), intent(in) :: var_names

! !LOCAL VARIABLES:
    integer :: kw, var
#ifdef _NETCDF
    integer                            :: nt, nt_id, nw, nw_id
    integer, dimension(:), allocatable :: var_id
#endif
!EOP
!BOC

    select case (get_file_type(file_id))
#ifdef _NETCDF
      case (NETCDF_FILE_TYPE)
        nt = CVmix_vars%nlev
        nw = CVmix_vars%nlev+1
        call netcdf_check(nf90_def_dim(file_id, "nt", nt, nt_id))
        call netcdf_check(nf90_def_dim(file_id, "nw", nw, nw_id))
        allocate(var_id(size(var_names)))
        do var=1,size(var_names)
          if (trim(var_names(var)).eq."Rrho") then
            call netcdf_check(nf90_def_var(file_id, var_names(var), NF90_DOUBLE, &
                                           (/nt_id/), var_id(var)))
          else
            call netcdf_check(nf90_def_var(file_id, var_names(var), NF90_DOUBLE, &
                                           (/nw_id/), var_id(var)))
          end if
        end do
        call netcdf_check(nf90_enddef(file_id))
        do var=1,size(var_names)
          select case(trim(var_names(var)))
            case ("depth")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             CVmix_vars%z_iface(:)))
            case ("Ri")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             CVmix_vars%Ri_iface(:)))
            case ("visc")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             CVmix_vars%visc_iface(:)))
            case ("diff")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             CVmix_vars%diff_iface(:,1)))
            case ("Rrho")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                             CVmix_vars%strat_param_num(:)/ &
                                             CVmix_vars%strat_param_denom(:)))
            case DEFAULT
              print*, "ERROR: unable to write variable ", var_names(var)
              stop
          end select
        end do
#endif
      case (ASCII_FILE_TYPE)
        do kw=1,CVmix_vars%nlev+1
          do var=1,size(var_names)
            select case(trim(var_names(var)))
              case ("depth")
                write(file_id,"(E24.17E2)",advance='no') &
                      CVmix_vars%z_iface(kw)
              case ("Ri")
                write(file_id,"(E24.17E2)",advance='no') &
                      CVmix_vars%Ri_iface(kw)
              case ("visc")
                write(file_id,"(E24.17E2)",advance='no') &
                      CVmix_vars%visc_iface(kw)
              case ("diff")
                write(file_id,"(E24.17E2)",advance='no') &
                      CVmix_vars%diff_iface(kw,1)
              case ("Rrho")
                if (kw<CVmix_vars%nlev+1) then
                  write(file_id,"(E24.17E2)",advance='no') &
                        CVmix_vars%strat_param_num(kw) / &
                        CVmix_vars%strat_param_denom(kw)
                else
                  write(file_id,"(E24.17E2)",advance='no') 0.0
                end if
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

  end subroutine cvmix_output_write_single_col

!BOP

! !IROUTINE: cvmix_output_write_multi_col
! !INTERFACE:

  subroutine cvmix_output_write_multi_col(file_id, CVmix_vars, var_names)

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
    integer,                             intent(in) :: file_id
    type(cvmix_data_type), dimension(:), intent(in) :: CVmix_vars
    character(len=*),      dimension(:), intent(in) :: var_names

! !LOCAL VARIABLES:
    integer :: ncol, nw, nt, icol, kw, var
    logical :: z_err
#ifdef _NETCDF
    integer                                          :: nt_id, nw_id, ncol_id
    integer,             dimension(:),   allocatable :: var_id
    real(kind=cvmix_r8), dimension(:,:), allocatable :: lcl_visc, lcl_diff, lcl_Rrho
#endif
!EOP
!BOC

    z_err = .false.
    ncol = size(CVmix_vars)
    nt = CVmix_vars(1)%nlev
    nw = CVmix_vars(1)%nlev+1
    ! Make sure all levels are the same
    do icol=2,ncol
      if (CVmix_vars(icol)%nlev+1.ne.nw) then
        z_err = .true.
      else
        ! Make sure z_iface lines up for Bryan-Lewis case
        if (associated(CVmix_vars(1)%z_iface)) then
          if (any(CVmix_vars(icol)%z_iface.ne.CVmix_vars(icol-1)%z_iface)) then
            z_err = .true.
          end if
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
        call netcdf_check(nf90_def_dim(file_id, "nt",   nt,   nt_id))
        call netcdf_check(nf90_def_dim(file_id, "nw",   nw,   nw_id))
        call netcdf_check(nf90_def_dim(file_id, "ncol", ncol, ncol_id))
        allocate(var_id(size(var_names)))
        do var=1,size(var_names)
          if (trim(var_names(var)).eq."depth") then
            call netcdf_check(nf90_def_var(file_id, var_names(var),          &
                                NF90_DOUBLE, (/nw_id/), var_id(var)))
          else
            if (trim(var_names(var)).eq."Rrho") then
              call netcdf_check(nf90_def_var(file_id, var_names(var),        &
                                  NF90_DOUBLE, (/ncol_id,nt_id/), var_id(var)))
            else
              call netcdf_check(nf90_def_var(file_id, var_names(var),        &
                                  NF90_DOUBLE, (/ncol_id,nw_id/), var_id(var)))
            end if
          end if
          if (trim(var_names(var)).eq."visc") then
            allocate(lcl_visc(ncol,nw))
            do icol=1,ncol
              lcl_visc(icol,:) = CVmix_vars(icol)%visc_iface
            end do
          endif
          if (trim(var_names(var)).eq."diff") then
            allocate(lcl_diff(ncol,nw))
            do icol=1,ncol
              lcl_diff(icol,:) = CVmix_vars(icol)%diff_iface(:,1)
            end do
          endif
          if (trim(var_names(var)).eq."Rrho") then
            allocate(lcl_Rrho(ncol,nt))
            do icol=1,ncol
              lcl_Rrho(icol,:) = CVmix_vars(icol)%strat_param_num(:) / &
                                 CVmix_vars(icol)%strat_param_denom(:)
            end do
          endif
        end do
        call netcdf_check(nf90_enddef(file_id))
        do var=1,size(var_names)
          select case(trim(var_names(var)))
            case("depth")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                CVmix_vars(1)%z_iface(:)))
            case("Ri")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                CVmix_vars(1)%Ri_iface(:)))
            case("visc")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                lcl_visc))
              deallocate(lcl_visc)
            case("diff")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                lcl_diff))
              deallocate(lcl_diff)
            case("Rrho")
              call netcdf_check(nf90_put_var(file_id, var_id(var), &
                                lcl_Rrho))
              deallocate(lcl_Rrho)
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
                      CVmix_vars(1)%z_iface(kw)
              case ("Ri")
                write(file_id,"(E24.17E2)",advance='no') &
                      CVmix_vars(1)%Ri_iface(kw)
              case ("visc")
                do icol=1,ncol
                  write(file_id,"(E24.17E2)",advance='no') &
                        CVmix_vars(icol)%visc_iface(kw)
                  if (icol.ne.ncol) write(file_id, "(1X)", advance='no')
                end do
              case ("diff")
                do icol=1,ncol
                  write(file_id,"(E24.17E2)",advance='no') &
                        CVmix_vars(icol)%diff_iface(kw,1)
                  if (icol.ne.ncol) write(file_id, "(1X)", advance='no')
                end do
              case ("Rrho")
                do icol=1,ncol
                  if (kw.ne.nw) then
                    write(file_id,"(E24.17E2)",advance='no')     &
                          CVmix_vars(icol)%strat_param_num(kw) / &
                          CVmix_vars(icol)%strat_param_denom(kw)
                  else
                    write(file_id,"(E24.17E2)",advance='no') 0.0
                  end if
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

  end subroutine cvmix_output_write_multi_col

!BOP

! !IROUTINE: cvmix_io_close
! !INTERFACE:

  subroutine cvmix_io_close_all

! !DESCRIPTION:
!  Routine to close all files open (meant to be called prior to an abort)
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !LOCAL VARIABLES:
    integer :: fid

!EOP
!BOC

    write(*,"(A)") "Closing all open files..."
    do while (allocated(file_database))
      fid = file_database(1)%file_id
      write(*, "(A,1X,A)") "...", trim(get_file_name(fid))
      call cvmix_io_close(fid)
    end do
    write(*,"(A)") "All files closed."
!EOC
  end subroutine cvmix_io_close_all

!BOP

! !IROUTINE: cvmix_io_close
! !INTERFACE:

  subroutine cvmix_io_close(file_id)

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
    type(cvmix_file_entry), pointer :: ifile, file_to_close
    logical                         :: file_found
    integer                         :: file_type
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
        file_database(1)%file_name = file_to_close%file_name
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

  end subroutine cvmix_io_close

!BOP

! !IROUTINE: get_file_name
! !INTERFACE:

  function get_file_name(file_id)

! !DESCRIPTION:
!  Returns the name of the file associated with a given file\_id. If the file
!  is not in the database, returns FILE\_NOT\_FOUND.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    integer, intent(in) :: file_id

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: get_file_name

! !LOCAL VARIABLES:
    type(cvmix_file_entry), pointer :: ifile
!EOP
!BOC

    ifile => file_database(1)
    if (ifile%file_id.eq.file_id) then
      get_file_name = ifile%file_name
      return
    end if
    do while(associated(ifile%next))
      ifile => ifile%next
      if (ifile%file_id.eq.file_id) then
        get_file_name = ifile%file_name
        return
      end if
    end do
    get_file_name = "FILE_NOT_FOUND"
!EOC

  end function get_file_name

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
    type(cvmix_file_entry), pointer :: ifile
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

! Routine to handle errors returned from netcdf
  subroutine netcdf_check(status)

    integer, intent(in) :: status

#ifdef _NETCDF
    if (status.ne.nf90_noerr) then
      print*, "netCDF error: ", trim(nf90_strerror(status))
      stop 1
    end if
#else
    print*, "ERROR: can not call netcdf_check unless compiling -D_NETCDF"
    print*, "The status you passed in = ", status
    stop
#endif

  end subroutine netcdf_check

! DEBUGGING ROUTINE
  subroutine print_open_files()

    type(cvmix_file_entry), pointer :: ifile

    if (.not.allocated(file_database)) then
      print*, "No Open files"
    else
      ifile => file_database(1)
      do while (associated(ifile%next))
        print*, "file id: ", ifile%file_id, ifile%file_type, trim(ifile%file_name)
        ifile => ifile%next
      end do
      print*, "file id: ", ifile%file_id, ifile%file_type, trim(ifile%file_name)
    end if
    print*, "----"

  end subroutine print_open_files

end module cvmix_io

