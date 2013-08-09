module cvmix_background

!BOP
!\newpage
! !MODULE: cvmix_background
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  time independent static background mixing coefficients.  It specifies
!  either a scalar, 1D, or 2D field for viscosity and diffusivity. It also
!  calculates the background diffusivity using the Bryan-Lewis method.
!  It then sets the viscosity and diffusivity to the specified value.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

  use cvmix_kinds_and_types, only : cvmix_PI,                  &
                                    cvmix_r8,                  &
                                    cvmix_data_type,           &
                                    cvmix_global_params_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_bkgnd
   public :: cvmix_coeffs_bkgnd
   public :: cvmix_bkgnd_lvary_horizontal
   public :: cvmix_bkgnd_static_diff
   public :: cvmix_bkgnd_static_visc
   public :: cvmix_put_bkgnd
   public :: cvmix_get_bkgnd_real_2D

  interface cvmix_init_bkgnd
    module procedure cvmix_init_bkgnd_scalar
    module procedure cvmix_init_bkgnd_1D
    module procedure cvmix_init_bkgnd_2D
    module procedure cvmix_init_bkgnd_BryanLewis
  end interface cvmix_init_bkgnd

  interface cvmix_put_bkgnd
    module procedure cvmix_put_bkgnd_real
    module procedure cvmix_put_bkgnd_real_1D
    module procedure cvmix_put_bkgnd_real_2D
  end interface cvmix_put_bkgnd

! !PUBLIC TYPES:

  ! cvmix_bkgnd_params_type contains the necessary parameters for background
  ! mixing. Background mixing fields can vary from level to level as well as
  ! over latitude and longitude.
  type, public :: cvmix_bkgnd_params_type
      private
      real(cvmix_r8), allocatable :: static_visc(:,:) ! ncol, nlev+1
      real(cvmix_r8), allocatable :: static_diff(:,:) ! ncol, nlev+1

      ! Note: need to include some logic to avoid excessive memory use
      !       when static_visc and static_diff are constant or 1-D
      logical                     :: lvary_vertical   ! True => second dim not 1
      logical                     :: lvary_horizontal ! True => first dim not 1
  end type cvmix_bkgnd_params_type

!EOP

  type(cvmix_bkgnd_params_type), target :: CVmix_bkgnd_params_saved

contains

!BOP

! !IROUTINE: cvmix_init_bkgnd_scalar
! !INTERFACE:

  subroutine cvmix_init_bkgnd_scalar(bkgnd_diff, bkgnd_visc, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for static background mixing coefficients. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  scalar constants.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), intent(in) :: bkgnd_diff
    real(cvmix_r8), intent(in) :: bkgnd_visc

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), optional, target, intent(inout) :: &
                                              CVmix_bkgnd_params_user
!EOP
!BOC

    type(cvmix_bkgnd_params_type), pointer :: CVmix_bkgnd_params_out

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    if (.not.allocated(CVmix_bkgnd_params_out%static_visc)) then
      CVmix_bkgnd_params_out%lvary_vertical   = .false.
      CVmix_bkgnd_params_out%lvary_horizontal = .false.
      allocate(CVmix_bkgnd_params_out%static_visc(1,1))
      allocate(CVmix_bkgnd_params_out%static_diff(1,1))

      ! Set static_visc and static_diff in background_input_type
      CVmix_bkgnd_params_out%static_visc(1,1) = bkgnd_visc
      CVmix_bkgnd_params_out%static_diff(1,1) = bkgnd_diff
    end if
    ! else error out... can't call init twice!

!EOC

  end subroutine cvmix_init_bkgnd_scalar

!BOP

! !IROUTINE: cvmix_init_bkgnd_1D
! !INTERFACE:

  subroutine cvmix_init_bkgnd_1D(bkgnd_diff, bkgnd_visc, ncol,                &
                                 CVmix_params_user, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for static background mixing coefficients. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  1D field. If field varies horizontally, need to include ncol!
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(:),   intent(in) :: bkgnd_diff
    real(cvmix_r8), dimension(:),   intent(in) :: bkgnd_visc
    integer, optional,              intent(in) :: ncol
    type(cvmix_global_params_type), optional, target, intent(in) :: &
                                                  CVmix_params_user

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  optional, target, intent(inout) :: &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    ! local vars
    integer :: nlev
    type(cvmix_global_params_type), pointer :: CVmix_params_in
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    nullify(CVmix_params_in)
    if (present(CVmix_params_user)) then
      CVmix_params_in => CVmix_params_user
    endif

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    ! NOTE: need to verify that bkgnd_visc and bkgnd_diff are ncol x 1 or
    !       1 x nlev+1

    if (.not.allocated(CVmix_bkgnd_params_out%static_visc)) then
      if (present(ncol)) then
        CVmix_bkgnd_params_out%lvary_vertical   = .false.
        CVmix_bkgnd_params_out%lvary_horizontal = .true.
        allocate(CVmix_bkgnd_params_out%static_visc(ncol,1))
        allocate(CVmix_bkgnd_params_out%static_diff(ncol,1))

        ! Set static_visc and static_diff in background_input_type
        CVmix_bkgnd_params_out%static_visc(:,1) = bkgnd_visc(:)
        CVmix_bkgnd_params_out%static_diff(:,1) = bkgnd_diff(:)
      else
        nlev = CVmix_params_in%max_nlev
        CVmix_bkgnd_params_out%lvary_vertical   = .true.
        CVmix_bkgnd_params_out%lvary_horizontal = .false.
        allocate(CVmix_bkgnd_params_out%static_visc(1,nlev+1))
        allocate(CVmix_bkgnd_params_out%static_diff(1,nlev+1))

        ! Set static_visc and static_diff in background_input_type
        CVmix_bkgnd_params_out%static_visc(1,:) = bkgnd_visc(:)
        CVmix_bkgnd_params_out%static_diff(1,:) = bkgnd_diff(:)
      end if
    end if
    ! else error out... can't call init twice!

!EOC

  end subroutine cvmix_init_bkgnd_1D

!BOP

! !IROUTINE: cvmix_init_bkgnd_2D
! !INTERFACE:

  subroutine cvmix_init_bkgnd_2D(bkgnd_diff, bkgnd_visc, ncol,                &
                                 CVmix_params_user, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for static background mixing coefficients. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  2D field.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_global_params_type), target, optional, intent(in) :: &
                                                  CVmix_params_user
    real(cvmix_r8), dimension(:,:), intent(in) :: bkgnd_diff
    real(cvmix_r8), dimension(:,:), intent(in) :: bkgnd_visc
    integer,                        intent(in) :: ncol

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  target, optional, intent(inout) :: &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    ! local vars
    integer :: nlev
    type(cvmix_global_params_type), pointer :: CVmix_params_in
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    nullify(CVmix_params_in)
    if (present(CVmix_params_user)) then
      CVmix_params_in => CVmix_params_user
    endif

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    ! NOTE: need to verify that bkgnd_visc and bkgnd_diff are ncol x nlev+1

    if (.not.allocated(CVmix_bkgnd_params_out%static_visc)) then
      CVmix_bkgnd_params_out%lvary_vertical   = .true.
      CVmix_bkgnd_params_out%lvary_horizontal = .true.
      nlev = CVmix_params_in%max_nlev
      allocate(CVmix_bkgnd_params_out%static_visc(ncol,nlev+1))
      allocate(CVmix_bkgnd_params_out%static_diff(ncol,nlev+1))

      ! Set static_visc and static_diff in background_input_type
      CVmix_bkgnd_params_out%static_visc(:,:) = bkgnd_visc(:,:)
      CVmix_bkgnd_params_out%static_diff(:,:) = bkgnd_diff(:,:)
    end if
 
!EOC

  end subroutine cvmix_init_bkgnd_2D

!BOP

! !IROUTINE: cvmix_init_bkgnd_BryanLewis
! !INTERFACE:

  subroutine cvmix_init_bkgnd_BryanLewis(CVmix_vars, bl1, bl2, bl3, bl4,      &
                                         CVmix_params_user,                   &
                                         CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Initialization routine for Bryan-Lewis diffusivity/viscosity calculation.
!  For each column, this routine sets the static viscosity \& diffusivity
!  based on the specified parameters. Note that the units of these parameters
!  must be consistent with the units of viscosity and diffusivity -- either
!  cgs or mks, but do not mix and match!
!  \\
!  \\
!  The Bryan-Lewis parameterization is based on the following:
!  \begin{eqnarray*}
!  \kappa_{BL} &=& \textrm{bl1} + \frac{\textrm{bl2}}{\pi}\tan^{-1}\bigg(
!                  \textrm{bl3}(|z|-\textrm{bl4})\bigg)\\
!  \nu_{BL} &=& \textrm{Pr}\cdot\kappa_{BL}
!  \end{eqnarray*}
!  This method is based on the following paper:
!  \begin{quote}
!  \emph{A Water Mass Model of the World Ocean}\\
!  K. Bryan and L. J. Lewis\\
!  Journal of Geophysical Research, vol 84 (1979), pages 2503-2517.
!  \end{quote}
!
!  In that paper, they recommend the parameters 
!  \begin{itemize}
!  \item[] bl1 $= 8 \cdot 10^{-5}$ m$^2/$s
!  \item[] bl2 $= 1.05 \cdot 10^{-4}$ m$^2/$s
!  \item[] bl3 $= 4.5 \cdot 10^{-3}$ m$^{-1}$
!  \item[] bl4 $= 2500$ m
!  \end{itemize}
!  However, more recent usage of their scheme may warrant different settings.
!    
! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    ! Contains depth and nlev
    type(cvmix_data_type),          intent(in) :: CVmix_vars
    ! Contains Prandtl
    type(cvmix_global_params_type), target, optional, intent(in)  :: &
                                                   CVmix_params_user

    ! Units are first column if CVmix_data%depth is m, second if cm
    real(cvmix_r8), intent(in) :: bl1,     &! m^2/s or cm^2/s
                                  bl2,     &! m^2/s or cm^2/s
                                  bl3,     &! 1/m   or 1/cm
                                  bl4       ! m     or cm

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type),  target, optional, intent(inout) :: &
                                               CVmix_bkgnd_params_user
!EOP
!BOC

    ! Pointers to parameter data type
    type(cvmix_global_params_type), pointer :: CVmix_params_in
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_out

    ! Local index
    integer :: nlev  ! max number of levels

    ! Local copies to make code easier to read
    real(cvmix_r8), dimension(:), allocatable :: visc, diff, zw

    nullify(CVmix_params_in)
    if (present(CVmix_params_user)) then
      CVmix_params_in => CVmix_params_user
    endif

    CVmix_bkgnd_params_out => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_out => CVmix_bkgnd_params_user
    end if

    nlev = CVmix_params_in%max_nlev
    allocate(zw(nlev+1))
    allocate(visc(nlev+1))
    allocate(diff(nlev+1))
    
    ! Set static_visc and static_diff in background_input_type
    zw   = abs(CVmix_vars%zw_iface)
    diff = bl1 + (bl2/cvmix_PI)*atan(bl3*(zw-bl4))
    visc = CVmix_params_in%prandtl*diff

    call cvmix_put_bkgnd(CVmix_bkgnd_params_out, "static_diff", diff, nlev=nlev)
    call cvmix_put_bkgnd(CVmix_bkgnd_params_out, "static_visc", visc, nlev=nlev)
    deallocate(zw, visc, diff)

!EOC

  end subroutine cvmix_init_bkgnd_BryanLewis

!BOP

! !IROUTINE: cvmix_coeffs_bkgnd
! !INTERFACE:

  subroutine cvmix_coeffs_bkgnd(CVmix_vars, colid, CVmix_bkgnd_params_user)

! !DESCRIPTION:
!  Computes vertical tracer and velocity mixing coefficients for static
!  background mixing. This routine simply copies viscosity / diffusivity
!  values from CVmix\_bkgnd\_params to CVmix\_vars.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    ! Need to know column for pulling data from static_visc and _diff
    integer,                               optional, intent(in) :: colid
    type(cvmix_bkgnd_params_type), target, optional, intent(in) ::            &
                                           CVmix_bkgnd_params_user

! !INPUT/OUTPUT PARAMETERS:

    type(cvmix_data_type), intent(inout) :: CVmix_vars
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: nlev    ! Number of vertical levels
    type(cvmix_bkgnd_params_type),  pointer :: CVmix_bkgnd_params_in
       
    CVmix_bkgnd_params_in => CVmix_bkgnd_params_saved
    if (present(CVmix_bkgnd_params_user)) then
      CVmix_bkgnd_params_in => CVmix_bkgnd_params_user
    end if

    if (CVmix_bkgnd_params_in%lvary_horizontal.and.(.not.present(colid))) then
      print*, "ERROR: background parameters vary in horizontal so you must ", &
              "provide colid to cvmix_coeffs_bkgnd!"
      stop 1
    end if

    nlev = CVmix_vars%nlev
    if (CVmix_bkgnd_params_in%lvary_horizontal) then
      if (CVmix_bkgnd_params_in%lvary_vertical) then
        CVmix_vars%visc_iface(:)   =                            &
                  CVmix_bkgnd_params_in%static_visc(colid,1:nlev+1)
        CVmix_vars%diff_iface(:,1) =                            &
                  CVmix_bkgnd_params_in%static_diff(colid,1:nlev+1)
      else
        CVmix_vars%visc_iface(:)   =                            &
                  CVmix_bkgnd_params_in%static_visc(colid,1)
        CVmix_vars%diff_iface(:,1) =                            &
                  CVmix_bkgnd_params_in%static_diff(colid,1)
      end if
    else
      if (CVmix_bkgnd_params_in%lvary_vertical) then
        CVmix_vars%visc_iface(:)   =                          &
                  CVmix_bkgnd_params_in%static_visc(1,1:nlev+1)
        CVmix_vars%diff_iface(:,1) =                          &
                  CVmix_bkgnd_params_in%static_diff(1,1:nlev+1)
      else
        CVmix_vars%visc_iface(:)   =                          &
                  CVmix_bkgnd_params_in%static_visc(1,1)
        CVmix_vars%diff_iface(:,1) =                          &
                  CVmix_bkgnd_params_in%static_diff(1,1)
      end if
    end if

!EOC

  end subroutine cvmix_coeffs_bkgnd

!BOP

! !IROUTINE: cvmix_bkgnd_lvary_horizontal
! !INTERFACE:

  function cvmix_bkgnd_lvary_horizontal(CVmix_bkgnd_params_test)

! !DESCRIPTION:
!  Returns whether the background viscosity and diffusivity are
!  varying with horizontal position.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(in) :: CVmix_bkgnd_params_test

! !OUTPUT PARAMETERS:
  logical :: cvmix_bkgnd_lvary_horizontal
!EOP
!BOC

    cvmix_bkgnd_lvary_horizontal = CVmix_bkgnd_params_test%lvary_horizontal

!EOC

  end function cvmix_bkgnd_lvary_horizontal

!BOP

! !IROUTINE: cvmix_bkgnd_static_diff
! !INTERFACE:

  function cvmix_bkgnd_static_diff(CVmix_bkgnd_params_user,kw,colid)

! !DESCRIPTION:
!  Obtain the background diffusivity value at a position in a water column.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(in) :: CVmix_bkgnd_params_user
    integer, optional, intent(in) :: kw, colid

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_bkgnd_static_diff
!EOP
!BOC

    ! Error checks
    if (CVmix_bkgnd_params_user%lvary_horizontal) then
      if (.not.present(colid)) then
        print*, "ERROR: need to pass colid when static_diff varies across", &
                " columns."
        stop 1
      end if
    end if

    if (CVmix_bkgnd_params_user%lvary_vertical) then
      if (.not.present(kw)) then
        print*, "ERROR: need to pass kw (level id) when static_diff varies", &
                "across levels columns."
        stop 1
      end if
    end if

    if (CVmix_bkgnd_params_user%lvary_horizontal) then
      if (CVmix_bkgnd_params_user%lvary_vertical) then
        cvmix_bkgnd_static_diff = CVmix_bkgnd_params_user%static_diff(colid, kw)
      else
        cvmix_bkgnd_static_diff = CVmix_bkgnd_params_user%static_diff(colid, 1)
      end if
    else
      if (CVmix_bkgnd_params_user%lvary_vertical) then
        cvmix_bkgnd_static_diff = CVmix_bkgnd_params_user%static_diff(1, kw)
      else
        cvmix_bkgnd_static_diff = CVmix_bkgnd_params_user%static_diff(1, 1)
      end if
    end if

!EOC

  end function cvmix_bkgnd_static_diff

!BOP

! !IROUTINE: cvmix_bkgnd_static_visc
! !INTERFACE:

  function cvmix_bkgnd_static_visc(CVmix_bkgnd_params_user,kw,colid)

! !DESCRIPTION:
!  Obtain the background viscosity value at a position in a water column.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(in) :: CVmix_bkgnd_params_user
    integer, optional, intent(in) :: kw, colid

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_bkgnd_static_visc
!EOP
!BOC

    ! Error checks
    if (CVmix_bkgnd_params_user%lvary_horizontal) then
      if (.not.present(colid)) then
        print*, "ERROR: need to pass colid when static_visc varies across", &
                " columns."
        stop 1
      end if
    end if

    if (CVmix_bkgnd_params_user%lvary_vertical) then
      if (.not.present(kw)) then
        print*, "ERROR: need to pass kw (level id) when static_visc varies", &
                "across levels columns."
        stop 1
      end if
    end if

    if (CVmix_bkgnd_params_user%lvary_horizontal) then
      if (CVmix_bkgnd_params_user%lvary_vertical) then
        cvmix_bkgnd_static_visc = CVmix_bkgnd_params_user%static_visc(colid, kw)
      else
        cvmix_bkgnd_static_visc = CVmix_bkgnd_params_user%static_visc(colid, 1)
      end if
    else
      if (CVmix_bkgnd_params_user%lvary_vertical) then
        cvmix_bkgnd_static_visc = CVmix_bkgnd_params_user%static_visc(1, kw)
      else
        cvmix_bkgnd_static_visc = CVmix_bkgnd_params_user%static_visc(1, 1)
      end if
    end if

!EOC

  end function cvmix_bkgnd_static_visc

!BOP

! !IROUTINE: cvmix_put_bkgnd_real
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real(CVmix_bkgnd_params_put, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(inout) :: CVmix_bkgnd_params_put
!EOP
!BOC

    select case (trim(varname))
      case ('static_visc')
        if (.not.allocated(CVmix_bkgnd_params_put%static_visc)) then
          allocate(CVmix_bkgnd_params_put%static_visc(1,1))
          CVmix_bkgnd_params_put%lvary_horizontal=.false.
          CVmix_bkgnd_params_put%lvary_vertical=.false.
        else
          print*, "WARNING: overwriting static_visc in cvmix_bkgnd_params_type."
        end if
        CVmix_bkgnd_params_put%static_visc(:,:) = val

      case ('static_diff')
        if (.not.allocated(CVmix_bkgnd_params_put%static_diff)) then
          allocate(CVmix_bkgnd_params_put%static_diff(1,1))
          CVmix_bkgnd_params_put%lvary_horizontal=.false.
          CVmix_bkgnd_params_put%lvary_vertical=.false.
        else
          print*, "WARNING: overwriting static_diff in cvmix_bkgnd_params_type."
        end if
        CVmix_bkgnd_params_put%static_diff(:,:) = val

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_bkgnd_real

!BOP

! !IROUTINE: cvmix_put_bkgnd_real_1D
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real_1D(CVmix_bkgnd_params_put, varname, val, &
                                    ncol, nlev)

! !DESCRIPTION:
!  Write an array of real values into a cvmix\_bkgnd\_params\_type variable.
!  You must use \verb|opt='horiz'| to specify that the field varies in the
!  horizontal direction, otherwise it is assumed to vary in the vertical.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    real(cvmix_r8), dimension(:), intent(in) :: val
    integer, optional,            intent(in) :: ncol, nlev

! !OUTPUT PARAMETERS:
    type (cvmix_bkgnd_params_type), intent(inout) :: CVmix_bkgnd_params_put
!EOP
!BOC

    ! Local vars
    integer, dimension(2) :: dims
    integer               :: data_dims
    logical               :: lvary_horizontal

    ! Error checking to make sure dimension is specified
    if ((.not.present(ncol)).and.(.not.present(nlev))) then
      print*, "ERROR: when putting 1D data in cvmix_bkgnd_params_type ", &
              "you must specify nlev or ncol!"
      stop 1
    end if

    if ((present(ncol)).and.(present(nlev))) then
      print*, "ERROR: when putting 1D data in cvmix_bkgnd_params_type ", &
              "you can not specify both nlev or ncol!"
      stop 1
    end if

    data_dims = size(val)
    if (present(ncol)) then
      if (data_dims.gt.ncol) then
        print*, "ERROR: data array is bigger than number of columns specified."
        stop 1
      end if
      lvary_horizontal=.true.
      dims(1) = ncol
      dims(2) = 1
    else
      if (data_dims.gt.nlev+1) then
        print*, "ERROR: data array is bigger than number of levels specified."
        stop 1
      end if
      lvary_horizontal=.false.
      dims(1) = 1
      dims(2) = nlev+1
    end if

    select case (trim(varname))
      case ('static_visc')
        if (.not.allocated(CVmix_bkgnd_params_put%static_visc)) then
          allocate(CVmix_bkgnd_params_put%static_visc(dims(1),dims(2)))
          CVmix_bkgnd_params_put%lvary_horizontal = lvary_horizontal
          CVmix_bkgnd_params_put%lvary_vertical = .not.lvary_horizontal
        else
          print*, "WARNING: overwriting static_visc in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params_put%static_visc).ne.dims)) then
          print*, "ERROR: dimensions of static_visc do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        if (lvary_horizontal) then
          CVmix_bkgnd_params_put%static_visc(:,1)           = 0_cvmix_r8
          CVmix_bkgnd_params_put%static_visc(1:data_dims,1) = val
        else
          CVmix_bkgnd_params_put%static_visc(1,:)           = 0_cvmix_r8
          CVmix_bkgnd_params_put%static_visc(1,1:data_dims) = val
        end if

      case ('static_diff')
        if (.not.allocated(CVmix_bkgnd_params_put%static_diff)) then
          allocate(CVmix_bkgnd_params_put%static_diff(dims(1),dims(2)))
          CVmix_bkgnd_params_put%lvary_horizontal = lvary_horizontal
          CVmix_bkgnd_params_put%lvary_vertical = .not.lvary_horizontal
        else
          print*, "WARNING: overwriting static_diff in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params_put%static_diff).ne.dims)) then
          print*, "ERROR: dimensions of static_diff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        if (lvary_horizontal) then
          CVmix_bkgnd_params_put%static_diff(:,1)           = 0_cvmix_r8
          CVmix_bkgnd_params_put%static_diff(1:data_dims,1) = val
        else
          CVmix_bkgnd_params_put%static_diff(1,:)           = 0_cvmix_r8
          CVmix_bkgnd_params_put%static_diff(1,1:data_dims) = val
        end if

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_bkgnd_real_1D

!BOP

! !IROUTINE: cvmix_put_bkgnd_real_2D
! !INTERFACE:

  subroutine cvmix_put_bkgnd_real_2D(CVmix_bkgnd_params_put, varname, val, &
                                    ncol, nlev)

! !DESCRIPTION:
!  Write a 2D array of real values into a cvmix\_bkgnd\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    character(len=*),               intent(in) :: varname
    real(cvmix_r8), dimension(:,:), intent(in) :: val
    integer,                        intent(in) :: ncol, nlev

! !OUTPUT PARAMETERS:
    type (cvmix_bkgnd_params_type), intent(out) :: CVmix_bkgnd_params_put
!EOP
!BOC

    ! Local vars
    integer, dimension(2) :: dims, data_dims

    dims      = (/ncol, nlev+1/)
    data_dims = shape(val)

    if (any(data_dims.gt.dims)) then
      print*, "ERROR: data being put in cvmix_bkgnd_params_type is larger ", &
              "than (ncol, nlev+1)"
      stop 1
    end if

    select case (trim(varname))
      case ('static_visc')
        if (.not.allocated(CVmix_bkgnd_params_put%static_visc)) then
          allocate(CVmix_bkgnd_params_put%static_visc(dims(1),dims(2)))
          CVmix_bkgnd_params_put%lvary_horizontal=.true.
          CVmix_bkgnd_params_put%lvary_vertical=.true.
        else
          print*, "WARNING: overwriting static_visc in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params_put%static_visc).ne.dims)) then
          print*, "ERROR: dimensions of static_visc do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        CVmix_bkgnd_params_put%static_visc = 0.0_cvmix_r8
        CVmix_bkgnd_params_put%static_visc(1:data_dims(1), 1:data_dims(2)) = val

      case ('static_diff')
        if (.not.allocated(CVmix_bkgnd_params_put%static_diff)) then
          allocate(CVmix_bkgnd_params_put%static_diff(dims(1),dims(2)))
          CVmix_bkgnd_params_put%lvary_horizontal=.true.
          CVmix_bkgnd_params_put%lvary_vertical=.true.
        else
          print*, "WARNING: overwriting static_diff in cvmix_bkgnd_params_type."
        end if
        if (any(shape(CVmix_bkgnd_params_put%static_diff).ne.dims)) then
          print*, "ERROR: dimensions of static_diff do not match what was ", &
                  "sent to cvmix_put"
          stop 1
        end if
        CVmix_bkgnd_params_put%static_diff = 0.0_cvmix_r8
        CVmix_bkgnd_params_put%static_diff(1:data_dims(1), 1:data_dims(2)) = val

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end subroutine cvmix_put_bkgnd_real_2D

!BOP

! !IROUTINE: cvmix_get_bkgnd_real_2D
! !INTERFACE:

  function cvmix_get_bkgnd_real_2D(CVmix_bkgnd_params_get, varname)

! !DESCRIPTION:
!  Read the real values of a cvmix\_bkgnd\_params\_type 2D array variable.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type(cvmix_bkgnd_params_type), intent(in) :: CVmix_bkgnd_params_get
    character(len=*),              intent(in) :: varname

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(size(CVmix_bkgnd_params_get%static_visc,1),    &
                              size(CVmix_bkgnd_params_get%static_visc,2)) :: &
                    cvmix_get_bkgnd_real_2D
!EOP
!BOC

    select case (trim(varname))
      case ('static_visc')
        cvmix_get_bkgnd_real_2D = CVmix_bkgnd_params_get%static_visc(:,:)

      case ('static_diff')
        cvmix_get_bkgnd_real_2D = CVmix_bkgnd_params_get%static_diff(:,:)

      case DEFAULT
        print*, "ERROR: ", trim(varname), " not a valid choice!"
        stop 1
      
    end select

!EOC

  end function cvmix_get_bkgnd_real_2D


end module cvmix_background

