module vmix_background

!BOP
!\newpage
! !MODULE: vmix_background
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  background mixing (either specifying a scalar, 1D, or 2D field for
!  viscosity and diffusivity coefficients or calculating these coefficients
!  using the Bryan-Lewis method) and to set the viscosity and diffusivity
!  coefficients to this specified value.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use vmix_kinds_and_types, only : vmix_PI,                  &
                                    vmix_r8,                  &
                                    vmix_data_type,           &
                                    vmix_global_params_type,  &
                                    vmix_bkgnd_params_type
!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
   public :: vmix_init_bkgnd
   public :: vmix_coeffs_bkgnd

  interface vmix_init_bkgnd
    module procedure vmix_init_bkgnd_scalar
    module procedure vmix_init_bkgnd_1D
    module procedure vmix_init_bkgnd_2D
    module procedure vmix_init_bkgnd_BryanLewis
  end interface vmix_init_bkgnd
!EOP
!BOC

contains

!BOP

! !IROUTINE: vmix_init_bkgnd_scalar
! !INTERFACE:

  subroutine vmix_init_bkgnd_scalar(Vmix_bkgnd_params, bkgnd_visc, bkgnd_diff)

! !DESCRIPTION:
!  Initialization routine for background mixing with a static field. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  scalar constants.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    real(vmix_r8), intent(in) :: bkgnd_visc
    real(vmix_r8), intent(in) :: bkgnd_diff

! !OUTPUT PARAMETERS:
    type (vmix_bkgnd_params_type), intent(out) :: Vmix_bkgnd_params
!EOP
!BOC

    if (.not.allocated(Vmix_bkgnd_params%static_visc)) then
      Vmix_bkgnd_params%lvary_vertical   = .false.
      Vmix_bkgnd_params%lvary_horizontal = .false.
      allocate(Vmix_bkgnd_params%static_visc(1,1))
      allocate(Vmix_bkgnd_params%static_diff(1,1))

      ! Set static_visc and static_diff in background_input_type
      Vmix_bkgnd_params%static_visc(1,1) = bkgnd_visc
      Vmix_bkgnd_params%static_diff(1,1) = bkgnd_diff
    end if
    ! else error out... can't call init twice!

!EOC

  end subroutine vmix_init_bkgnd_scalar

!BOP

! !IROUTINE: vmix_init_bkgnd_1D
! !INTERFACE:

  subroutine vmix_init_bkgnd_1D(Vmix_params, Vmix_bkgnd_params, &
                                bkgnd_visc, bkgnd_diff, ncol)

! !DESCRIPTION:
!  Initialization routine for background mixing with a static field. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  1D field. If field varies horizontally, need to include ncol!
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type (vmix_global_params_type), intent(in) :: Vmix_params
    real(vmix_r8), dimension(:),    intent(in) :: bkgnd_visc
    real(vmix_r8), dimension(:),    intent(in) :: bkgnd_diff
    integer, optional,              intent(in) :: ncol

! !OUTPUT PARAMETERS:
    type (vmix_bkgnd_params_type), intent(out) :: Vmix_bkgnd_params
!EOP
!BOC

    ! local vars
    integer :: nlev

    ! NOTE: need to verify that bkgnd_visc and bkgnd_diff are ncol x 1 or
    !       1 x nlev+1

    if (.not.allocated(Vmix_bkgnd_params%static_visc)) then
      if (present(ncol)) then
        Vmix_bkgnd_params%lvary_vertical   = .false.
        Vmix_bkgnd_params%lvary_horizontal = .true.
        allocate(Vmix_bkgnd_params%static_visc(ncol,1))
        allocate(Vmix_bkgnd_params%static_diff(ncol,1))

        ! Set static_visc and static_diff in background_input_type
        Vmix_bkgnd_params%static_visc(:,1) = bkgnd_visc(:)
        Vmix_bkgnd_params%static_diff(:,1) = bkgnd_diff(:)
      else
        nlev = Vmix_params%max_nlev
        Vmix_bkgnd_params%lvary_vertical   = .true.
        Vmix_bkgnd_params%lvary_horizontal = .false.
        allocate(Vmix_bkgnd_params%static_visc(1,nlev+1))
        allocate(Vmix_bkgnd_params%static_diff(1,nlev+1))

        ! Set static_visc and static_diff in background_input_type
        Vmix_bkgnd_params%static_visc(1,:) = bkgnd_visc(:)
        Vmix_bkgnd_params%static_diff(1,:) = bkgnd_diff(:)
      end if
    end if
    ! else error out... can't call init twice!

!EOC

  end subroutine vmix_init_bkgnd_1D

!BOP

! !IROUTINE: vmix_init_bkgnd_2D
! !INTERFACE:

  subroutine vmix_init_bkgnd_2D(Vmix_params, Vmix_bkgnd_params, &
                                bkgnd_visc, bkgnd_diff, ncol)

! !DESCRIPTION:
!  Initialization routine for background mixing with a static field. For each
!  column, this routine sets the static viscosity / diffusivity to the given
!  2D field.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type (vmix_global_params_type), intent(in) :: Vmix_params
    real(vmix_r8), dimension(:,:),  intent(in) :: bkgnd_visc
    real(vmix_r8), dimension(:,:),  intent(in) :: bkgnd_diff
    integer,                        intent(in) :: ncol

! !OUTPUT PARAMETERS:
    type (vmix_bkgnd_params_type), intent(out) :: Vmix_bkgnd_params
!EOP
!BOC

    ! local vars
    integer :: nlev

    ! NOTE: need to verify that bkgnd_visc and bkgnd_diff are ncol x nlev+1

    if (.not.allocated(Vmix_bkgnd_params%static_visc)) then
      Vmix_bkgnd_params%lvary_vertical   = .true.
      Vmix_bkgnd_params%lvary_horizontal = .true.
      nlev = Vmix_params%max_nlev
      allocate(Vmix_bkgnd_params%static_visc(ncol,nlev+1))
      allocate(Vmix_bkgnd_params%static_diff(ncol,nlev+1))

      ! Set static_visc and static_diff in background_input_type
      Vmix_bkgnd_params%static_visc(:,:) = bkgnd_visc(:,:)
      Vmix_bkgnd_params%static_diff(:,:) = bkgnd_diff(:,:)
    end if
    ! else error out... can't call init twice!
 
!EOC

  end subroutine vmix_init_bkgnd_2D

!BOP

! !IROUTINE: vmix_init_bkgnd_BryanLewis
! !INTERFACE:

  subroutine vmix_init_bkgnd_BryanLewis(Vmix_vars, Vmix_params,         &
                                        Vmix_bkgnd_params, colid, ncol, &
                                        bl1, bl2, bl3, bl4)

! !DESCRIPTION:
!  Initialization routine for background mixing with a Bryan-Lewis mixing.
!  For each column, this routine sets the static viscosity / diffusivity
!  based on the specified parameters. Note that the units of these parameters
!  must be consistent with the units of viscosity and diffusivity -- either
!  cgs or mks, but don't mix and match!
!  \\
!  \\
!  The Bryan-Lewis parameterization uses the following:
!  \begin{eqnarray*}
!  \kappa_{BL} &=& \textrm{bl1} + \frac{\textrm{bl2}}{\pi}\tan^{-1}\bigg(
!                  \textrm{bl3}(z-\textrm{bl4})\bigg)\\
!  \nu_{BL} &=& \textrm{Pr}\cdot\kappa_{BL}
!  \end{eqnarray*}
!  This is all based on the following paper:
!  \begin{quote}
!  \emph{A Water Mass Model of the World Ocean}\\
!  K. Bryan and L. J. Lewis\\
!  Journal of Geophysical Research, vol 84 (1979), pages 2503-2517.
!  \end{quote}
!
!  In that paper,
!  \begin{itemize}
!  \item[] bl1 $= 8 \cdot 10^{-5}$ m$^2/$s
!  \item[] bl2 $= 1.05 \cdot 10^{-4}$ m$^2/$s
!  \item[] bl3 $= 4.5 \cdot 10^{-3}$ m$^{-1}$
!  \item[] bl4 $= 2500$ m
!  \end{itemize}

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:
    type (vmix_data_type),          intent(in) :: Vmix_vars   ! Depth, nlev
    type (vmix_global_params_type), intent(in) :: Vmix_params ! Prandtl
    integer,                        intent(in) :: colid, ncol

    ! Units are first column if Vmix_data%depth is m, second if cm
    real(vmix_r8), intent(in)         :: bl1,     &! m^2/s or cm^2/s
                                         bl2,     &! m^2/s or cm^2/s
                                         bl3,     &! 1/m   or 1/cm
                                         bl4       ! m     or cm

! !OUTPUT PARAMETERS:
    type (vmix_bkgnd_params_type), intent(out) :: Vmix_bkgnd_params
!EOP
!BOC

    ! Local index
    integer                :: kw,  &! vertical cell index
                              nlev  ! max number of levels

    ! Local copies to make code easier to read
    real(vmix_r8) :: visc, diff, z

    if (.not.allocated(Vmix_bkgnd_params%static_visc)) then
      Vmix_bkgnd_params%lvary_vertical   = .true.
      Vmix_bkgnd_params%lvary_horizontal = .true.
      nlev = Vmix_params%max_nlev
      allocate(Vmix_bkgnd_params%static_diff(ncol,nlev+1))
      allocate(Vmix_bkgnd_params%static_visc(ncol,nlev+1))
    end if
    Vmix_bkgnd_params%static_diff(colid,:) = 0.0d0
    Vmix_bkgnd_params%static_visc(colid,:) = 0.0d0
    
    ! Set static_visc and static_diff in background_input_type
    do kw = 1, Vmix_vars%nlev+1
       z    = Vmix_vars%z_iface(kw)
       diff = bl1 + (bl2/vmix_PI)*atan(bl3*(z-bl4))
       visc = Vmix_params%prandtl*diff

       Vmix_bkgnd_params%static_diff(colid, kw) = diff
       Vmix_bkgnd_params%static_visc(colid, kw) = visc
    end do

!EOC

  end subroutine vmix_init_bkgnd_BryanLewis

!BOP

! !IROUTINE: vmix_coeffs_bkgnd
! !INTERFACE:

  subroutine vmix_coeffs_bkgnd(Vmix_vars, Vmix_bkgnd_params, colid)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for static mixing. This routine
!  simply copies viscosity / diffusivity values from Vmix\_bkgnd\_params to
!  Vmix\_vars.
!\\
!\\

! !USES:
!  Only those used by entire module. 

! !INPUT PARAMETERS:

    type (vmix_bkgnd_params_type), intent(in) :: Vmix_bkgnd_params
    ! Need to know column for pulling data from static_visc and _diff
    integer,                       intent(in) :: colid

! !INPUT/OUTPUT PARAMETERS:

    type (vmix_data_type), intent(inout) :: Vmix_vars
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: nlev    ! Number of vertical levels
       
    nlev = Vmix_vars%nlev
    if (Vmix_bkgnd_params%lvary_horizontal) then
      if (Vmix_bkgnd_params%lvary_vertical) then
        Vmix_vars%visc_iface(:)   =                            &
                  Vmix_bkgnd_params%static_visc(colid,1:nlev+1)
        Vmix_vars%diff_iface(:,1) =                            &
                  Vmix_bkgnd_params%static_diff(colid,1:nlev+1)
      else
        Vmix_vars%visc_iface(:)   =                            &
                  Vmix_bkgnd_params%static_visc(colid,1)
        Vmix_vars%diff_iface(:,1) =                            &
                  Vmix_bkgnd_params%static_diff(colid,1)
      end if
    else
      if (Vmix_bkgnd_params%lvary_vertical) then
        Vmix_vars%visc_iface(:)   =                          &
                  Vmix_bkgnd_params%static_visc(1,1:nlev+1)
        Vmix_vars%diff_iface(:,1) =                          &
                  Vmix_bkgnd_params%static_diff(1,1:nlev+1)
      else
        Vmix_vars%visc_iface(:)   =                          &
                  Vmix_bkgnd_params%static_visc(1,1)
        Vmix_vars%diff_iface(:,1) =                          &
                  Vmix_bkgnd_params%static_diff(1,1)
      end if
    end if

!EOC

  end subroutine vmix_coeffs_bkgnd

end module vmix_background

