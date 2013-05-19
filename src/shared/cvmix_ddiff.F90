!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_ddiff

!BOP
!\newpage
! !MODULE: cvmix_ddiff
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  double diffusion mixing and to set the diffusivity coefficient
!  accordingly.
!\\
!\\
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !USES:

   use cvmix_kinds_and_types, only : one,                     &
                                     cvmix_r8,                &
                                     cvmix_data_type,         &
                                     cvmix_ddiff_params_type
   use cvmix_put_get, only : cvmix_put
!EOP

   implicit none
   private
   save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

   public :: cvmix_init_ddiff
   public :: cvmix_coeffs_ddiff
!EOP

 contains

!BOP

! !IROUTINE: cvmix_init_ddiff
! !INTERFACE:

  subroutine cvmix_init_ddiff(CVmix_ddiff_params, units, strat_param_max, &
                              kappa_ddiff_t, kappa_ddiff_s, ddiff_exp1, &
                              ddiff_exp2, mol_diff, kappa_ddiff_param1, &
                              kappa_ddiff_param2, kappa_ddiff_param3)

! !DESCRIPTION:
!  Initialization routine for double diffusion mixing. This mixing technique
!  looks for two unstable cases in a column - salty water over fresher
!  water and colder water over warmer water - and computes different
!  diffusivity coefficients in each of these two locations. The parameter
!  \begin{eqnarray*}
!  R_\rho = \frac{\alpha (\partial \Theta / \partial z)}
!                {\beta (\partial S / \partial z)}
!  \end{eqnarray*}
!  to determine as a stratification parameter. If $(\partial S / \partial z)$
!  is positive and $1 < R_\rho < R_\rho^0$ then salt water sits on top
!  of fresh water and the diffusivity is given by
!  \begin{eqnarray*}
!  \kappa = \kappa^0 \left[ 1 - \left(\frac{R_\rho - 1}{R_\rho^0 - 1} \right)^{p_1}\right]^{p_2}
!  \end{eqnarray*}
!  The user must specify which set of units to use, either \verb|'mks'| or \verb|'cgs'|.
!  By default, $R_\rho^0 = 2.55$, but that can be changed by setting 
!  \verb|strat_param_max| in the code. Similarly, by default $p_1 = 1$ 
! (\verb|ddiff_exp1|), $p_2 = 3$ (\verb|ddiff_exp2|), and
!  \begin{eqnarray*}
!  \kappa^0 = \left\{ \begin{array}{r l}
!             7 \cdot 10^{-5}\ \textrm{m}^2\textrm{/s} & \textrm{for temperature}
!             \ (\verb|kappa_ddiff_t|\ \textrm{in this routine})\\
!             10^{-4}\ \textrm{m}^2\textrm{/s} & \textrm{for salinity and other tracers}
!             \ (\verb|kappa_ddiff_s|\ \textrm{in this routine}).
!                     \end{array} \right.
!  \end{eqnarray*}
!  On the other hand, if $(\partial \Theta / \partial z)$ is negative and
!  $0 < R_\rho < 1$ then cold water sits on warm warm water and the
!  diffusivity for temperature is given by
!  \begin{eqnarray*}
!  \kappa = \nu_\textrm{molecular} \cdot 0.909\exp\left\{ 4.6\exp\left[
!           -0.54\left( \frac{1}{R_\rho} - 1 \right) \right] \right\}
!  \end{eqnarray*}
!  where $\nu_\textrm{molecular}$ Is the molecular viscosity of water. By default it
!  is set to $1.5 \cdot 10^{-6}\ \textrm{m}^2\textrm{/s}$, but it can be changed
!  through \verb|mol_diff| in the code. Similarly, 0.909, 4.6, and -0.54 are the
!  default values of \verb|kappa_ddiff_param1|, \verb|kappa_ddiff_param2|, and
!  \verb|kappa_ddiff_param3|, respectively.\\
!\\
!  For salinity and other tracers, $\kappa$ above is multiplied by the factor
!  \begin{eqnarray*}
!  \textrm{factor} = \left\{ \begin{array}{c l}
!                    0.15R_\rho & R_\rho < 0.5\\
!                    1.85R_\rho - 0.85 & 0.5 \le R_\rho < 1\\
!                     \end{array} \right.
!  \end{eqnarray*}
!  $\kappa$ is stored in \verb|CVmix_vars%diff_iface(:,1)|, while the modified value
!  for non-temperature tracers is stored in \verb|CVmix_vars%diff_iface(:,2)|.\\
!\\
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: units ! "mks" or "cgs"
    real(cvmix_r8),   optional, intent(in) :: strat_param_max, &
                                              kappa_ddiff_t, &
                                              kappa_ddiff_s, &
                                              ddiff_exp1, &
                                              ddiff_exp2, &
                                              mol_diff, &
                                              kappa_ddiff_param1, &
                                              kappa_ddiff_param2, &
                                              kappa_ddiff_param3

! !OUTPUT PARAMETERS:
    type(cvmix_ddiff_params_type), intent(inout) :: CVmix_ddiff_params
!EOP
!BOC

    ! Unitless parameters
    if (present(strat_param_max)) then
      call cvmix_put(CVmix_ddiff_params, "strat_param_max", strat_param_max)
    else
      call cvmix_put(CVmix_ddiff_params, "strat_param_max", 2.55_cvmix_r8)
    end if
    if (present(ddiff_exp1)) then
      call cvmix_put(CVmix_ddiff_params, "ddiff_exp1", ddiff_exp1)
    else
      call cvmix_put(CVmix_ddiff_params, "ddiff_exp1", 1.0_cvmix_r8)
    end if
    if (present(ddiff_exp2)) then
      call cvmix_put(CVmix_ddiff_params, "ddiff_exp2", ddiff_exp2)
    else
      call cvmix_put(CVmix_ddiff_params, "ddiff_exp2", 3.0_cvmix_r8)
    end if
    if (present(kappa_ddiff_param1)) then
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_param1", kappa_ddiff_param1)
    else
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_param1", 0.909_cvmix_r8)
    end if
    if (present(kappa_ddiff_param2)) then
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_param2", kappa_ddiff_param2)
    else
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_param2", 4.6_cvmix_r8)
    end if
    if (present(kappa_ddiff_param3)) then
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_param3", kappa_ddiff_param3)
    else
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_param3", -0.54_cvmix_r8)
    end if

    ! Parameters with units
    if (present(kappa_ddiff_t)) then
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_t", kappa_ddiff_t)
    end if
    if (present(kappa_ddiff_s)) then
      call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_s", kappa_ddiff_s)
    end if
    if (present(mol_diff)) then
      call cvmix_put(CVmix_ddiff_params, "mol_diff", mol_diff)
    end if

    select case (trim(units))
      case ('mks')
        if (.not.present(kappa_ddiff_t)) then
          call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_t", 7d-5)
        end if
        if (.not.present(kappa_ddiff_s)) then
          call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_s", 1d-4)
        end if
        if (.not.present(mol_diff)) then
          call cvmix_put(CVmix_ddiff_params, "mol_diff", 1.5d-6)
        end if
      case ('cgs')
        if (.not.present(kappa_ddiff_t)) then
          call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_t", 7d-1)
        end if
        if (.not.present(kappa_ddiff_s)) then
          call cvmix_put(CVmix_ddiff_params, "kappa_ddiff_s", 1.0_cvmix_r8)
        end if
        if (.not.present(mol_diff)) then
          call cvmix_put(CVmix_ddiff_params, "mol_diff", 1.5d-2)
        end if

      case DEFAULT
        print*, "ERROR: ", trim(units), " is not a valid choice for double ", &
                "diffusion mixing. Only 'mks' and 'cgs' are supported."
        stop 1

    end select

!EOC

  end subroutine cvmix_init_ddiff

!***********************************************************************
!BOP
! !IROUTINE: cvmix_coeffs_ddiff
! !INTERFACE:

  subroutine cvmix_coeffs_ddiff(CVmix_vars, CVmix_ddiff_params)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for the double diffusion mixing
!  parameterizatiion.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_ddiff_params_type), intent(in) :: CVmix_ddiff_params

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

! !LOCAL VARIABLES:
    integer :: k ! column index
    real(cvmix_r8) :: ddiff, Rrho

!EOP
!BOC

    ! Determine coefficients based on units requested
    CVmix_vars%diff_iface = 0_cvmix_r8
    do k = 1, CVmix_vars%nlev
      if ((CVmix_vars%strat_param_num(k).gt.CVmix_vars%strat_param_denom(k)).and.&
          (CVmix_vars%strat_param_denom(k).gt.0)) then
        ! Rrho > 1 and dS/dz < 0 => Salt fingering
        Rrho = CVmix_vars%strat_param_num(k) / CVmix_vars%strat_param_denom(k)
        if (Rrho.lt.CVmix_ddiff_params%strat_param_max) then
          ddiff = (one-((Rrho-one)/(CVmix_ddiff_params%strat_param_max-one))** &
                  CVmix_ddiff_params%ddiff_exp1)**CVmix_ddiff_params%ddiff_exp2
          CVmix_vars%diff_iface(k,1) = CVmix_ddiff_params%kappa_ddiff_t*ddiff
          CVmix_vars%diff_iface(k,2) = CVmix_ddiff_params%kappa_ddiff_s*ddiff
        end if
      end if
      if ((CVmix_vars%strat_param_num(k).gt.CVmix_vars%strat_param_denom(k)).and.&
          (CVmix_vars%strat_param_num(k).lt.0)) then
        ! Rrho < 1 and dT/dz > 0 => Diffusive convection
        Rrho = CVmix_vars%strat_param_num(k) / CVmix_vars%strat_param_denom(k)
        ddiff = CVmix_ddiff_params%mol_diff*CVmix_ddiff_params%kappa_ddiff_param1*&
                exp(CVmix_ddiff_params%kappa_ddiff_param2*exp(&
                CVmix_ddiff_params%kappa_ddiff_param3*(one/Rrho-one)))
        CVmix_vars%diff_iface(k,1) = ddiff
        if (Rrho.lt.0.5_cvmix_r8) then
          CVmix_vars%diff_iface(k,2) = 0.15_cvmix_r8*Rrho*ddiff
        else
          CVmix_vars%diff_iface(k,2) = (1.85_cvmix_r8*Rrho-0.85_cvmix_r8)*ddiff
        end if
      end if
    end do
    CVmix_vars%diff_iface(CVmix_vars%nlev+1,:) = 0.0_cvmix_r8

!EOC
  end subroutine cvmix_coeffs_ddiff

end module cvmix_ddiff
