   subroutine cvmix_test(nlev, ncol, ocn_depth, col1_vdc1, col1_vdc2, col1_linv, col1_dpth, col2_vdc1, col2_vdc2, col2_linv, col2_dpth, iface_depth, viscosity, diffusivity)
     use cvmix_kinds_and_types
     use cvmix_background
     use cvmix_put_get

     Implicit None
   
     type(cvmix_data_type)         , dimension(2) :: CVmix_vars
     type(cvmix_global_params_type)               :: CVmix_params
     type(cvmix_bkgnd_params_type) , dimension(2) :: CVmix_BL_params

     ! "Namelist" variables
     ! 1) General mixing parameters
     integer,        intent(in)         :: ncol      ! number of columns
     integer,        intent(in)         :: nlev      ! number of levels for column
     real(cvmix_r8), intent(in)         :: ocn_depth ! Depth of ocn
     ! 2) Bryan-Lewis mixing parameters for column 1
     real(cvmix_r8), intent(in)         :: col1_vdc1, col1_vdc2, col1_linv, col1_dpth
     ! 3) Bryan-Lewis mixing parameters for column 2
     real(cvmix_r8), intent(in)         :: col2_vdc1, col2_vdc2, col2_linv, col2_dpth

   
     ! Will use 2 columns, viscosity will be 2 x nlev+1 and diffusivity will 
     ! be 2 x nlev+1 x 1 (diffusivity is 2D in column)
     ! iface_depth is the depth of each interface;  same in both columns
     real(cvmix_r8), dimension(nlev+1),     intent(in), target :: iface_depth
     real(cvmix_r8), dimension(ncol, nlev+1),  intent(out), target :: viscosity
     real(cvmix_r8), dimension(ncol, nlev+1, 1), intent(out), target :: diffusivity
   
     ! Global parameter
     ! array indices
     integer :: icol,kw
   
     ! Initialization for CVMix data types
     call cvmix_put(CVmix_params,  'max_nlev', nlev)
     call cvmix_put(CVmix_params,  'prandtl',  0.0_cvmix_r8)
     do icol=1,ncol
       call cvmix_put(CVmix_vars(icol), 'nlev',     nlev)
       ! Point CVmix_vars values to memory allocated above
       CVmix_vars(icol)%visc_iface => viscosity(icol,:)
       CVmix_vars(icol)%diff_iface => diffusivity(icol,:,:)
       CVmix_vars(icol)%z_iface => iface_depth
     end do
   
     ! Read / set B-L parameters for column 1
     call cvmix_init_bkgnd(CVmix_vars(1), CVmix_params, CVmix_BL_params(1), &
                           col1_vdc1, col1_vdc2, col1_linv, col1_dpth)
     call cvmix_coeffs_bkgnd(CVmix_vars(1), CVmix_BL_params(1), 1)
     
     ! Read / set B-L parameters for column 2
     call cvmix_init_bkgnd(CVmix_vars(2), CVmix_params, CVmix_BL_params(2), &
                           col2_vdc1, col2_vdc2, col2_linv, col2_dpth)
     call cvmix_coeffs_bkgnd(CVmix_vars(2), CVmix_BL_params(2), 1)
     
  end subroutine cvmix_test
