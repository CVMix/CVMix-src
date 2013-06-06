!BOP
!\newpage
! !ROUTINE: cvmix_tidal_driver

! !DESCRIPTION: A routine to test the Simmons implementation of tidal mixing.
!\\
!\\

! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$URL$

! !INTERFACE:

Subroutine cvmix_tidal_driver()

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                 &
                                    cvmix_strlen,             &
                                    cvmix_data_type,          &
                                    cvmix_global_params_type, &
                                    cvmix_tidal_params_type
  use cvmix_tidal,           only : cvmix_init_tidal,         &
                                    cvmix_coeffs_tidal
  use cvmix_put_get,         only : cvmix_put
  use cvmix_io,              only : cvmix_io_open,            &
                                    cvmix_input_read,         &
                                    cvmix_output_write,       &
                                    cvmix_output_write_att,   &
                                    cvmix_io_close

  Implicit None

!EOP
!BOC

  ! CVMix datatypes
  type(cvmix_data_type), dimension(:,:), allocatable :: CVmix_vars
  type(cvmix_global_params_type) :: CVmix_params
  type(cvmix_tidal_params_type)  :: CVmix_Simmons_params

  real(cvmix_r8), dimension(:,:,:,:), allocatable, target :: diffusivity

  ! file index
  integer :: fid

  ! Namelist variables
  character(len=cvmix_strlen) :: grid_file, physics_file, energy_flux_file, &
                                 energy_flux_var
  integer :: lon_out, lat_out

  ! Local variables
  real(cvmix_r8), dimension(:,:,:), allocatable :: buoy
  real(cvmix_r8), dimension(:,:),   allocatable :: ocn_depth, energy_flux, &
                                                   lat, lon
  integer,        dimension(:,:),   allocatable :: ocn_levels
  real(cvmix_r8), dimension(:),     allocatable :: depth_iface, depth
  real(cvmix_r8)                                :: FillVal, this_lon, this_lat
  character(len=cvmix_strlen) :: lonstr, latstr
  integer :: i, j, k, nlon, nlat, nlev, max_nlev

  ! Namelists that may be read in, depending on desired mixing scheme
  namelist/Simmons_nml/grid_file, physics_file, energy_flux_file,             &
                       energy_flux_var, lon_out, lat_out

  ! Hardcode in file dimensions (to do: read in from netCDF)
  nlon = 320
  nlat = 384
  max_nlev = 60
  ! Allocate memory for CVmix columns
  allocate(CVmix_vars(nlon, nlat))
  allocate(lat(nlon, nlat),lon(nlon, nlat))

  ! Read namelist variables 
  grid_file = "none"
  physics_file = "none"
  energy_flux_file = "none"
  energy_flux_var = "none"
  lon_out = 35
  lat_out = 345
  read(*, nml=Simmons_nml)

  ! Allocate memory for energy flux, ocean depth, number of ocean levels,
  ! depth of each level / interface, and buoyancy frequency
  allocate(energy_flux(nlon, nlat), ocn_depth(nlon, nlat))
  allocate(buoy(nlon, nlat,max_nlev+1))
  allocate(ocn_levels(nlon, nlat))
  allocate(depth(max_nlev), depth_iface(max_nlev+1))
  ! Set buoyancy frequency = 0 at top interface (POP doesn't store these
  ! zeroes and the input data set is coming from POP output)
  buoy(:,:,1) = 0.0_cvmix_r8

  ! Allocate memory to store diffusivity values
  allocate(diffusivity(nlon, nlat, max_nlev+1,1))
  ! Set diffusivity to _FillValue
  FillVal     = 100000.0_cvmix_r8
  diffusivity = FillVal

  ! Read in global data from grid file, physics file, and energy flux file
  call cvmix_io_open(fid, trim(grid_file), 'nc', read_only=.true.)
  call cvmix_input_read(fid, 'lon', lon)
  call cvmix_input_read(fid, 'lat', lat)
  call cvmix_input_read(fid, 'zw', depth_iface)
  call cvmix_input_read(fid, 'H', ocn_depth)
  call cvmix_input_read(fid, 'H_index', ocn_levels)
  call cvmix_io_close(fid)
  call cvmix_io_open(fid, trim(physics_file), 'nc', read_only=.true.)
  call cvmix_input_read(fid, 'Nsqr', buoy(:,:,2:max_nlev+1))
  call cvmix_io_close(fid)
  call cvmix_io_open(fid, trim(energy_flux_file), 'nc', read_only=.true.)
  call cvmix_input_read(fid, trim(energy_flux_var), energy_flux)
  call cvmix_io_close(fid)

  ! Compute center of each layer (maybe this should be stored in grid file?)
  do k=1, max_nlev
    depth(k) = 0.5_cvmix_r8*(depth_iface(k)+depth_iface(k+1))
  end do

  ! Initialize tidal mixing parameters
  call cvmix_init_tidal(CVmix_Simmons_params, 'Simmons', 'mks', &
                        local_mixing_frac=0.33_cvmix_r8,        &
                        max_coefficient=0.01_cvmix_r8)
  print*, "Namelist variables"
  print*, "------------------"
  print*, "mix_scheme = ", trim(CVmix_Simmons_params%mix_scheme)
  print*, "efficiency = ", CVmix_Simmons_params%efficiency
  print*, "vertical_decay_scale = ", CVmix_Simmons_params%vertical_decay_scale
  print*, "max_coefficient = ", CVmix_Simmons_params%max_coefficient
  print*, "local_mixing_frac = ", CVmix_Simmons_params%local_mixing_frac
  print*, "depth_cutoff = ", CVmix_Simmons_params%depth_cutoff

  ! For starters, using column from 353.9634 E, 58.84838 N)
  ! That's i=35, j=345 (compare result to KVMIX(0, :, 344, 34) in NCL)
  do i=1,nlon
    do j=1,nlat
      nlev = ocn_levels(i,j)

      ! Initialization for CVMix data types
      call cvmix_put(CVmix_vars(i,j), 'nlev', nlev)
      if (nlev.gt.0) then
        call cvmix_put(CVmix_vars(i,j),  'surf_hgt',          0.0_cvmix_r8)
        call cvmix_put(CVmix_vars(i,j),  'zw_iface', depth_iface(1:nlev+1))
        call cvmix_put(CVmix_vars(i,j),        'zw',         depth(1:nlev))
        call cvmix_put(CVmix_vars(i,j),      'buoy',    buoy(i,j,1:nlev+1))
        call cvmix_put(CVmix_vars(i,j), 'ocn_depth',        ocn_depth(i,j))

        call cvmix_put(CVmix_params, 'max_nlev',        max_nlev)
        call cvmix_put(CVmix_params,   'fw_rho', 1000.0_cvmix_r8)
        ! Point CVmix_vars values to memory allocated above
        CVmix_vars(i,j)%diff_iface => diffusivity(i,j,1:nlev+1,:)

        call cvmix_coeffs_tidal(CVmix_vars(i,j), CVmix_Simmons_params, &
                                CVmix_params, energy_flux(i,j))

      end if

      ! Output
      if ((i.eq.lon_out).and.(j.eq.lat_out)) then
        if (nlev.gt.0) then
          this_lon = lon(lon_out, lat_out)
          ! Need this_lon between -180 and 180
          do while(this_lon.lt.-180.0_cvmix_r8)
            this_lon = this_lon + 360.0_cvmix_r8
          end do
          do while(this_lon.gt.180.0_cvmix_r8)
            this_lon = this_lon - 360.0_cvmix_r8
          end do
          this_lat = lat(lon_out, lat_out)
          call cvmix_io_open(fid, "single_col.nc", "nc")
          call cvmix_output_write(fid, CVmix_vars(i,j), (/"depth", "diff "/))
          if (this_lon.ge.0) then
            write(lonstr,"(F6.2,1X,A)") this_lon, "E"
          else
            write(lonstr,"(F6.2,1X,A)") abs(this_lon), "W"
          end if
          if (this_lat.ge.0) then
            write(latstr,"(F6.2,1X,A)") this_lat, "N"
          else
            write(latstr,"(F6.2,1X,A)") abs(this_lat), "S"
          end if
          ! Global Attributes
          call cvmix_output_write_att(fid, "column_lon", lonstr)
          call cvmix_output_write_att(fid, "column_lat", latstr)

          ! Variable Attributes
          call cvmix_output_write_att(fid, "long_name", "tracer diffusivity", &
                                      var_name="diff")
          call cvmix_output_write_att(fid, "units", "m^2/s", var_name="diff")
          call cvmix_output_write_att(fid, "long_name", "depth to interface", &
                                      var_name="depth")
          call cvmix_output_write_att(fid, "positive", "up", var_name="depth")
          call cvmix_output_write_att(fid, "units", "m", var_name="depth")
          call cvmix_io_close(fid)
        else
          print*, "ERROR: column requested for output is a land cell."
          stop 1
        end if
      end if

    end do
  end do

  ! Write diffusivity field to netcdf
  call cvmix_io_open(fid, "diff.nc", "nc")
  call cvmix_output_write(fid, "diff", (/"nlon  ", "nlat  ", "niface"/),      &
                          diffusivity(:,:,:,1), FillVal=FillVal)
  call cvmix_output_write_att(fid, "long_name", "tracer diffusivity",         &
                              var_name="diff")
  call cvmix_output_write_att(fid, "units", "m^2/s", var_name="diff")
  call cvmix_io_close(fid)

  ! memory cleanup
  deallocate(CVmix_vars)
  deallocate(energy_flux, ocn_depth)
  deallocate(buoy)
  deallocate(ocn_levels)
  deallocate(depth, depth_iface)

!EOC

End Subroutine cvmix_tidal_driver
