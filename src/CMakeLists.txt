# CVMix library
add_library(cvmix_objects OBJECT .)
set_property( TARGET cvmix_objects PROPERTY FOLDER cvmixlib )
target_sources( cvmix_objects PRIVATE
   shared/cvmix_kinds_and_types.F90
   shared/cvmix_background.F90
   shared/cvmix_convection.F90
   shared/cvmix_ddiff.F90
   shared/cvmix_kpp.F90
   shared/cvmix_math.F90
   shared/cvmix_put_get.F90
   shared/cvmix_shear.F90
   shared/cvmix_tidal.F90
   shared/cvmix_utils.F90
)

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cvmix_version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/cvmix_version.F90)

# CVMix driver dependencies
if(CVMIX_BUILD_DRIVER)

   add_library(cvmix_io STATIC .)
   set_property( TARGET cvmix_io PROPERTY FOLDER driver )
   target_sources( cvmix_io PRIVATE
      cvmix_io.F90
   )
#   target_link_libraries(cvmix_io PRIVATE cvmix_static PUBLIC netcdff )
   target_link_libraries(cvmix_io PRIVATE cvmix_static PUBLIC $ENV{NetCDF_LIBRARIES} )

   add_library(cvmix_drivers STATIC .)
   set_property( TARGET cvmix_drivers PROPERTY FOLDER driver )
   target_sources( cvmix_drivers PRIVATE
      drivers/cvmix_bgrnd_BL.F90
      drivers/cvmix_ddiff_drv.F90
      drivers/cvmix_kpp_drv.F90
      drivers/cvmix_shear_drv.F90
      drivers/cvmix_tidal_Simmons.F90
   )
   target_link_libraries( cvmix_drivers PRIVATE cvmix_io )

endif()
