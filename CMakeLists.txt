cmake_minimum_required( VERSION 3.9 )

project( cvmix VERSION 4.0.1 LANGUAGES Fortran )

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Use standard GNU installation directories
if ( NOT WIN32 )
   include( GNUInstallDirs )
endif()

# Configuration options
option( CVMIX_USE_NetCDF "Enable NetCDF format" OFF )
if(CVMIX_USE_NetCDF)
   add_compile_definitions(_NETCDF)
   include_directories($ENV{NetCDF_INCLUDE})
endif(CVMIX_USE_NetCDF)
option( CVMIX_BUILD_STATIC_LIBS "Build static library" ON )
option( CVMIX_BUILD_SHARED_LIBS "Build shared library" ON )
option( CVMIX_BUILD_DRIVER "Build CVMix example/test driver" OFF )

# Check if shared or static builds are turned off
if( NOT CVMIX_BUILD_STATIC_LIBS )
   message( STATUS "Turning STATIC builds OFF" )
endif()
if( NOT CVMIX_BUILD_SHARED_LIBS )
   message( STATUS "Turning SHARED builds OFF" )
endif()

# Set default build type to Release if not specified
if( NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE )
   message( STATUS "Setting default build type to 'Release'.  Set CMAKE_BUILD_TYPE variable to change build types." )
   set_property( CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release" )
endif()

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

add_subdirectory( src )

# Note - KB:
# Building static and shared libs require an ugly construct to build on both
# Windows and Linux. Windows will not create the libraries unless a Fortran
# file is explicitly given - and Linux complains (I think due to a race
# condition) if it is. Therefor the construct below - a proper solution is
# most welcome.

# Add static lib target
if( CVMIX_BUILD_STATIC_LIBS )
   if ( WIN32 )
	   add_library( cvmix_static STATIC ${CMAKE_SOURCE_DIR}/src/dummy.F90 $<TARGET_OBJECTS:cvmix_objects> )
   else()
      add_library( cvmix_static STATIC $<TARGET_OBJECTS:cvmix_objects> )
   endif()
   list( APPEND LIB_TARGETS cvmix_static )
endif()

# Add shared lib target
if( CVMIX_BUILD_SHARED_LIBS )
   if ( WIN32 )
	   add_library( cvmix_shared SHARED ${CMAKE_SOURCE_DIR}/src/dummy.F90 $<TARGET_OBJECTS:cvmix_objects> )
   else()
      add_library( cvmix_shared SHARED $<TARGET_OBJECTS:cvmix_objects> )
   endif()
   list( APPEND LIB_TARGETS cvmix_shared )
endif()

# Set common lib target properties
foreach( _lib IN LISTS LIB_TARGETS )
   target_compile_definitions( ${_lib} PUBLIC "${PUBLIC_FLAGS}" )
   target_include_directories( ${_lib}
      PUBLIC
         $<INSTALL_INTERFACE:include/cvmix>
      PRIVATE
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${INCLUDE_DIR}> )
   set_target_properties( ${_lib} PROPERTIES OUTPUT_NAME cvmix)
endforeach()

if(CVMIX_BUILD_DRIVER)

   add_executable( cvmix_driver . )
   target_sources( cvmix_driver PRIVATE src/cvmix_driver.F90 )
   set_property( TARGET cvmix_driver PROPERTY FOLDER driver )
   target_include_directories( cvmix_driver PRIVATE ${CMAKE_BINARY_DIR}/modules )
   target_link_libraries( cvmix_driver PRIVATE cvmix_drivers )

endif()

# Install
if(CVMIX_BUILD_DRIVER)
install( TARGETS cvmix_driver EXPORT cvmixConfig
	 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif()
install( TARGETS ${LIB_TARGETS} EXPORT cvmixConfig
	 ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/modules/ EXPORT cvmixConfig
	 DESTINATION  include/cvmix
	 FILES_MATCHING
	 PATTERN "*.mod" )
install(EXPORT cvmixConfig DESTINATION cmake)

