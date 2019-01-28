Last update: July 23, 2014
--------------------------


CVMix is a transparent, robust, flexible, and well-documented open source
Fortran library for parameterizing ocean vertical mixing processes. It contains
the following first-order closures:
1. Background [time-independent] mixing, including the Bryan and Lewis (1979)
   parameterization
2. Convective mixing (both density-based and the Brunt-Vaisala scheme)
3. Double diffusion mixing, as described in Large, et al. (1994) and modified
   according to Danabasoglu, et al. (2006)
4. The KPP boundary layer mixing scheme from Large, et al. 1994
5. Shear [Richardson number-based] mixing, including both Pacanowski and
   Philander (1981) and the shear mixing used in Large, et al. (1994)
6. Tidal mixing, which currently only offers the method described in Simmons,
   et al. (2004)

The library can be built on its own and linked in to a general ocean
circulation model. Additionally, CVMix is packaged with a very light-weight
stand-alone driver that offers the ability to run quick (typically single
column) tests to ensure everything is working properly.

Community input is both welcomed and encouraged. Documentation online at

https://github.com/CVMix/CVMix-src/wiki

goes into details about how to contribute, but basically we ask three things:
1. Follow the CVMix coding policies to ensure all modules use the same
   formatting and all documentation is consistent.
2. Run the provided suite of regression tests to ensure you did not
   accidentally change something in another module.
3. Provide a stand-alone driver for your method so that others can ensure they
   do not accidentally change something in your module.


INSTALLATION NOTES
------------------

The src directory contains a Makefile and a simple 'make' should be sufficient
to build the standalone driver. The first time you build, the 'cvmix_setup'
utility will run and prompt you for compiler and netcdf information - it will
only run once, and the info is saved in bld/.CVMix_env.

To build with netcdf, run 'make netcdf'.

The default executable is $CVMix/bin/cvmix, but it can be overwritten with
'make EXE=[executable]'.

To build just the library, run the following (from src/shared):

$ make FC=[compiler]                        \
       FCFLAGS=[compiler flags]             \
       INC_DIR=[place to save .mod files]   \
       LIB_DIR=[place to create libcvmix.a] \
       OBJ_DIR=[place to save .o files]

And then use -I$(INC_DIR) -L$(LIB_DIR) -lcvmix when you build software using
the CVMix library.


DIRECTORY STRUCTURE
-------------------

bin/ -- Default location for the cvmix executable.

bld/ -- Contains auxiliary files needed by the build system. CompileFlags.mak
        has default compile flags for 5 different compilers -- gfortran, 
        pgf90, ifort, xlf90, and nagfor, as well as ftn (the Cray wrapper for
        pgf90). At this time, no other compilers are supported on Cray systems.
        cvmix_setup is a python script that saves information about what
        compiler to use and where netCDF is installed in .CVMix_env.

  bld/obj -- Where .o  and .mod files for the stand-alone drivers are stored.

doc/ -- Contains documentation. At this point, it only has a PDF of the latest
        protex in-code notes and a script to generate said PDF.

git_config/ -- Contains files for recommended git setup of repository. At this
               time, only file is git_commit_template.txt, and running

               $ git config commit.template git_config/git_commit_template.txt

               Will add comments to default commit log to help organize the
               commit log.

include/ -- Default location for the .o and .mod files that need to be included
            with libcvmix.a.

inputdata/ -- The stand-alone tidal driver requires some large netCDF files to
              provide initial conditions and grid information. These files are
              downloaded the first time you run Simmons-test.sh (due to the
              file size, they are not kept in the same repository as the source
              code).

lib/ -- Default location for libcvmix.a, the library that can be linked in
        to an ocean model.

reg_tests/ -- Contains functional tests for all mixing methods except
              convective. Each test directory contains a shell script that
              calls functions defined in common/ to determine what settings to
              use for the test. All scripts have a "-h" (or "--help") flag to
              print usage options and a "-nc" (or "--netcdf") flag to turn on
              netCDF output. Note that the "-nc" is required for the tidal
              mixing test because it reads in netcdf data.

              All directories contain one (or many) NCL script(s) to produce a
              plot from the output; by default these plots are in the PDF
              format.

  reg_tests/Bryan-Lewis/ -- location of the Bryan-Lewis background mixing test
              $ ./BL_test.sh

              * plot_diff_coeffs.ncl -- makes a plot with depth on the y-axis
                and diffusivity (for both columns) on the x-axis.

  reg_tests/double_diff/ -- location of the double diffusion test
              $ ./double_diff-test.sh

              * plot_diff_coeffs.ncl -- makes two plots
                1. ddiff-diffuse: shows the interior diffusivity for potential
                   temperature (normalized by the molecular diffusivity)
                   in the diffusive convection regime plotted against the
                   inverse stratification parameter.
                2. ddiff-salt: shows the interior diffusivity for salinity
                   (normalized by its maximum) in the salt fingering regime
                   plotted against the stratification parameter.
                Note: these plots are comparable to the two plots in Figure 4
                      of Large et al.

  reg_tests/kpp/ -- location of the KPP mixing test
              $ ./kpp-test.sh

              * compare_interpolation_methods.ncl -- makes three plots
                1. bldepth: compares different interpolation techniques
                   (linear, quadratic, and cubic) when it comes to estimating
                   the boundary layer depth based on the depth of the mix
                   layer (analytically, HBL = HMIX+2)
                2. bldepth_error: the error of the three different
                   interpolation techniques compared to the analytic solution
                3. single_col22.05.pdf: the estimated boundary layer depth
                   using each interpolation method when HMIX=22.05 (shows
                   bulk Richardson number on x-axis, HBL is point where bulk
                   Richardson = 0.3)
              * plot_bulk_Rich.ncl -- makes the KPP-bulk_Rich plot, which
                shows buoyancy, zonal wind, and bulk Richardson values at
                each layer center. Compare to Figure C1 of Large et al.
              * plot_flux_profiles.ncl -- makes the KPP-flux_profile plot,
                which shows the non-dimensional flux profiles for momentum
                and scalars as a function of the stability parameter. Compare
                to Figure B1 of Large et al.

  reg_tests/shear-KPP/ -- location of the shear mixing test
              $ ./Large_test.sh

              * plot_diff_coeffs.ncl -- makes a plot with the local gradient
                Richardson number on the x-axis and the temperature
                diffusivity coefficient (normalized by its maximum) on the
                y-axis. Compare to Figure 3 of Large et al.

  reg_tests/tidal-Simmons/ -- location of the tidal mixing test
              $ ./Simmons_test.sh -nc

              * plot_diff_coeffs.ncl -- shows how the tidal diffusivity
                changes over depth for a specified column (column is set in
                the Simmons_nml namelist)
              * plot_tidal_energy_map.ncl -- produces a global contour map
                for the tidal energy flux (taken from the inputdata/
                directory)

  reg_tests/common/ -- a set of bash shell scripts that are used by the
                       individual scripts mentioned above to reduce the amount
                       of duplicate code from test to test.

src/ -- Contains the source code, organized as follows. The top directory
        contains modules needed by the stand-alone driver (output, for example)
        as well as the driver routine itself. Also contains the Makefile used
        to build the cvmix executable.

  src/drivers/ -- Subroutines called by the driver (one per test).

  src/shared/  -- Where all the modules that are needed to use CVMix with an 
                  outside model are stored. Also contains the Makefile used to
                  build the libcvmix.a library.

templates/ -- Contains a template F90 file for adding a new mixing module to
              the src/shared/ directory. See the CVMix webpage for directions
              for requesting to have your module included in the main CVMix
              repository.


ABOUT THE STAND-ALONE CVMIX DRIVER
----------------------------------

There are five options for the stand-alone driver which can be used to test two
different mixing methods. The first test is to output Bryan-Lewis mixing on two
columns (a high latitude column and a tropical column), and output the tracer
diffusivity coefficients at each level. This tests both the Bryan-Lewis portion
of the background mixing module and the ability for the CVMix data type (used
for the "wrapped" interface of all of the coefficient computation routines) to
use pointers or memory copies to store column data.

The second driver tests double diffusion mixing by setting up two columns, one
in the salt fingering regime and one in the diffusive convective instability
regime. In each case, the stratification parameter is set and then the
temperature and salinity diffusivity coefficients are computed.

The third driver runs a set of 6 tests for KPP boundary layer mixing:
  1. Compute the boundary layer depth for a column with a given HMIX (HBLT
     should = |HMIX|+2)
  2. Compute shape function coefficients with G(0) = G(1) = G'(1) = 0 and G'(0)
     = 1; should be sigma*(sigma-1)^2
  3. Compute nondimensional flux profile over range of stability parameter
     values
  4. Compute diffusivity in the boundary layer when boundary layer is (a) below
     cell center and (b) above cell center
  5. Compute bulk Richardson number from buoyancy and zonal wind data
  6. Compute velocity scales (a) when surface buoyancy flux is 0 and friction
     is positive; and (b) when surface buoyancy flux is negative and friction
     velocity is 0.

The fourth driver sets up a single column using the shear mixing formula found in
the Large, et al. KPP paper. For this test, each level is a different local
gradient Richardson number rather than a different depth to show how tracer
diffusivity varies with different Richardson numbers.

The fifth driver reads in tidal energy flux and buoyancy frequency data from an
input data set and computes diffusivities in a column based on these physical
values.
