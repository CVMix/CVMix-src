# CVMix

This is a copy of [CVMix](https://github.com/CVMix/CVMix-src) to test Langmuir turbulence parameterizations in KPP. New features as compared to the standard version include the following.

- Langmuir turbulence enhanced entrainment according to [Li and Fox-Kemper (2017)](https://doi.org/10.1175%2FJPO-D-17-0085.1).

- New functions in `cvmix_kpp.F90`:

  - `cvmix_kpp_efactor_model(u10, ustar, hbl, CVmix_params_in)` to estimate the enhancement factor from empirical wave spectra ([Li et al., 2017](https://doi.org/10.1016%2Fj.ocemod.2017.03.016)). This function has already been included in the standard version of CVMix.

  - `cvmix_kpp_ustokes_SL_model(u10, hbl, CVmix_params_in)` to estimate the surface layer averaged Stokes drift from empirical wave spectra for the Langmuir turbulence enhanced entrainment.  

  - `cvmix_kpp_efactor_read(infile, lon, lat, time)` to read the enhanncement factor climatology from a file ([Li et al., 2017](https://doi.org/10.1016%2Fj.ocemod.2017.03.016)).

- A driver and test cases for testing Langmuir turbulence parameterizations in KPP. _Under development._