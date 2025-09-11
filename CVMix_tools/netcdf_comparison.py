#!/usr/bin/env python

"""
    Usage:
        $ ./netcdf_comparison.py --baseline BASELINE_FILE --new-file NEW_FILE
                                 --strict {exact,loose} [-r RTOL] [-a ATOL] [-t THRES]

    Use xarray and numpy to compare two netcdf files.
    For each variable, flag
    1. Variables that are present in one file but not the other
    2. Variables where the data type doesn't match across files
    3. Variables where the dimensions don't match across files
    4. Variables where the missing values are not aligned
    5. Variables that differ in one of two ways (user specifies which strictness level to use):
       a. Variables that are not exactly the same (--strict exact)
       b. Variables that are not "close" to each other (--strict loose)
          -- For values very close to 0 (<THRES), variables that differ by
             more than ATOL are flagged
          -- For values larger than THRES variables with a relative difference
             of more than RTOL are flagged

    There is some knowledge of unit equivalence between mks and cgs.
"""

import logging

##################

# Store default values of rtol, atol, and thres in a global dictionary
# to make it easy to update the default values if necessary
# rtol = 1e-11 fails the cgs vs mks comparison
#              (POC_REMIN_DIC and PON_REMIN_NH4 have rel errors of ~1.4e-11)
DEFAULT_TOLS = {'rtol' : 1e-9, 'atol' : 1e-16, 'thres' : 1e-16}

##################

def ds_comparison_exact(ds_base, ds_new):
    """
        Compare baseline to new_file using xarray
        Will only pass if two files are identical (bfb)
        (run with --strict=exact)
    """
    # Compare remaining variables
    return ds_comparison_loose(ds_base, ds_new, rtol=0, thres=0)

##################

def ds_comparison_loose(ds_base, ds_new, rtol=DEFAULT_TOLS['rtol'], atol=DEFAULT_TOLS['atol'],
                        thres=DEFAULT_TOLS['thres']):
    """
        Compare baseline to new_file using xarray
        Will pass if two files are within specified tolerance
        (run with --strict=loose)
    """
    header_fail, ds_base, ds_new = _reduce_to_matching_variables(ds_base, ds_new)

    # Compare remaining variables
    return  _variable_check_loose(ds_base, ds_new, rtol, atol, thres) or header_fail

##################

def _open_files(baseline, new_file):
    """
    Reads two netCDF files, returns xarray datasets
    """
    import xarray as xr

    logger = logging.getLogger(__name__)
    logger.info("Comparing %s to the baseline %s\n", new_file, baseline)

    return xr.open_dataset(baseline), xr.open_dataset(new_file)

##################

def _reduce_to_matching_variables(ds_base, ds_new):
    """
        Are variables and dimensions the same in both datasets? Check for:
            1. Variables in one dataset but not the other
            2. Variables that are different types
            3. Variables that are different dimensions
               * check both name and size of dimension
        Returns header_fail, ds_base, ds_new
            * header_fail is True if any of three checks above fail
            * ds_base and ds_new are same as inputs, except any variables that do not match the
              above criteria are dropped from the dataset
    """
    logger = logging.getLogger(__name__)
    header_fail = False
    failed_vars = []

    # 1. Any variables in one file but not the other?
    base_vars = ds_base.variables
    new_vars = ds_new.variables
    common_vars = set(base_vars) & set(new_vars)
    base_vars = list(set(base_vars) - common_vars)
    new_vars = list(set(new_vars) - common_vars)
    if base_vars:
        header_fail = True
        ds_base = ds_base.drop_vars(base_vars)
        logger.info("The following variables are in the baseline file but not the new file:")
        base_vars.sort()
        for var in base_vars:
            logger.info("* %s", var)
        logger.info("")
    if new_vars:
        header_fail = True
        ds_new = ds_new.drop_vars(new_vars)
        logger.info("The following variables are in the new file but not the baseline file:")
        new_vars.sort()
        for var in new_vars:
            logger.info("* %s", var)
        logger.info("")

    # 2. Can variables be compared?
    common_vars = list(common_vars)
    common_vars.sort()
    for var in common_vars:
        err_messages = []
        # (a) Are they the same type?
        if ds_base[var].dtype != ds_new[var].dtype:
            err_messages.append(
                "Variable is {} in baseline and {} in new file".format(
                    ds_base[var].dtype, ds_new[var].dtype
                )
            )
        # (b) Do the dimension names match?
        if ds_base[var].dims != ds_new[var].dims:
            err_messages.append(
                "Baseline dimensions are {} and new file dimensions are {}".format(
                    ds_base[var].dims, ds_new[var].dims
                )
            )
        # (c) Do the dimension sizes match?
        if ds_base[var].data.shape != ds_new[var].data.shape:
            err_messages.append(
                "Baseline dimensions are {} and new file dimensions are {}".format(
                    ds_base[var].data.shape, ds_new[var].data.shape
                )
            )

        # Report errors
        if _report_errs(var, err_messages):
            header_fail = True
            failed_vars.append(var)

    return header_fail, ds_base.drop_vars(failed_vars), ds_new.drop_vars(failed_vars)

##################

def _get_conversion_factor(ds_base, ds_new, var):
    unit_conversion = {key: {} for key in
                       ['cm/s', 'cm','nmol/cm^3',  'nmol/cm^3/s', 'nmol/cm^2/s',
                        '(nmol/cm^3)^-1 s^-1', 'g/cm^3/s', 'g/cm^2/s', 'meq/m^3', 'meq/m^3/s',
                        'neq/cm^2/s', 'meq/m^3 cm/s', 'mg/m^3 cm/s']}
    unit_conversion['cm/s']['m/s'] = 0.01 # cm/s -> m/s
    unit_conversion['cm']['m'] = 0.01 # cm -> m
    unit_conversion['nmol/cm^3']['mmol/m^3'] = 1 # nmol/cm^3 -> mmol/m^3
    unit_conversion['nmol/cm^3/s']['mmol/m^3/s'] = 1 # nmol/cm^3/s -> mmol/m^3/s
    unit_conversion['nmol/cm^2/s']['mmol/m^2/s'] = 0.01 # nmol/cm^2/s -> mmol/m^2/s
    unit_conversion['(nmol/cm^3)^-1 s^-1']['(mmol/m^3)^-1 s^-1'] = 1 # same as nmol/cm^3 -> mmol/m^3
    unit_conversion['g/cm^3/s']['kg/m^3/s'] = 1000 # g/cm^3/s -> kg/m^3/s
    unit_conversion['g/cm^2/s']['kg/m^2/s'] = 10 # g/cm^2/s -> kg/m^2/s
    unit_conversion['meq/m^3']['neq/cm^3'] = 1 # meq/m^3 -> neq/cm^3
    unit_conversion['meq/m^3/s']['neq/cm^3/s'] = 1 # meq/m^3/s -> neq/cm^3/s
    unit_conversion['neq/cm^2/s']['meq/m^2/s'] = 0.01 # meq/m^3 cm/s -> meq/m^2/s
    unit_conversion['meq/m^3 cm/s']['neq/cm^2/s'] = 1 # meq/m^3 cm/s -> neq/cm^2/s
    unit_conversion['meq/m^3 cm/s']['meq/m^2/s'] = 0.01 # meq/m^3 cm/s -> meq/m^2/s
    unit_conversion['mg/m^3 cm/s']['mg/m^2/s'] = 0.01 # mg/m^3 cm/s -> mg/m^2/s

    conversion_factor = 1.
    if('units' in ds_base[var].attrs and ds_new[var].attrs):
        old_units = ds_base[var].attrs['units']
        new_units = ds_new[var].attrs['units']
        if old_units != new_units:
            found=False
            if new_units in unit_conversion and old_units in unit_conversion[new_units]:
                conversion_factor = unit_conversion[new_units][old_units]
                found = True
            if not found:
                if old_units in unit_conversion and new_units in unit_conversion[old_units]:
                    conversion_factor = 1. / unit_conversion[old_units][new_units]
                    found = True
            if not found:
                raise KeyError(f'Can not convert from {new_units} to {old_units} for {var}')
    return conversion_factor

def _variable_check_loose(ds_base, ds_new, rtol, atol, thres):
    """
        Assumes both datasets contain the same variables with the same dimensions
        Checks:
        1. Are NaNs in the same place?
        2. If baseline value = 0, then |new value| must be <= atol
        3. Absolute vs relative error:
           i. If 0 < |baseline value| <= thres then want absolute difference < atol
           ii. If |baseline value| > thres, then want relative difference < rtol
           Note: if thres = 0 and rtol = 0, then this reduces to an exact test
                 (Are non-NaN values identical?)
    """
    import numpy as np

    error_checking = {'var_check_count': 0}
    common_vars = list(set(ds_base.variables) & set(ds_new.variables))
    common_vars.sort()

    for var in common_vars:
        error_checking['messages'] = []

        # (0) Update units of new file to match baseline
        conversion_factor = _get_conversion_factor(ds_base, ds_new, var)

        # (1) Are NaNs in the same place?
        if var.lower() == 'time':
            continue
        mask = np.isfinite(ds_base[var].data)
        if np.any(mask ^ np.isfinite(ds_new[var].data)):
            error_checking['messages'].append('NaNs are not in same place')

        # (2) compare everywhere that baseline is 0
        if np.any(np.where(ds_base[var].data[mask] == 0,
                           np.abs(ds_new[var].data[mask]) > atol,
                           False)):
            error_checking['messages'].append(
                'Baseline is 0 at some indices where abs value ' +
                'of new data is larger than {}'.format(atol))

        # (3i) Compare everywhere that 0 < |baseline| <= thres
        base_data = np.where((ds_base[var].data[mask] != 0) &
                             (np.abs(ds_base[var].data[mask]) <= thres),
                             ds_base[var].data[mask], 0)
        new_data = np.where((ds_base[var].data[mask] != 0) &
                            (np.abs(ds_base[var].data[mask]) <= thres),
                            conversion_factor*ds_new[var].data[mask], 0)
        abs_err = np.abs(new_data - base_data)
        if np.any(abs_err > atol):
            error_checking['messages'].append("Max absolute error ({}) exceeds {}".format(
                np.max(abs_err), atol))

        # (3ii) Compare everywhere that |baseline| is > thres
        base_data = np.where(np.abs(ds_base[var].data[mask]) > thres,
                             ds_base[var].data[mask],
                             0)
        new_data = np.where(np.abs(ds_base[var].data[mask]) > thres,
                            conversion_factor*ds_new[var].data[mask],
                            0)
        rel_err = _compute_rel_err(ds_base[var], base_data, new_data, mask)
        if np.any(rel_err > rtol):
            if rtol == 0:
                abs_err = np.abs(new_data - base_data)
                error_checking['messages'].append("Values are not the same everywhere\n{}".format(
                    "    Max relative error: {}\n    Max absolute error: {}".format(
                        np.max(rel_err), np.max(abs_err))
                ))
            else:
                error_checking['messages'].append("Max relative error ({}) exceeds {}".format(
                    np.max(rel_err), rtol))

        error_checking['var_check_count'] += _report_errs(var, error_checking['messages'])

    return error_checking['var_check_count']>0

##################

def _compute_rel_err(da_base, base_data, new_data, mask):
    # denominator for relative error is local max value (3-point stencil)
    # note the assumption that column is first dimension
    import numpy as np

    if 'num_levels' in da_base.dims:
        # For variables with a depth dimension, we use a 3-point stencil in depth to find the
        # maximum value of the baseline value to use in the denominator of the relative error.
        # For the top (bottom) of the column, we use a 2-point stencil with the value
        # below (above) the given level
        rel_denom = da_base.rolling(num_levels=3, center=True, min_periods=2).max().data[mask]
    else:
        # For variables without a depth dimension, the denominator is the baseline value
        rel_denom = da_base.data[mask]
    rel_denom = np.where(np.isfinite(rel_denom), rel_denom, 0)
    return np.where(base_data != 0, np.abs(new_data - base_data), 0) / \
           np.where(rel_denom != 0, rel_denom, 1)

##################

def _report_errs(var, err_messages):
    """
        err_messages is a list of all accumulated errors
    """
    logger = logging.getLogger(__name__)
    if err_messages:
        logger.info("Variable: %s", var)
        for err in err_messages:
            logger.info("... %s", err)
        logger.info("")
        return True
    return False

##################

def _parse_args():
    """ Parse command line arguments
    """

    import argparse

    parser = argparse.ArgumentParser(description="Compare two netCDF files using xarray and numpy",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Baseline for comparison
    parser.add_argument('-b', '--baseline', action='store', dest='baseline', required=True,
                        help='Baseline for comparison')

    # File to compare to baseline
    parser.add_argument('-n', '--new-file', action='store', dest='new_file', required=True,
                        help="File to compare to baseline")

    # Tolerances
    parser.add_argument('--strict', choices=['exact', 'loose'], required=True,
                        help="Should files be bit-for-bit [exact] or within some tolerance [loose]")

    parser.add_argument('-r', '--rtol', action='store', dest='rtol',
                        default=DEFAULT_TOLS['rtol'], type=float,
                        help="Maximum allowable relative tolerance (only if strict=loose)")

    parser.add_argument('-a', '--atol', action='store', dest='atol',
                        default=DEFAULT_TOLS['atol'], type=float,
                        help="Maximum allowable absolute tolerance (only if strict=loose)")

    parser.add_argument('-t', '--thres', action='store', dest='thres',
                        default=DEFAULT_TOLS['thres'], type=float,
                        help="Threshold to switch from abs tolerance to rel (only if strict=loose)")

    return parser.parse_args()

##################

if __name__ == "__main__":
    import os

    # Set up logging
#    logging.basicConfig(format='%(levelname)s (%(funcName)s): %(message)s', level=logging.INFO)
    logging.basicConfig(format='%(message)s', level=logging.INFO)
    LOGGER = logging.getLogger(__name__)

    args = _parse_args()
    ds_base_in, ds_new_in = _open_files(args.baseline, args.new_file)
    if args.strict == 'loose':
        if ds_comparison_loose(ds_base_in, ds_new_in, args.rtol, args.atol, args.thres):
            LOGGER.error("Differences found between files!")
            os.sys.exit(1)
        LOGGER.info("PASS: All variables match and are within specified tolerance.")
    if args.strict == 'exact':
        if ds_comparison_exact(ds_base_in, ds_new_in):
            LOGGER.error("Differences found between files!")
            os.sys.exit(1)
        LOGGER.info("PASS: All variables match and have exactly the same values.")
