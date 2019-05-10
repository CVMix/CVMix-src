#!/usr/bin/env python

"""
This script flags lines that do not conform to CVMix's coding practices
"""

import os
import sys
import logging
from collections import deque # Faster pop / append than standard lists
from collections import OrderedDict
UNIBLOCK = u"\u2588"

##############

class ConsistencyTestClass(object):
    """
    This class contains all the Fortran consistency check tests.
    The logs dictionary contains a deque object with a list of all source lines
    that fail a particular test (given as the logs key)
    """
    def __init__(self):
        self.logs = OrderedDict()

    ##############

    def process(self):
        """
        Check results from all tests that have been run:
        1. For each test:
            i.  log (as info) test description and number of lines of code that fail the test
            ii. log (as info) all lines of code that do not conform to standards
        2. return total number of errors across all tests
        """
        tot_err_cnt = 0
        while self.logs:
            desc, log = self.logs.popitem(last=False)
            err_cnt = len(log)
            LOGGER.info("* %s: %d error(s) found", desc, err_cnt)
            while log:
                msg = log.popleft()
                LOGGER.info("  %s", msg)
            tot_err_cnt += err_cnt
        return tot_err_cnt
    ##############

    def init_log(self, description):
        """
        If self.logs does not already have a key for description, set
        self.logs[description] to an empty deque list.
        """
        if description not in self.logs.keys():
            self.logs[description] = deque([])

    ##############

    def check_for_hard_tabs(self, file_and_line, line):
        """
        Log any lines containing hard tabs
        """
        test_desc = 'Check for hard tabs'
        self.init_log(test_desc)
        if "\t" in line:
            self.logs[test_desc].append("%s: %s" % (file_and_line,
                                                    line.replace("\t", 2*UNIBLOCK)))

    ##############

    def check_for_trailing_whitespace(self, file_and_line, line):
        """
        Log any lines containing trailing whitespace
        """
        test_desc = 'Check for trailing white space'
        self.init_log(test_desc)
        full_line_len = len(line)
        no_trailing_space_len = len(line.rstrip(" "))
        if no_trailing_space_len < full_line_len:
            self.logs[test_desc].append("%s: %s" % (file_and_line,
                                                    line.rstrip(" ") +
                                                    (full_line_len - no_trailing_space_len) *
                                                    UNIBLOCK))

    ##############

    def check_line_length(self, file_and_line, line, comment_char="!", max_len=132):
        """
        Log any lines exceeding max_len
        Currently ignores comment characters / anything following
        """
        test_desc = 'Check length of lines'
        self.init_log(test_desc)
        line_len = len(line.split(comment_char)[0].rstrip(" "))
        if line_len > max_len:
            self.logs[test_desc].append("%s: %s" % (file_and_line,
                                                    line[:max_len] +
                                                    (line_len-max_len) * UNIBLOCK))

    ##############

    def check_case_in_module_statements(self, file_and_line, line, comment_char="!"):
        """
        The following module statements should be all lowercase:
        * implicit none
        * public
        * private
        * save
        Note that at some point we may want to remove "save" altogether,
        since it is implicit in a module
        """
        test_desc = 'Check for case sensitive statements'
        self.init_log(test_desc)
        statements = ["implicit none", "public", "private", "save"]
        # ignore comments, and strip all white space
        line_without_comments = line.split(comment_char)[0].strip(" ")
        if line_without_comments.lower() in statements:
            if line_without_comments not in statements:
                self.logs[test_desc].append("%s: %s" % (file_and_line, line))

    ##############

    def check_for_spaces(self, file_and_line, line, comment_char="!"):
        """
        The following statements should all include spaces:
        * else if
        * end if
        * else where
        * end where
        * end do
        * end module
        """
        test_desc = 'Check for spaces in statements'
        self.init_log(test_desc)
        statements = ["elseif", "endif", "elsewhere", "endwhere", "enddo", "endmodule"]
        # ignore comments, and strip all white space
        line_without_comments = line.split(comment_char)[0].strip(" ")
        for bad_statement in statements:
            if line_without_comments.lower().startswith(bad_statement):
                self.logs[test_desc].append("%s: %s" % (file_and_line, line_without_comments))
                break

    ##############

    def check_for_double_quotes(self, file_and_line, line):
        """
        All Fortran strings should appear as 'string', not "string"
        """
        test_desc = 'Check for double quotes in statements'
        self.init_log(test_desc)
        if '"' in line:
            self.logs[test_desc].append("%s: %s" % (file_and_line, line))

    ##############

    def check_logical_statements(self, file_and_line, line):
        """
        Use symbols, not words, for logical operators:
        * >=, not .ge.
        * >, not .gt.
        * <=, not .le.
        * <, not .lt.
        * ==, not .eq.
        * /=, not .ne.
        """
        test_desc = 'Check for unwanted logical operators'
        self.init_log(test_desc)
        operators = ['.ge.', '.gt.', '.le.', '.lt.', '.eq.', '.ne.']
        for operator in operators:
            if operator in line:
                self.logs[test_desc].append("%s: %s" % (file_and_line, line))
                break

    ##############

    def check_r8_settings(self, file_and_line, line, comment_char="!"):
        """
        Make sure all real numbers are cast as r8
        """
        test_desc = 'Check for r8'
        self.init_log(test_desc)
        import re
        # Looking for decimal numbers that do not end in _r8
        # Edge cases:
        # 1. ignore decimals in comments
        # 2. Ignore decimals in format strings (e.g. E24.16)
        # 3. Allow 1.0e-2_r8
        line_without_comments = line.split(comment_char)[0].strip(" ")

        # Regex notes
        # 1. (?<!\w) -- do not match numbers immediately following a letter,
        #    such as E24 or I0 (these are Fortran format strings)
        #    [regex refers to this as a negative lookbehind]
        # 2. \d+\. -- match N consecutive decimal digits (for N>=1) followed
        #    by a decimal point
        # 3. (\d+[eE]([+-])?)? -- Optionally match a number followed by either
        #    e or E (optionally followed by + or -)
        # 4. \d+ -- match N consecutive decimal digits (for N>=1)
        # 5. (?!\d+|_|[eE]) -- do not match a decimal, an underscore, an e,
        #    or an E
        #    [regex refers to this as a negative lookahead]
        regex = r'(?<!\w)\d+\.(\d+[eE]([+-])?)?\d+(?!\d+|_|[eE])'
        valid = re.compile(regex)
        if valid.search(line_without_comments):
            self.logs[test_desc].append("%s: %s" % (file_and_line, line))

##############

if __name__ == "__main__":
    # CVMIX_ROOT is the top-level CVMix directory, which is a level above the
    # directory containing this script
    # * Do some string manipulation to make CVMIX_ROOT as human-readable as
    #   possible because it appears in the output if you run this script with
    #   the "-h" option
    SCRIPT_DIR = os.path.dirname(sys.argv[0])
    if SCRIPT_DIR == '.':
        CVMIX_ROOT = '..'
    elif SCRIPT_DIR.endswith('CVMix_tools'):
        CVMIX_ROOT = SCRIPT_DIR[:-12]
    else:
        CVMIX_ROOT = os.path.join(SCRIPT_DIR, '..')

    fortran_files = [] # pylint: disable=C0103
    python_files = []  # pylint: disable=C0103
    for root, dirs, files in os.walk(CVMIX_ROOT):
        for thisfile in files:
            if thisfile.endswith(".F90"):
                fortran_files.append(os.path.join(root, thisfile))
            elif thisfile.endswith(".py"):
                python_files.append(os.path.join(root, thisfile))

    # Use logging to write messages to stdout
    logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    LOGGER = logging.getLogger(__name__)

    Tests = ConsistencyTestClass() # pylint: disable=C0103

    # Fortran error checks
    LOGGER.info("Check Fortran files for coding standard violations:")
    for thisfile in fortran_files:
        with open(thisfile, "r") as fortran_file:
            line_cnt = 0
            for thisline in fortran_file.readlines():
                line_without_cr = thisline.rstrip("\n")
                line_cnt = line_cnt + 1
                file_and_line_number = "%s:%d" % (thisfile, line_cnt)
                Tests.check_for_hard_tabs(file_and_line_number, line_without_cr)
                Tests.check_for_trailing_whitespace(file_and_line_number, line_without_cr)
                Tests.check_line_length(file_and_line_number, line_without_cr)
                Tests.check_case_in_module_statements(file_and_line_number, line_without_cr)
                #Tests.check_for_spaces(file_and_line_number, line_without_cr)
                #Tests.check_for_double_quotes(file_and_line_number, line_without_cr)
                #Tests.check_logical_statements(file_and_line_number, line_without_cr)
                #Tests.check_r8_settings(file_and_line_number, line_without_cr)
    FORTRAN_ERROR_COUNT = Tests.process()
    LOGGER.info("Fortran errors found: %d", FORTRAN_ERROR_COUNT)

    # Python error checks
    LOGGER.info("\nCheck python files for coding standard violations:")
    for thisfile in python_files:
        with open(thisfile, "r") as python_file:
            line_cnt = 0
            for thisline in python_file.readlines():
                line_without_cr = thisline.rstrip("\n")
                line_cnt = line_cnt + 1
                file_and_line_number = "%s:%d" % (thisfile, line_cnt)
                Tests.check_for_hard_tabs(file_and_line_number, line_without_cr)
                Tests.check_for_trailing_whitespace(file_and_line_number, line_without_cr)
    PYTHON_ERROR_COUNT = Tests.process()
    LOGGER.info("Python errors found: %d", PYTHON_ERROR_COUNT)

    if FORTRAN_ERROR_COUNT + PYTHON_ERROR_COUNT > 0:
        LOGGER.info("\nTotal error count: %d", FORTRAN_ERROR_COUNT + PYTHON_ERROR_COUNT)
        sys.exit(1)

    LOGGER.info("\nNo errors found!")
