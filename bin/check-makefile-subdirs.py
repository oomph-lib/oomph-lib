#!/usr/bin/env python

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path
import glob

from os.path import join as pjoin


import parallel_self_test
from parallel_self_test import variable_from_makefile
from parallel_self_test import NoMakefileError

# ??ds todo: deal with mpi/eigenproblems/hypre being disabled in a smart
# way?


def recursive_find_makefile_ams(base_dirs):
    """Construct a list of directories to check in by searching for Makefile.am s."""

    sub_dirs = []
    for base in base_dirs:
        for root, dirs, files in os.walk(base):
            if 'Makefile.am' in files:
                sub_dirs.append(root)

    return sub_dirs


def list_subdirs_with_makefile_ams(base_dir):
    """Get a list of the sub directories of base_dir which contain a file
    called Makefile.am
    """

    sub_dirs = []
    for f in glob.glob(pjoin(base_dir, "*", "Makefile.am")):
        d = os.path.dirname(f)
        rel_d = os.path.relpath(d, base_dir)
        sub_dirs.append(rel_d)

    return sub_dirs


def check_makefileam(direc):
    """Compare subdirs variable in makefile in direc with actual sub
    directories of direc containing a makefile.am (i.e. which can be built).
    """
    
    # Extract the subdirs variable from the Makefile if one exists
    try:
        subdirs_var_string = variable_from_makefile("SUBDIRS",
                                                    pjoin(direc, "Makefile"))
    except NoMakefileError:
        print("Makefile not found in", direc)
        subdirs_var_string = ""

    subdirs_var = subdirs_var_string.split()

    # Check which subdirs really contain makefile.am s
    subdirs_actual = list_subdirs_with_makefile_ams(direc)

    # If theres a difference then report it
    if set(subdirs_var) != set(subdirs_actual):
        print("Not the same in:", direc, "\n")
        print("Makefile.am SUBDIRS variable is:", set(subdirs_var), "\n")
        print("But the actual list of compilable subdirs is:",
              set(subdirs_actual), "\n")

        return False

    else:
        return True


def main():
    """Compare subdirs variable in makefiles found recursively from the given
    root directory with actual sub directories of their directory
    containing a makefile.am (i.e. which can be built).
    """
    
    # Parse arguments
    # ============================================================
    
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     
    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--oomph-root', '-C', action = "store",
                        default="/home/david/oomph-lib")
    
    args = parser.parse_args()

    
    # Create list of folders to check
    makefile_ams_to_check = recursive_find_makefile_ams([args.oomph_root])

    # Folders where it is not expected to be true due to various tricks
    skiplist = [pjoin(args.oomph_root, "external_src"),
                pjoin(args.oomph_root, "external_distributions"),
                ]

    # Loop over directories and check the makefileams
    successes = map(check_makefileam,
                    [r for r in makefile_ams_to_check if r not in skiplist])


    # Return 0 if all ok, 1 otherwise    
    if all(successes):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
