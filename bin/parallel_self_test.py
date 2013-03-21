#!/usr/bin/env python


# TODO:

# Different message for dummy oks due to missing code? But we can't see
# what's missing from outside the code...

# Interrupt handling? Would be nice for C-c to kill all processes not just
# the current one.

# How do hypre/trilinos interact with the tests?

# Caveats:

# Eigenproblems fail because there is no good way to determine if we have
# arpack or not.

# Doesn't record successes until it gets to the end. Would perhaps be nice
# if it did it immediately.


import subprocess as subp
from multiprocessing import Pool
import sys
import argparse
import os
import os.path
import multiprocessing
import itertools as it

from datetime import datetime as dt
from functools import partial as pt


class Colours:
    """Very simple class for coloring in text output. Stolen from
    http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python

    If we were using python 3 we could use the termcolor package instead."""
    Header = '\033[96m'
    Okgreen = '\033[92m'
    MakeFail = '\033[93m'
    TestFail = '\033[91m'
    Endc = '\033[0m'

    def disable(self):
        self.Header = ''
        self.Okgreen = ''
        self.MakeFail = ''
        self.TestFail = ''
        self.Endc = ''

# Make a GLOBAL colours object to tell print functions how to make things
# pretty. Easier to make it global because we want to be able to disable it
# from inside main() without passing around a Colours object all over the
# place...
COLOURS = Colours()

# Various functions for printing with pretty colours


def highlight(directorypath):
    return os.path.dirname(directorypath) + "/" + COLOURS.Header + \
        os.path.basename(directorypath) + COLOURS.Endc


def build_fail_message(directory):
    print COLOURS.MakeFail+"[BUILD FAIL] " + highlight(directory) + COLOURS.Endc


def check_fail_message(directory):
    print COLOURS.TestFail+"[FAILED]     " + highlight(directory) + COLOURS.Endc


def check_success_message(directory):
    print COLOURS.Okgreen + "[OK]         " + highlight(directory) + COLOURS.Endc


def no_check_message(directory):
    print COLOURS.Okgreen + "[NO CHECK]   " + highlight(directory) + COLOURS.Endc


def no_mpi_message(directory):
    print COLOURS.Okgreen + "[NO MPI]     " + highlight(directory) + COLOURS.Endc


# Some general utility functions
# ============================================================

def partition(pred, iterable):
    """Apply a predicate function "pred" to each item in "iterable". Output
    two lists, the first with all items for which the function returned
    true and the second list containing the rest.
    """
    trues = []
    falses = []
    for item in iterable:
        if pred(item):
            trues.append(item)
        else:
            falses.append(item)
    return trues, falses


def move_to_front(front_list, full_list):
    """Make a new list with anything in front_list first.

    Not at all efficient, but ok for now.
    """
    new_list = []
    for item in full_list:
        if any(s in full_list for s in front_list):
            new_list.insert(0, item)
        else:
            new_list.append(item)


def mygrep(file, searchstring):
    """(F)Grep for a line in file, return None if not found."""
    for line in open(file):
        if searchstring in line:
            return line
    return None

# Storing/reading "previous passes" lists
# ============================================================


def get_previous_passes(previous_run_filename):
    """Read list of previously passed tests from
    "$ROOT/demo_drivers/passed_tests_list"
    """

    try:
        with open(previous_run_filename, 'r') as f:
            passes_list = f.readlines()

    # If there is no file there then the list is empty. Use try/catch
    # rather than explcitly checking for a file to avoid race conditions.
    except IOError as e:
        passes_list = []

    return [dir.strip() for dir in passes_list]


def write_previous_passes(passes, previous_run_filename):
    """Write list of tests that passed to
    "$ROOT/demo_drivers/passed_tests_list"
    """

    with open(previous_run_filename, 'w') as f:
        for item in passes:
            f.write("%s\n" % item)

    return


# Validation functions
# ============================================================

def find_validate_dirs(base_dirs):
    """Construct a list of validation directories by searching for
    validate.sh scripts."""

    # Identify the really slow drivers so they can be pushed to the front of
    # the list.
    slow_driver_dirs = ['pseudo_solid_collapsible_tube',
                        'fsi_channel_seg_and_precond',
                        'fsi_collapsible_channel',
                        'vmtk_fsi']
    slow_mpi_driver_dirs = []

    all_validation_dirs = []
    for base in base_dirs:
        for root, dirs, files in os.walk(base):
            if 'validate.sh' in files:
                all_validation_dirs.append(root)

    # Seperate into mpi directory lists and other
    mpi_dirs, validation_dirs = partition(lambda x: "mpi" in x.lower(),
                                          all_validation_dirs)

    return validation_dirs, mpi_dirs


def mpi_cores_used(oomph_root):
    """ Find how many cores oomph-lib is configured to use for each mpi
    self test."""

    # Find the line in the Makefile which specifies the mpi run command.
    makefilepath = os.path.join(oomph_root, "Makefile")
    line = mygrep(makefilepath, "MPI_RUN_COMMAND")

    # Get out the command in list of strings format
    try:
        mpicommand = line.rstrip().rsplit(" = ", 1)[1].split()
    except IndexError:
        # If we get an index error there is no mpi command
        return None

    # Parse the mpi command for the -np argument
    parser = argparse.ArgumentParser(description='mpirun parser')
    parser.add_argument('--n', '-np', '-c', '-n', dest='np')
    args, unknown = parser.parse_known_args(mpicommand)

    # Check the mpi run command is ok (it must be 2)
    if int(args.np) != 2:
        print """##### MPIRUN COMMAND MUST TWO CORES (-np 2) #####"""
        print """Otherwise the meshes used are not valid."""
        sys.exit(11)

    # Return as an integer
    return int(args.np)


# The function doing the bulk of the actual work (called many times in
# parallel by main).
def make_check_in_dir(directory, previous_run_passes_list):
    """
    Rebuild binaries in the directory using make if needed. Then run the
    tests if we did a rebuild or if the tests failed previously.

    Since everything important is written to validation logs just summarise
    passes/fails on stdout.

    Output from compilation is sent to /dev/null since it is trivial to
    rerun make and get it again.
    """

    # Check if anything needs to rebuild (it has to be 'make check' because
    # 'make' alone does nothing with demo drivers).
    rebuild_return = subp.call(['make', 'check', '-q', '-C', str(directory)],
                               stdout=open(os.devnull, 'w'),
                               stderr=open(os.devnull, 'w'))

    # 'make -q' returns non-zero if we need to rebuild anything, get a bool
    # saying if it was non-zero.
    rebuilt = (rebuild_return != 0)

    # Rebuild but don't run the test yet.
    build_return = subp.call(['make', 'check', '-k',
                              'TESTS_ENVIRONMENT=true',
                              '-C', str(directory)],
                             stdout=open(os.devnull, 'w'),
                             stderr=open(os.devnull, 'w'))

    # If it failed then return a failure immediately
    if build_return != 0:
        build_fail_message(directory)
        return False

    # Check if previous run of the tests failed
    if directory not in previous_run_passes_list:
        previous_run_failed = True
    else:
        previous_run_failed = False

    # If we rebuilt anything or if this test didn't pass earlier then run
    # 'make check'.
    if rebuilt or previous_run_failed:
        test_result = subp.call(['make', 'check', '-C', str(directory)],
                                stdout=open(os.devnull, 'w'),
                                stderr=open(os.devnull, 'w'))
        if test_result == 0:
            check_success_message(directory)
            return directory
        else:
            check_fail_message(directory)
            return False

    # Otherwise no test is needed
    else:
        no_check_message(directory)
        return directory


# Function dealing with parsing of arguments, getting everything set up and
# dispatching jobs to processes.
def main():
    """ Run self tests in parallel (one per core for serial tests, one per
        two cores for mpi). Only run the tests if a rebuild was needed or
        if they failed last time.

        Note: to abort C-c will not work (due to limitations of
        python's multiprocessing module). If you want to abort you will
        probably have to close the terminal emulator entirely.

        Also note: eigensolver tests will fail if you don't have arpack
        installed. I haven't found a way to not run them without arpack
        yet.
        """

    # Parse inputs
    # ============================================================
    parser = argparse.ArgumentParser(description=main.__doc__)

    # make uses -C for "run make in this directory"
    parser.add_argument('-C', default="/home/david/oomph-lib",
                        dest='oomph_root',
                        help='Set the root directory of oomph-lib.')

    parser.add_argument('-a', action='store_true', dest='makeclean',
                        help='Run make clean on test folders before starting.')

    parser.add_argument('-n', '-j', dest='ncores',
                        help='Specifiy how many cores should be used altogether \
                        (taking mpi runs into account). By default use all cores.',
                        default=multiprocessing.cpu_count())

    parser.add_argument('--no-colour', action='store_true', dest='no_colours',
                        help='Disable colours in output.')

    args = parser.parse_args()

    if args.no_colours:
        COLOURS.disable()

    # Construct the filename where data on what passed last time is stored.
    previous_run_filename = os.path.join(
        args.oomph_root, "demo_drivers/previous_run_passes_list")

    # Gather directory lists etc.
    # ============================================================
    # Create a list of absolute paths to demo driver directories
    base_dirs = [os.path.join(args.oomph_root, base_dir)
                 for base_dir in ["demo_drivers", "self_test"]]

    # Clean up from past runs if requested
    if args.makeclean:
        for dir in base_dirs:
            subp.check_call(['make', 'clean', '-k',
                             '-C', str(dir), '-j', str(args.ncores)])

        # Kill the old "previous_run_passes_list" in case we cancel
        # before writing out a new one. It's ok if it doesn't exist.
        try:
            os.remove(previous_run_filename)
        except OSError:
            pass


    # Grab data from past runs and wrap into test function
    previous_run_passes_list = get_previous_passes(previous_run_filename)
    check_func = pt(make_check_in_dir,
                    previous_run_passes_list=previous_run_passes_list)

    # Construct a list of validation directories
    validation_dirs, mpi_dirs = find_validate_dirs(base_dirs)

    # Run tests
    # ============================================================
    # Start all non-mpi tests first. Run "make_check_in_dir" on all
    # directories in "validation_dirs".
    pool = Pool()
    test_results = pool.map_async(check_func, validation_dirs).get(None)
    pool.close()                # Tell python there are no more jobs coming

    # Now run the MPI tests
    basempicores = mpi_cores_used(args.oomph_root)
    if basempicores is not None:
        # Figure out how many cores to use. Each mpi test uses some number
        # of cores per test already.
        mpi_ncores = int(int(args.ncores) / basempicores)

        # Now run the mpi checks
        mpipool = Pool(mpi_ncores)
        mpi_test_results = mpipool.map_async(check_func, mpi_dirs).get(None)
        mpipool.close()         # Tell python there are no more jobs coming
        mpipool.join()          # Wait for everything to finish

    else:
        # Otherwise just print a list of tests that were not run (with a
        # useful message).
        map(no_mpi_message, mpi_dirs)
        mpi_test_results = []

    pool.join()       # Wait for everything to finish

    # Record what passed (collect all results and remove any "False" then
    # write to file).
    all_results = it.ifilter(lambda x: x,
                             it.chain(test_results, mpi_test_results))
    write_previous_passes(all_results, previous_run_filename)

    # Done!
    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
