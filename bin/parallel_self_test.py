#!/usr/bin/env python


# TODO:

# Different message for dummy oks due to missing code? But we can't see
# what's missing from outside the code...

# Interrupt handling? Would be nice for C-c to kill all processes not just
# the current one.

# How do hypre/trilinos interact with the tests?

# Replace searching for mpi string in dir name with checking for
# $MPI_RUN_COMMAND in validate.sh?

# Print timing information?

# Caveats:

# Eigenproblems fail because there is no good way to determine if we have
# arpack or not.




# ASSUMPTIONS:

# All validate.sh scripts return an appropriate exit status. Otherwise this
# script cannot know if the test failed!

# All mpi drivers have "mpi" in their path name. If this is not true it
# should still work but things might be slow because we will have more
# processes than cores.

# All mpi drivers use the number of cores specified by MPI_RUN_COMMAND.

# MPI_RUN_COMMAND has a '-np' argument specifying the number of cores.


import subprocess as subp
from multiprocessing import Pool
import sys
import argparse
import os
import os.path
import multiprocessing
import itertools as it

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
# place.
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


def split_validation_dirs_mpi(all_validation_dirs):
    """ Split a list of validation directories into those containing "mpi"
    (case-insensitive) and the rest. e.g.

    serial_dirs, mpi_dirs = split_validation_dirs_mpi(all_validation_dirs)
    """
    return partition(lambda x: not "mpi" in x.lower(), all_validation_dirs)


def mygrep(file, searchstring):
    """(F)Grep for a line in file, return None if not found."""
    for line in open(file):
        if searchstring in line:
            return line
    return None


def variable_from_makefile(variable_name, makefile_path="Makefile"):
    """Extract a variable from a makefile (using make)."""

    if not os.path.isfile(makefile_path):
        raise IOError

    # Run make using both a real makefile and a dummy one we are about to
    # create. Call the "print-var" command which will be in the dummy
    # makefile.
    process = subp.Popen(["make", "-f", "-", "-f", makefile_path, "print-var"],
                         stdin=subp.PIPE, stdout=subp.PIPE)

    # Send in a dummy makefile with a print command, output should be the
    # value of the variable.
    stdout, _ = process.communicate("print-var:; @echo $(" + variable_name + ")")

    # Check that make exited with a success code (0)
    returncode = process.wait()
    if returncode != 0:
        raise subp.CalledProcessError

    # Get rid of trailing newline/whitespace and return
    return stdout.rstrip()


def error(*args):
    """Write an error message to stderr."""
    sys.stderr.write("\nERROR:\n" + "\n".join(args) + "\n")
    sys.exit(2)

# Validation functions
# ============================================================

def find_validate_dirs(base_dirs):
    """Construct a list of validation directories by searching for
    validate.sh scripts."""

    all_validation_dirs = []
    for base in base_dirs:
        for root, dirs, files in os.walk(base):
            if 'validate.sh' in files:
                all_validation_dirs.append(root)

    # Seperate into mpi directory lists and other
    return split_validation_dirs_mpi(all_validation_dirs)


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
def make_check_in_dir(directory):
    """
    Rebuild binaries in the directory using make if needed then run the
    tests.

    Since everything important is written to validation logs so just
    summarise passes/fails on stdout.

    Output from compilation is sent to /dev/null since it is trivial to
    rerun make and get it again.
    """

    # Rebuild (but don't run the test yet).
    build_return = subp.call(['make', 'check',
                              'TESTS_ENVIRONMENT=true'],
                              cwd = directory,
                              stdout=open(os.devnull, 'w'),
                              stderr=open(os.devnull, 'w'))

    # If it failed then return a failure immediately
    if build_return != 0:
        build_fail_message(directory)
        return

    # Run make check
    test_result = subp.call(['make', 'check'],
                            cwd = directory,
                            stdout=open(os.devnull, 'w'),
                            stderr=open(os.devnull, 'w'))
    if test_result == 0:
        check_success_message(directory)
        return
    else:
        check_fail_message(directory)
        return

# Function dealing with parsing of arguments, getting everything set up and
# dispatching jobs to processes.
def main():
    """
    Run self tests in parallel (one per core for serial tests, one per
    two cores for mpi).

    Note: aborting with C-c will not work due to limitations of python's
    multiprocessing module. The easiest way to abort is probably to kill
    the terminal emulator itself. Alternatively background the python
    process with C-z then kill it (i.e. run "kill %%").

    Also note: eigensolver tests will fail if you don't have arpack
    installed. I haven't found a way to not run them without arpack
    yet.

    Typical usage is just:

        parallel_self_test.py

    from within an oomph-lib directory, or

        parallel_self_test.py -C /abs/path/to/oomph/root

    from anywhere.


    For this script to work correctly ALL validate.sh scripts must return
    an exit status. A simple command to check for this is:

       find -name "validate.sh" | xargs grep -i -L "exit" | xargs grep -i -L "set -o errexit"
    """

    # Parse inputs
    # ============================================================
    parser = argparse.ArgumentParser(
        description=main.__doc__,

        # Don't mess up my formating in the help message
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    # make uses -C for "run make in this directory"
    parser.add_argument('-C', dest='oomph_root',
                        help='Set the root directory of oomph-lib, by default try to \
    extract it from a Makefile.')

    parser.add_argument('-a', action='store_true', dest='makeclean',
                        help='Run make clean on test folders before starting.')

    parser.add_argument('-j', '-n', dest='ncores',
                        help='Specifiy how many cores should be used altogether \
                        (taking mpi runs into account). By default use all cores.',
                        default=multiprocessing.cpu_count())

    parser.add_argument('--no-colour', action='store_true', dest='no_colours',
                        help='Disable colours in output.')

    args = parser.parse_args()

    if args.no_colours:
        COLOURS.disable()


    # Attempt to get the oomph-lib root dir from a Makefile in the current
    # directory.
    if args.oomph_root is None:
        try:
            print("Trying to extract the root dir from a Makefile in the pwd.")
            args.oomph_root = variable_from_makefile("abs_top_srcdir")
            print("Extracted oomph_root = " + args.oomph_root)

            # If no makefile exists we are stuck:
        except IOError:
            error("You must either run this script from within an oomph-lib",
                  "directory or specify the path to oomph_root using -C.")


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

    # Construct a list of validation directories
    validation_dirs, mpi_dirs = find_validate_dirs(base_dirs)

    # Run tests
    # =========================================================================
    # Start all non-mpi tests first. Run "make_check_in_dir" on all
    # directories in "validation_dirs". Set chunksize to 1 (i.e. each "make
    # check" call is in its own "chunk of work") to avoid the situation
    # where multiple slow "make check"s end up in the same chunk and we
    # have to wait ages for it to finish.
    pool = Pool()
    pool.map(make_check_in_dir, validation_dirs, 1)
    pool.close()                # Tell python there are no more jobs coming
    pool.join()       # Wait for everything to finish

    # Now run the MPI tests
    basempicores = mpi_cores_used(args.oomph_root)
    if basempicores is not None:
        # Figure out how many cores to use. Each mpi test uses some number
        # of cores per test already.
        mpi_ncores = int(int(args.ncores) / basempicores)

        # Now run the mpi checks
        mpipool = Pool(mpi_ncores)
        mpipool.map(make_check_in_dir, mpi_dirs, 1)
        mpipool.close()         # Tell python there are no more jobs coming
        mpipool.join()          # Wait for everything to finish

    else:
        # Otherwise just print a list of tests that were not run (with a
        # useful message).
        map(no_mpi_message, mpi_dirs)

    # Done!
    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
