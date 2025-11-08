#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import sys

# Handles un-indenting multiline strings
from textwrap import dedent


def run_fpdiff(oomph_lib_root_path: str,
               output_data_file: str,
               validata_file: str,
               fpdiff_tol: float,
               fpdiff_reltold: float,
               log_file: "TextIOWrapper"):

    subprocess([f"{oomph_lib_root_path}/scripts/fpdiff.py",
                validata_file,
                output_data_file,
                f"{fpdiff_tol}",
                f"{fpdiff_reltold}"],
               stdout=log_file)


def run_compare_file_length_with_tolerance(
        oomph_lib_root_path: str,
        output_data_file: str,
        validata_file: str,
        threshold_for_number_of_iterations: int,
        log_file: "TextIOWrapper"):

    subprocess([f"{oomph_lib_root_path}/scripts/compare_file_length_with_tolerance.bash",
                output_data_file,
                validata_file,
                f"{fpdiff_tol}",
                f"{fpdiff_reltold}"],
               stdout=log_file)


def test_two_d_linear_elasticity_with_simple_block_diagonal_preconditioner(
        oomph_lib_root_path: str,
        test_directory: str,
        log_file: "TextIOWrapper",
        have_fpdiff: str) -> None:
    """Validation for simple block preconditioner for linear elasticity."""
    # If we're not in the test directory, go there now
    if os.getcwd() != test_directory:
        os.chdir(test_directory)

    print(F"TEST DIRECTORY: {test_directory}")

    print("Simple block preconditioner for 2D linear elasticity")

    # If the result directory already exists, wipe it, then make a new one
    if os.path.isdir(test_subdirectory + "/RESLT"):
        shutil.rmtree(test_subdirectory + "/RESLT")
    os.mkdir("RESLT")

    # Test info
    executable_name = "two_d_linear_elasticity_with_simple_block_diagonal_preconditioner"
    fpdiff_tol = 1.0e-12
    fpdiff_reltold = 0.1
    threshold_for_number_of_iterations = 3

    # Run the executable and redirect the output
    subprocess(f"../{executable_name} > RESLT/OUTPUT",
               shell=True)
    print("Done")

    # Output info to the log file and remove the indent from the multiline
    # string using dedent()
    print(dedent(f"""
    2D Linear elasticity simple preconditioner validation
    -----------------------------------------------------
    Validation directory:

        {test_directory}

    """), file=log_file)

    # Copy the files we need to inspect to the current directory
    shutil.copy(
        test_directory + "/RESLT/soln.dat",
        test_directory + "/linear_elasticity_simple_prec_results.dat")
    shutil.copy(
        test_directory + "/RESLT/iterative_solver_convergence.dat",
        test_directory + "/lin_elast_simple_iterative_solver_convergence.dat")

    # Run fpdiff on the files (if we can) to compare results
    if have_fpdiff:
        run_fpdiff(oomph_lib_root_path,
                   "linear_elasticity_simple_prec_results.dat",
                   "../validata/linear_elasticity_simple_prec_results.dat.gz",
                   fpdiff_tol, fpdiff_reltold, log_file)
    else:
        print(dedent("""
        dummy [OK] -- Can't run fpdiff.py because we don't have python or validata
        dummy [OK] -- Can't run fpdiff.py because we don't have python or validata
        """), file=log_file)

    # Compare number of iterations against reference data and append
    run_compare_file_length_with_tolerance(
        oomph_lib_root_path,
        "lin_elast_simple_iterative_solver_convergence.dat",
        "../validata/lin_elast_simple_iterative_solver_convergence.dat",
        threshold_for_number_of_iterations,
        log_file)

    # Rename the results directory
    os.rename("RESLT", "RESLT_linear_elasticity")
    return None


def get_validation_exit_status(validation_log_path: str, n_test: int) -> int:
    """Checks that we get the correct number of OKs. The function will return
    one of the following codes depending on the status of the tests:

            0 if all tests has passed.
            1 if some tests failed.
            2 if there are more 'OK' than expected.
    """
    # IMPLEMENT THIS
    EXIT_CODE = None

    with open("Validation/validation.log") as log_file:
        log_file_contents = log_file.read()

    # Use regular expression matching and list comprehension to create a list
    # of each location that "OK" occurs in the log file then find the length of
    # this list to find the number of occurrences
    n_match = len([m.start() for m in re.finditer('OK', log_file_contents)])

    # Output the test results to 'log_file' and set the exit code accordingly
    if n_match == n_test:
        output_string = f"""
        ======================================================================
        All tests in
            {os.getcwd()}
        passed successfully.
        ======================================================================
        """
        EXIT_CODE = 0

    elif n_match < n_test:
        output_string = f"""
        ======================================================================
        Only {n_match} of {n_test} test(s) passed; see
            {os.getcwd() + "/Validation/validation.log"}
        for details.
        ======================================================================
        """
        EXIT_CODE = 1

    else:
        output_string = f"""
        ======================================================================
        More OKs than tests! Need to update NUM_TESTS in
            {os.getcwd() + "/validate.sh"}
        ======================================================================
        """
        EXIT_CODE = 2

    # Print the multiline string to the log file. Rather than printing it
    # directly, we use this funky split-then-join method to make sure the string
    # isn't indented in the log file
    print(dedent(output_string), file=log_file)

    return EXIT_CODE


def main(argv) -> int:
    """Handle the execution of the demo driver tests and returns an exit code
    describing the outcome of the tests.
    """
    # Read in the oomph-lib root directory; we use this to find the location of
    # the fpdiff.py and validate_ok_count script. We also use this location to
    # store the master validation.log file.
    OOMPH_ROOT_DIR = argv[1]

    # Assess whether we have fpdiff.py
    HAVE_FPDIFF = True
    if (len(sys.argv) > 2) and (argv[2] == "no_fpdiff"):
        HAVE_FPDIFF = False

    # Set the number of tests to be checked
    NUM_TESTS = 89

    # Initialise the exit code value
    EXIT_CODE = 0

    # Get the full path to the subdirectory we'll do all of our testing in
    test_subdirectory = os.path.realpath("Validation")

    # Delete the Validation directory if it already exists then create a new one
    if os.path.isdir(test_subdirectory):
        shutil.rmtree(test_subdirectory)
    os.mkdir(test_subdirectory)

    # Enter the test directory
    os.chdir(test_subdirectory)

    # What is the path to the validata/ folder from here?
    path_to_validata = "../validata"

    # The file to output the contents of each test to
    log_file_name = "validation.log"

    # Use a Context Manager to safely handle opening/closing the log file for us
    with open(log_file_name, "w") as log_file:
        # Run the tests
        test_two_d_linear_elasticity_with_simple_block_diagonal_preconditioner(
            OOMPH_ROOT_DIR, test_subdirectory, path_to_validata, log_file, HAVE_FPDIFF)

    # Use os.path.join() to deduce paths in a platform-independent manner
    validation_log_path = os.path.join(test_subdirectory, log_file_name)

    # Get the exit code for this test from the log file
    return get_validation_exit_status(validation_log_path, NUM_TESTS)


# Execute main() if the script is run from the commandline and do nothing
# otherwise (e.g. if the script is imported by another script)
if __name__ == "__main__":
    sys.exit(main(sys.argv))
