#!/bin/sh

set -o errexit
set -o nounset

# This tests that all the time steppers can run on a simple ODE example. So
# far it does NOT test that the adaptivity works.

# To test adaptivity: 1) run on a polynomial example that should be
# integrated exactly by all methods, if time steps are not enormous then
# the adaptivity has messed up. 2) Run some more complex example and
# roughly compare nsteps at the end, all methods should be within a factor
# of ~2 of each other.


# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=16

log_file="Validation/validation.log"

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Fixed time step tests 
# ============================================================

fixed_step_test()
{
    dir="Validation/fixed_step_$1_$2_$3"
    mkdir -p $dir

    # We fix the number of steps to 50 always so that it is easy to check
    # the errors against a file full of zeros with fpdiff.py.

    command="./basic_ode -max-step 50 -ts $1 -element-type $2 -ode $3 -outdir $dir"

    echo >> $log_file
    echo "Running $command" >> $log_file

    if $command > $dir/OUTPUT; then 
        # everything ran ok

        # Extract error norm (=|exact - approx|) to a file and check that
        # it's less than $4.
        cut -d\; -f 4 $dir/trace > $dir/error_norms
        fpdiff.py $dir/error_norms zeros 0.0 $4 >> $log_file
    else
        echo "FAIL: did not run successfully" >> $log_file
        echo >> $log_file
    fi
}

# Test each time stepper (except "real" imr) with normal ode elements
fixed_step_test bdf2 normal sin 4e-3
fixed_step_test imr normal sin 1e-3
fixed_step_test tr normal sin 1e-3

# Test each time stepper (including "real" imr) with ode elements modified
# to allow "real" imr.
fixed_step_test bdf2 imr_element sin 4e-3
fixed_step_test imr imr_element sin 1e-3
fixed_step_test real-imr imr_element sin 1e-3
fixed_step_test tr imr_element sin 1e-3

# Now run some other test odes, bdf2 is heavily used elsewhere so don't
# bother testing it as much here.
fixed_step_test imr normal poly2 1e-6
fixed_step_test tr normal poly2 1e-6
fixed_step_test real-imr imr_element poly2 1e-6

fixed_step_test imr normal stiff_test 1e-3
fixed_step_test tr normal stiff_test 1e-3
fixed_step_test real-imr imr_element stiff_test 1e-3

fixed_step_test imr normal damped_oscillation 5e-3
fixed_step_test tr normal damped_oscillation 5e-3
fixed_step_test real-imr imr_element damped_oscillation 5e-3


# Adaptive time step tests
# ============================================================

adaptive_step_test()
{
    dir="Validation/adaptive_step_$1_$2_$3"
    mkdir -p $dir

    # We fix the number of steps to 50 always so that it is easy to check
    # the errors against a file full of zeros with fpdiff.py.

    command="./basic_ode -max-step 50 -tol 1e-4 -ts $1 -element-type $2 -ode $3 -outdir $dir"
    echo "Running $command" >> $log_file

    if $command > $dir/OUTPUT; then 
        # everything ran ok

        # Extract error norm (=|exact - approx|) to a file and check that
        # it's less than $4.
        cut -d\; -f 4 $dir/trace > $dir/error_norms
        fpdiff.py $dir/error_norms zeros 0.0 $4 >> $log_file
    else
        echo "FAIL: did not run successfully" >> $log_file
        echo >> $log_file
    fi
}

# adaptive_step_test bdf2 normal poly2 1e-6

# ??ds fix these tests!

# Append local log file to global log file
cat $log_file >> ../../../validation.log



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
