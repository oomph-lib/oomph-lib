#!/bin/bash

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

# Set the number of tests to be checked, computed as we go along.
NUM_TESTS=0

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

        # Extract error norm (=|exact - approx|), delete header line and
        # write to a file. Then check that all norms are less than $4.
        cut -d\; -f 4 $dir/trace | sed '1d' > $dir/error_norms
        ../../../bin/fpdiff.py $dir/error_norms zeros 0.0 $4 >> $log_file
    else	
        echo "FAIL: did not run successfully" >> $log_file
        echo >> $log_file
    fi
    
}

# Test each time stepper (except "real" imr) with normal ode elements
fixed_step_test bdf2 normal sin 4e-3 # bdf2 is less accurate
fixed_step_test imr normal sin 1e-3
fixed_step_test tr normal sin 1e-3
NUM_TESTS=$(expr $NUM_TESTS + 3)

# Test each time stepper (including "real" imr) with ode elements modified
# to allow "real" imr.
fixed_step_test bdf2 imr-element sin 4e-3 
fixed_step_test imr imr-element sin 1e-3
fixed_step_test real-imr imr-element sin 1e-3
fixed_step_test tr imr-element sin 1e-3
NUM_TESTS=$(expr $NUM_TESTS + 4)

# Now run some other test odes, bdf2 is heavily used elsewhere so don't
# bother testing it as much here.
fixed_step_test imr normal poly2 1e-6
fixed_step_test tr normal poly2 1e-6
fixed_step_test real-imr imr-element poly2 1e-6

fixed_step_test imr normal stiff_test 1e-3
fixed_step_test tr normal stiff_test 1e-3
fixed_step_test real-imr imr-element stiff_test 1e-3

fixed_step_test imr normal damped_oscillation 5e-3
fixed_step_test tr normal damped_oscillation 5e-3
fixed_step_test real-imr imr-element damped_oscillation 5e-3

NUM_TESTS=$(expr $NUM_TESTS + 9)


# Adaptive time step tests
# ============================================================

adaptive_step_test()
{
    dir="Validation/adaptive_step_$1_$2_$3"
    mkdir -p $dir

    # We fix the number of steps to 50 always so that it is easy to check
    # the errors against a file full of zeros with fpdiff.py.

    command="./basic_ode -max-step 100 -tmax $6 -tol 1e-5 -ts $1 -element-type $2 -ode $3 -outdir $dir"
    echo "Running $command" >> $log_file

    if $command > $dir/OUTPUT; then 
        # everything ran ok

        # Extract error norm (=|exact - approx|), delete header line and
        # write to a file. Then check that all norms are less than $4.
        cut -d\; -f 4 $dir/trace | sed '1d' > $dir/error_norms
        if ../../../bin/fpsmall.py $dir/error_norms $4 > $dir/fpsmall_output; then
            nsteps="$(wc -l < $dir/trace)" # Have to use '<' or we will get
                                           # filename in the output of wc.
            if [[ $nsteps -lt "$5" ]]; then
                echo "  [OK]" >> $log_file
            else
                echo "  FAIL: too many steps" >> $log_file
            fi
        else
            echo "  FAIL: errors too large" >> $log_file
            cat $dir/fpsmall_output >> $log_file
        fi
    else
        echo "  FAIL: did not run successfully" >> $log_file
        echo >> $log_file
    fi
}

# Test that the adaptivity recognises that these methods exactly integrate
# a second order polynomial (and so only need a few steps). Also tests that
# adaptive versions can run without crashing.
adaptive_step_test bdf2 normal poly2 1e-6 10 20.0
adaptive_step_test tr normal poly2 1e-6 10 20.0
adaptive_step_test imr normal poly2 1e-6 10 20.0
adaptive_step_test real-imr imr-element poly2 1e-6 10 20.0
NUM_TESTS=$(expr $NUM_TESTS + 4)

# Test adaptivity for a stiff problem
adaptive_step_test bdf2 normal stiff_test 2e-2 20 5.0 # bdf2 is less accurate again
adaptive_step_test tr normal stiff_test 1.1e-2 20 5.0
adaptive_step_test imr normal stiff_test 1.1e-2 20 5.0
adaptive_step_test real-imr imr-element stiff_test 1.1e-2 20 5.0
NUM_TESTS=$(expr $NUM_TESTS + 4)

# A harder adaptivity test with a "real" ode
adaptive_step_test bdf2 normal damped_oscillation 1e-2 55 1.0 # bdf2 is worse again...
adaptive_step_test tr normal damped_oscillation 1e-2 40 1.0
adaptive_step_test imr normal damped_oscillation 1e-2 40 1.0
adaptive_step_test real-imr imr-element damped_oscillation 1e-2 40 1.0
NUM_TESTS=$(expr $NUM_TESTS + 4)


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
