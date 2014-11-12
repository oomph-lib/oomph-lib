#!/bin/sh

set -o errexit
set -o nounset

# This tests that all the time steppers can run on a simple ODE example. So
# far it does NOT test: that the results are correct, adaptivity can run
# without crashing, adaptivity is appropriate.

# To test results: run with a few different -ode arguments, use cut to
# extract the solutions and exact solutions columns from the various trace
# files and use fpdiff.py to compare them.

# To test adaptivity: 1) run on a polynomial example that should be
# integrated exactly by all methods, if time steps are not enormous then
# the adaptivity has messed up. 2) Run some more complex example and
# roughly compare nsteps at the end, all methods should be within a factor
# of ~2 of each other.


# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=7

log_file="Validation/validation.log"

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

run_test()
{
    dir="Validation/$1_$2_$3"
    mkdir -p $dir

    command="./basic_ode -ts $1 -element-type $2 -ode $3 -outdir $dir"
    echo "Running $command" >> $log_file

    if $command > $dir/OUTPUT; then 
        # everything ran ok
        echo "OK" >> $log_file
        echo >> $log_file
    else
        echo "FAIL" >> $log_file
        echo >> $log_file
    fi
}

# Test each time stepper (except "real" imr) with normal ode elements
run_test bdf2 normal sin
run_test imr normal sin
run_test tr normal sin

# Test each time stepper (including "real" imr) with ode elements modified
# to allow "real" imr.
run_test bdf2 imr_element sin
run_test imr imr_element sin
run_test real-imr imr_element sin
run_test tr imr_element sin


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
