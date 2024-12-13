#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Elastic bretherton problem
#-------------------------------------
cd Validation

if [[ "$(uname)" == "Darwin" ]]; then
    echo " " >> validation.log
    echo "Elastic Bretherton validation" >> validation.log
    echo "------------------------" >> validation.log
    echo " " >> validation.log
    echo "Validation directory: " >> validation.log
    echo " " >> validation.log
    echo "  " `pwd` >> validation.log
    echo " " >> validation.log
    echo "dummy [OK] -- Not running on macOS as there is a long-standing bug in this test which fails non-deterministically." >> validation.log
    # Append log to main validation log
    cat validation.log >> ../../../../validation.log
else
    echo "Running elastic Bretherton problem validation "
    mkdir RESLT
    ../elastic_bretherton lalala > OUTPUT
    echo "done"
    echo " " >> validation.log
    echo "Elastic Bretherton validation" >> validation.log
    echo "------------------------" >> validation.log
    echo " " >> validation.log
    echo "Validation directory: " >> validation.log
    echo " " >> validation.log
    echo "  " `pwd` >> validation.log
    echo " " >> validation.log
    cat  RESLT/trace.dat  RESLT/soln2.dat > el_breth.dat

    if test "$1" = "no_fpdiff"; then
      echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    else
    ../../../../bin/fpdiff.py ../validata/el_breth.dat.gz  \
             el_breth.dat  0.5 9.0e-9 >> validation.log
    fi

    # Append log to main validation log
    cat validation.log >> ../../../../validation.log
fi

cd ..

#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
