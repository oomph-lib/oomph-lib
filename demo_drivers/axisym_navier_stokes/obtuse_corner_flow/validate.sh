#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=8

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation
cp default_parameters.dat Validation/parameters.dat

test_script()
{
    # Set up log file
    LOG="Validation/validation.log"
    echo " " >> $LOG 
    echo "Validation run" >> $LOG
    echo "---------------------------------------------" >> $LOG
    echo " " >> $LOG
    echo "Validation directory: " >> $LOG
    echo " " >> $LOG
    echo "  " `pwd` >> $LOG
    echo " " >> $LOG

    # If we are on mac then we can't run the tests yet
    if [[ "$(uname)" == "Darwin" ]]; then
        echo "dummy [OK] -- macOS is not supported for this demo yet" >> $LOG
    else
        # Run the test
        cd Validation
        mkdir $1
        var="../$2 > $2.out"
        echo $var
        eval $var
        echo "done"
        cd ../
        # Sorting here as the MPI runs have a different mesh ordering for some reason...
        sort Validation/$1/slip_surface0.dat > Validation/$2.dat
        if test "$1" = "no_fpdiff"; then
            echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> $LOG
        else
            ../../../bin/fpdiff.py validata/$2.dat.gz  \
               Validation/$2.dat 0.1 2e-14 >> $LOG
        fi
    fi
}

# This is needed to ensure that sort is locale independent
export LC_ALL=C

test_script RESLT_axi_no_fix structured_no_correction 
test_script RESLT_axi_no_fix_region structured_no_correction_region
test_script RESLT_axi_fix structured_with_correction 
test_script RESLT_axi_fix_region structured_with_correction_region
test_script RESLT_axi_sprittles_region sprittles_region
test_script RESLT_axi_no_fix_unstr unstructured_no_correction 
test_script RESLT_axi_fix_unstr unstructured_with_correction

LOG="Validation/validation.log"
echo " " >> $LOG 
echo "Validation run" >> $LOG
echo "---------------------------------------------" >> $LOG
echo " " >> $LOG
echo "Validation directory: " >> $LOG
echo " " >> $LOG
echo "  " `pwd` >> $LOG
echo " " >> $LOG
echo "./test_parameters"
if ./test_parameters -eq 0; then
    echo "done"
    echo "[OK] -- test_parameters passed" >> $LOG
else
    echo "done"
    echo "[FAILED] -- test_parameters failed" >> $LOG
fi

#######################################################################

# Append log to main validation log
cat $LOG >> ../../../validation.log

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
