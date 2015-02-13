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

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running BlockSelector validation "
../block_selector_test >> OUTPUT_WARNING_MESSAGE 2>&1
echo "done"
echo " " >> validation.log
echo "BlockSelector validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat OUTPUT OUTPUT_WARNING_MESSAGE >> blockselectoroutput
rm -rf OUTPUT OUTPUT_WARNING_MESSAGE
DIFF=$(diff blockselectoroutput ./../validata/block_selector_test_data.dat)
if [ "$DIFF" != "" ] 
then
    echo "Running diff on " >> validation.log
    echo "blockselectoroutput ./../validata/block_selector_test_data.dat" >> validation.log
    diff blockselectoroutput ./../validata/block_selector_test_data.dat >> validation.log
    echo "BlockSelector validation [FAILED]" >> validation.log
else
    echo "BlockSelector validation [OK]" >> validation.log
fi

# Append log to main validation log
cat validation.log >> ../../../validation.log

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
