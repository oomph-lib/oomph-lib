#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation
mkdir RESLT
var="../structured_no_correction > OUTPUT_no_correction"
echo $var
eval $var
echo "done"
cd ../
LOG="Validation/validation.log"
echo " " >> $LOG 
echo "Validation run" >> $LOG
echo "---------------------------------------------" >> $LOG
echo " " >> $LOG
echo "Validation directory: " >> $LOG
echo " " >> $LOG
echo "  " `pwd` >> $LOG
echo " " >> $LOG
cat Validation/RESLT/slip_surface0.csv > Validation/structured_no_correction.dat
mv Validation/RESLT/soln0.dat Validation/soln0.dat
rm -r Validation/RESLT
if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> $LOG
else
    ../../../bin/fpdiff.py validata/structured_no_correction.dat.gz  \
       Validation/structured_no_correction.dat >> $LOG
fi

cd Validation
mkdir -p RESLT
var="../structured_with_correction > OUTPUT_with_correction"
echo $var
eval $var
echo "done"
cd ../
LOG="Validation/validation.log"
echo " " >> $LOG 
echo "Validation run" >> $LOG
echo "---------------------------------------------" >> $LOG
echo " " >> $LOG
echo "Validation directory: " >> $LOG
echo " " >> $LOG
echo "  " `pwd` >> $LOG
echo " " >> $LOG
cat Validation/RESLT/slip_surface0.csv > Validation/structured_with_correction.dat
mv Validation/RESLT/soln0.dat Validation/soln1.dat
rm -r Validation/RESLT
if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> $LOG
else
    ../../../bin/fpdiff.py validata/structured_with_correction.dat.gz  \
       Validation/structured_with_correction.dat >> $LOG
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
