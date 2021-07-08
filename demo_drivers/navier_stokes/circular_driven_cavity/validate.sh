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

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running circular driven cavity validation "
mkdir RESLT RESLT2
../circular_driven_cavity > OUTPUT_circular_driven_cavity
echo "done"
echo " " >> validation.log
echo "Circular driven cavity validation" >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat  RESLT/soln0.dat   RESLT/soln1.dat > results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results1.dat.gz  \
    results1.dat >> validation.log
fi

../circular_driven_cavity2 > OUTPUT_circular_driven_cavity2
echo "done"
echo " " >> validation.log
echo "Circular driven cavity 2 validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat  RESLT2/soln0.dat   RESLT2/soln1.dat > results2.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results2.dat.gz  \
    results2.dat >> validation.log
fi

# Append log to main validation log
cat validation.log >> ../../../../validation.log

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
