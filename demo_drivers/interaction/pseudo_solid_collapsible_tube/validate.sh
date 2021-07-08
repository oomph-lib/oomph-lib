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

echo "Running pseudo-solid collapsible tube validation "
mkdir RESLT
../pseudo_solid_collapsible_tube bla > OUTPUT
echo "done"
echo " " >> validation.log
echo "Pseudo-solid collapsible tube validation" >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/solid_soln3.dat  RESLT/solid_soln4.dat \
    RESLT/fluid_soln3.dat  RESLT/fluid_soln4.dat > results.dat
echo "Running pseudo-solid collapsible tube preconditioner validation "
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz \
    results.dat 6.0 3.0e-5 >> validation.log
fi
mv RESLT RESLT_WITHOUT_PREC
mkdir RESLT
../pseudo_solid_collapsible_tube bla bla > OUTPUT
echo "done"
echo " " >> validation.log
echo "Pseudo-solid collapsible tube preconditioner validation" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/solid_soln3.dat  RESLT/solid_soln4.dat \
    RESLT/fluid_soln3.dat  RESLT/fluid_soln4.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz \
    results.dat  6.0 3.0e-5 >> validation.log
fi
mv RESLT RESLT_WITH_PREC



#Append log to main validation log
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
