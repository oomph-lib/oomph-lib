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

# Validation for rectangular driven cavity
#-----------------------------------------
cd Validation

echo "Running two-fluid spherical cap validation "
mkdir RESLT
cd RESLT
../../two_fluid_spherical_cap > ../OUTPUT_two_fluid_spherical_cap
cd ..
echo "done"
echo " " >> validation.log
echo "Two-fluid sperical cap validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/step0_spine.dat  RESLT/step1_spine.dat  RESLT/step2_spine.dat \
     RESLT/step3_spine.dat  RESLT/step4_spine.dat  RESLT/step5_spine.dat \
     > two_fluid_cap.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/cap.dat.gz  \
   two_fluid_cap.dat 0.1 2.5e-8 >> validation.log
fi

mv RESLT RESLT_two_fluid_spherical_cap

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
