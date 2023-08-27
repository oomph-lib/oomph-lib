#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo hele shaw
#----------------------------
cd Validation

echo "Running string quad validation "
mkdir RESLT
../string_quad > OUTPUT
echo "done"
echo " " >> validation.log
echo "Integral structured rectangle validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > string_quad.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/string_quad.dat.gz   \
    string_quad.dat 0.1 1e-13  >> validation.log
fi

#----------------------------
echo "Running backward step quad validation "
mkdir RESLT
../backward_step_quad > OUTPUT
echo "done"
echo " " >> validation.log
echo "Integral structured rectangle validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > backward_step_quad.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/backward_step_quad.dat.gz   \
    backward_step_quad.dat 0.1 1e-13  >> validation.log
fi

#----------------------------
echo "Running simple cubic quad validation "
mkdir RESLT
../simple_cubic_quad > OUTPUT
echo "done"
echo " " >> validation.log
echo "Integral structured rectangle validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > simple_cubic_quad.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/simple_cubic_quad.dat.gz   \
    simple_cubic_quad.dat 0.1 1e-13  >> validation.log
fi

#----------------------------
echo "Running semi circle tri validation "
mkdir RESLT
../semi_circle_tri > OUTPUT
echo "done"
echo " " >> validation.log
echo "Integral structured rectangle validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > semi_circle_tri.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/semi_circle_tri.dat.gz   \
    semi_circle_tri.dat 0.1 1e-13  >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
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
