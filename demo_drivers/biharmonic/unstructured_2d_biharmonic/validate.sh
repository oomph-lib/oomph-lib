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

# Validation for unstructured mesh 2d biharmonic with bell elements
#-----------------------------------------
cd Validation

echo "Running two d unstructured triangular biharmonic bell validation "
mkdir RESLT_bell
../unstructured_2d_biharmonic_bellelement ../Circle1.1.node ../Circle1.1.ele ../Circle1.1.poly > OUTPUT_biharmonic_bell
echo "done"
echo " " >> validation.log
echo "Two d unstructured triangular biharmonic bell validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/biharmonicbell.dat.gz  \
         RESLT_bell/soln0.dat 0.1 1.0e-11 >> validation.log
fi

# Validation for unstructured mesh 2d biharmonic with curved elements
#-----------------------------------------

echo "Running two d unstructured triangular biharmonic curved validation "
mkdir RESLT_curved
../unstructured_2d_biharmonic_curvedelement ../Circle1.1.node ../Circle1.1.ele ../Circle1.1.poly > OUTPUT_biharmonic_curved
echo "done"
echo " " >> validation.log
echo "Two d unstructured triangular biharmonic curved validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/biharmoniccurved.dat.gz  \
         RESLT_curved/soln0.dat 0.1 3.0e-10 >> validation.log
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
