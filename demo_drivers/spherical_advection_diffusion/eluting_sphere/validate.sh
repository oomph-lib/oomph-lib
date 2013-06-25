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

# Validation for adaptive spin up
#--------------------------------
cd Validation
mkdir RESLT
mkdir RESLT_adapt

echo "Running spherical polar eluting sphere validation "

../eluting_sphere > OUTPUT_eluting_sphere
echo "done"
echo " " >> validation.log
echo "Spherical polar eluting sphere validation" >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln10.dat > el_sphere.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/el_sphere.dat.gz  \
         el_sphere.dat 0.1 1.0e-14 >> validation.log
fi


echo "Running adaptive spherical polar eluting spehre validation "
../eluting_sphere_adapt > OUTPUT_eluting_sphere_adapt
echo "done"
echo " " >> validation.log
echo "Adaptive spherical polar eluting sphere validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_adapt/soln10.dat > el_sphere_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/el_sphere_adapt.dat.gz  \
         el_sphere_adapt.dat 0.1 1.0e-14 >> validation.log
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
