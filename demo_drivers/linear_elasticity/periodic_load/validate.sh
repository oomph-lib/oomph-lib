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

# Validation for periodically loaded linear elastic solid
#---------------------------------------------------------
cd Validation

echo "Running periodic load on linear elastic solid"
mkdir RESLT 
../periodic_load  > OUTPUT_periodic_load
mv RESLT RESLT_periodic_load
echo "done"
echo "Running periodic load on linear elastic solid (adaptive mesh)"
mkdir RESLT
../refineable_periodic_load  > OUTPUT_refineable_periodic_load
mv RESLT RESLT_refineable_periodic_load
echo "done"
echo " " >> validation.log
echo "Periodic load on linearly elastic solid" >> validation.log
echo "-------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_periodic_load/soln.dat > ./period.dat
cat  RESLT_refineable_periodic_load/soln.dat > ./adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Fixed mesh" >> validation.log
../../../../bin/fpdiff.py ../validata/period.dat.gz  \
         ./period.dat 0.1 2.0e-12 >> validation.log
echo "Adaptive mesh" >> validation.log
../../../../bin/fpdiff.py ../validata/adapt.dat.gz  \
         ./adapt.dat 0.1 5.0e-9 >> validation.log
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
