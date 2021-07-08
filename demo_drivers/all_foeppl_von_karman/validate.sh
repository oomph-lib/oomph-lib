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

cd Validation

# Validation for circular disk Foeppl von Karman
#-----------------------------------------------

echo "Running circular disk Foeppl von Karman validation "
mkdir RESLT
../circular_disk --validation > OUTPUT_circular_disk
echo "done"
echo " " >> validation.log
echo "Circular disk Foeppl von Karman validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > circular_disk_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../bin/fpdiff.py ../validata/circular_disk_results.dat.gz   \
  circular_disk_results.dat  >> validation.log
fi
mv RESLT RESLT_circular_disk



# Validation for displacement-based circular disk Foeppl von Karman
#------------------------------------------------------------------

echo "Running displacement-based circular disk Foeppl von Karman validation "
mkdir RESLT
../displacement_based_circular_disk --validation > OUTPUT_displacement_based_circular_disk
echo "done"
echo " " >> validation.log
echo "Displacement-based circular disk Foeppl von Karman validation" >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > displacement_based_circular_disk_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../bin/fpdiff.py ../validata/displacement_based_circular_disk_results.dat.gz   \
  displacement_based_circular_disk_results.dat  >> validation.log
fi
mv RESLT RESLT_displacement_based_circular_disk



echo "Running axisym Foeppl von Karman validation "
mkdir RESLT
../axisym_fvk > OUTPUT_axisym_fvk
echo "done"
echo " " >> validation.log
echo "Axisym Foeppl von Karman validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/sol_9.dat RESLT/w_centre.dat > axisym_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../bin/fpdiff.py ../validata/axisym_results.dat.gz   \
  axisym_results.dat  >> validation.log
fi
mv RESLT RESLT_axisym_fvk





echo "Running displacement based axisym Foeppl von Karman validation "
mkdir RESLT
../axisym_displ_based_fvk > OUTPUT_displacement_based_axisym_fvk
echo "done"
echo " " >> validation.log
echo "Displacement based axisym Foeppl von Karman validation" >> validation.log
echo "------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln9.dat RESLT/w_centre.dat > displacement_based_axisym_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../bin/fpdiff.py ../validata/displacement_based_axisym_results.dat.gz   \
  displacement_based_axisym_results.dat  >> validation.log
fi
mv RESLT RESLT_displacement_based_axisym_fvk



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
