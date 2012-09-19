#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for non-adaptive scattering
#-----------------------------------
cd Validation

echo "Running unstructured PML scattering"
mkdir RESLT
../unstructured_two_d_helmholtz > OUTPUT_non_adapt
echo "done"
echo " " >> validation.log
echo "Running unstructured PML scattering" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > pml_scattering_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pml_scattering_results.dat.gz   \
    pml_scattering_results.dat  0.9 1.0e-14 >> validation.log
fi
mv RESLT RESLT_non_adapt


# Validation for adaptive scattering
#-----------------------------------

echo "Running adaptive unstructured PML scattering"
mkdir RESLT
../unstructured_two_d_helmholtz_adapt > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "Running adaptive unstructured PML scattering" >> validation.log
echo "--------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > pml_adaptive_scattering_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pml_scattering_results.dat.gz   \
    pml_scattering_results.dat  0.9 1.0e-14 >> validation.log
fi
mv RESLT RESLT_adapt




# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..

#######################################################################


#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' Validation/validation.log`
if  [ $OK_COUNT -eq $NUM_TESTS ]; then
 echo " "
 echo "======================================================================"
 echo " " 
 echo "All tests in" 
 echo " " 
 echo "    `pwd`    "
 echo " "
 echo "passed successfully."
 echo " "
 echo "======================================================================"
 echo " " 
else
  if [ $OK_COUNT -lt $NUM_TESTS ]; then
   echo " "
   echo "======================================================================"
   echo " " 
   echo "Only $OK_COUNT of $NUM_TESTS test(s) passed; see"
   echo " " 
   echo "    `pwd`/Validation/validation.log"
   echo " " 
   echo "for details" 
   echo " " 
   echo "======================================================================"
   echo " "
  else 
   echo " "
   echo "======================================================================"
   echo " " 
   echo "More OKs than tests! Need to update NUM_TESTS in"
   echo " " 
   echo "    `pwd`/validate.sh"
   echo " "
   echo "======================================================================"
   echo " "
  fi
fi
