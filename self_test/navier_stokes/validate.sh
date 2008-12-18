#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running circular driven cavity validation "
mkdir RESLT
../circular_driven_cavity > OUTPUT_circular_driven_cavity
echo "done"
echo " " >> validation.log
echo "Circular driven cavity validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT
cat  soln0.dat   soln11.dat  soln14.dat  soln17.dat  soln2.dat  soln5.dat  \
     soln8.dat soln1.dat   soln12.dat  soln15.dat  soln18.dat  soln3.dat  \
     soln6.dat  soln9.dat soln10.dat  soln13.dat  soln16.dat  soln19.dat  \
     soln4.dat  soln7.dat > circular_driven_cavity_results.dat
cd ..



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/circular_driven_cavity_validation.dat.gz  \
         RESLT/circular_driven_cavity_results.dat >> validation.log
fi

mv RESLT RESLT_circular_driven_cavity



# Append log to main validation log
cat validation.log >> ../../../validation.log

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
