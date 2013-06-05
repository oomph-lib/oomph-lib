#! /bin/sh

set -o errexit
no_fpdiff=false # `test "$1" = "no_fpdiff"`
set -o nounset

#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo unsteady heat
#----------------------------------
cd Validation

echo "Running 2D unsteady heat midpoint elements with bdf2 time stepper validation"
mkdir RESLT
../two_d_unsteady_heat_midpoint "bdf2" > OUTPUT_bdf2
echo "done"
echo >> validation.log
echo "2D unsteady heat midpoint elements with bdf2 time stepper validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo >> validation.log
echo "Validation directory: " >> validation.log
echo >> validation.log
echo "  " `pwd` >> validation.log
echo >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_bdf2.dat
mv RESLT RESLT_bdf2

if $no_fpdiff; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result_bdf2.dat  >> validation.log
fi


echo "Running unsteady heat midpoint elements with midpoint time stepper validation"
mkdir RESLT
../two_d_unsteady_heat_midpoint "midpoint" > OUTPUT_midpoint
echo "done"
echo >> validation.log
echo "2D unsteady heat midpoint elements with midpoint time stepper validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo >> validation.log
echo "Validation directory: " >> validation.log
echo >> validation.log
echo "  " `pwd` >> validation.log
echo >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_midpoint.dat
mv RESLT RESLT_midpoint

if $no_fpdiff; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result_midpoint.dat  >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' Validation/validation.log`
if  [ $OK_COUNT -eq $NUM_TESTS ]; then
 echo 
 echo "======================================================================"
 echo  
 echo "All tests in" 
 echo  
 echo "    `pwd`    "
 echo 
 echo "passed successfully."
 echo 
 echo "======================================================================"
 echo  
 exit 0
else
  if [ $OK_COUNT -lt $NUM_TESTS ]; then
   echo 
   echo "======================================================================"
   echo  
   echo "Only $OK_COUNT of $NUM_TESTS test(s) passed; see"
   echo  
   echo "    `pwd`/Validation/validation.log"
   echo  
   echo "for details" 
   echo  
   echo "======================================================================"
   echo 
   exit 1
  else 
   echo 
   echo "======================================================================"
   echo  
   echo "More OKs than tests! Need to update NUM_TESTS in"
   echo  
   echo "    `pwd`/validate.sh"
   echo 
   echo "======================================================================"
   echo 
  exit 2
  fi
fi
# Never get here
exit 10
