#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Validation for demo unsteady heat with and without ALE
#-------------------------------------------------------
cd Validation

echo "Running 2D unsteady heat validation with and without ALE"
mkdir RESLT
mkdir RESLT_ALE
../two_d_unsteady_heat  > OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady heat validation with and without ALE" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT_ALE/soln0.dat RESLT_ALE/soln1.dat \
    > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
fi



# Validation for demo adaptive unsteady heat with and without ALE
#-----------------------------------------------------------------

echo "Running 2D adaptive unsteady heat validation with and without ALE"
mkdir RESLT_adapt
mkdir RESLT_adapt_ALE
../two_d_unsteady_heat_adapt  > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "2D adaptive unsteady heat validation with and without ALE" >> validation.log
echo "---------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_adapt/soln0.dat RESLT_adapt/soln1.dat \
    RESLT_adapt_ALE/soln0.dat RESLT_adapt_ALE/soln1.dat \
    > result_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/result_adapt.dat.gz \
    result_adapt.dat  >> validation.log
fi






# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../../validation.log


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
 exit 0
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
   exit 1
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
  exit 2
  fi
fi
# Never get here
exit 10
