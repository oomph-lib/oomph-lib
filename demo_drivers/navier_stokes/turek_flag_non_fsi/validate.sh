#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Validation for non-FSI Turek problem (algebraic node update)
#-------------------------------------------------------------
cd Validation

echo "Running non-FSI Turek problem (algebraic node update)"
mkdir RESLT

# Do validation run
../turek_flag_non_fsi_alg blabla > OUTPUT_alg
echo "done"
echo " " >> validation.log
echo "Non-FSI Turek problem (algebraic node update)" \
>> validation.log
echo "---------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat \
    > result_alg.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_alg.dat.gz \
    result_alg.dat >> validation.log
fi


mv RESLT RESLT_alg


# Validation for non-FSI Turek problem (Domain/MacroElement-based node update)
#-----------------------------------------------------------------------------

echo "Running non-FSI Turek problem (Domain-based node update)"
mkdir RESLT

# Do validation run
../turek_flag_non_fsi_macro blabla > OUTPUT_macro
echo "done"
echo " " >> validation.log
echo "Non-FSI Turek problem (Domain/MacroElement-based node update)" \
>> validation.log
echo "-------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat \
    > result_macro.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_macro.dat.gz \
    result_macro.dat >> validation.log
fi


mv RESLT RESLT_macro


# Append log to main validation log
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
