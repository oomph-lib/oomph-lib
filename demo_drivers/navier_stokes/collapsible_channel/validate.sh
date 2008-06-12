#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=6

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Validation for collapsible channel with prescribed wall motion
#---------------------------------------------------------------
cd Validation

echo "Running validation for collapsible channel with prescribed wall motion "
rm -rf RESLT_no_bl_squash
rm -rf RESLT
mkdir RESLT
# Do validation run
../collapsible_channel blabla > OUTPUT 
mv RESLT RESLT_no_bl_squash
mv OUTPUT OUTPUT_no_bl_squash
echo "done"
echo " " >> validation.log
echo "Collapsible channel with prescribed wall motion validation" \
>> validation.log
echo "----------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT_no_bl_squash/soln0.dat RESLT_no_bl_squash/soln1.dat  \
    RESLT_no_bl_squash/soln2.dat RESLT_no_bl_squash/soln3.dat  \
    > result_no_bl_squash.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/no_bl_squash.dat.gz \
    result_no_bl_squash.dat 0.1 1.0e-8 >> validation.log
fi



# Validation for collapsible channel with prescribed wall motion (BL squash)
#---------------------------------------------------------------------------

echo "Running validation for collapsible channel with prescribed  "
echo "wall motion with BL squashing in mesh"
rm -rf RESLT_bl_squash
rm -rf RESLT
mkdir RESLT
# Do validation run
../collapsible_channel_bl_squash blabla > OUTPUT 
mv RESLT RESLT_bl_squash
mv OUTPUT OUTPUT_bl_squash
echo "done"
echo " " >> validation.log
echo "Collapsible channel with prescribed wall motion validation (BL squash)" \
>> validation.log
echo "----------------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT_bl_squash/soln0.dat RESLT_bl_squash/soln1.dat  \
    RESLT_bl_squash/soln2.dat RESLT_bl_squash/soln3.dat  \
    > result_bl_squash.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/bl_squash.dat.gz \
    result_bl_squash.dat 0.1 1.0e-8 >> validation.log
fi






# Validation for alg collapsible channel with prescribed wall motion
#---------------------------------------------------------------
echo "Running validation for alg collapsible channel with prescribed wall motion "
rm -rf RESLT_alg_no_bl_squash
rm -rf RESLT
mkdir RESLT
# Do validation run
../collapsible_channel_algebraic blabla > OUTPUT 
mv RESLT RESLT_alg_no_bl_squash
mv OUTPUT OUTPUT_alg_no_bl_squash
echo "done"
echo " " >> validation.log
echo "Alg collapsible channel with prescribed wall motion validation" \
>> validation.log
echo "----------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT_alg_no_bl_squash/soln0.dat RESLT_alg_no_bl_squash/soln1.dat  \
    RESLT_alg_no_bl_squash/soln2.dat RESLT_alg_no_bl_squash/soln3.dat  \
    > result_alg_no_bl_squash.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg_no_bl_squash.dat.gz \
    result_alg_no_bl_squash.dat 0.1 1.0e-8 >> validation.log
fi



# Validation for alg collapsible channel with prescribed wall motion (BL squash)
#---------------------------------------------------------------------------

echo "Running validation for alg collapsible channel with prescribed  "
echo "wall motion with BL squashing in mesh"
rm -rf RESLT_alg_bl_squash
rm -rf RESLT
mkdir RESLT
# Do validation run
../collapsible_channel_algebraic_bl_squash blabla > OUTPUT 
mv RESLT RESLT_alg_bl_squash
mv OUTPUT OUTPUT_alg_bl_squash
echo "done"
echo " " >> validation.log
echo "Alg collapsible channel with prescribed wall motion validation (BL squash)" \
>> validation.log
echo "----------------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT_alg_bl_squash/soln0.dat RESLT_alg_bl_squash/soln1.dat  \
    RESLT_alg_bl_squash/soln2.dat RESLT_alg_bl_squash/soln3.dat  \
    > result_alg_bl_squash.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg_bl_squash.dat.gz \
    result_alg_bl_squash.dat 0.1 1.0e-8 >> validation.log
fi






# Validation for adapt alg collapsible channel with prescribed wall motion
#---------------------------------------------------------------
echo "Running validation for adapt alg collapsible channel with prescribed wall motion "
rm -rf RESLT_adapt_alg_no_bl_squash
rm -rf RESLT
mkdir RESLT
# Do validation run
../collapsible_channel_adaptive_algebraic blabla > OUTPUT 
mv RESLT RESLT_adapt_alg_no_bl_squash
mv OUTPUT OUTPUT_adapt_alg_no_bl_squash
echo "done"
echo " " >> validation.log
echo "Adaptive alg collapsible channel with prescribed wall motion validation" \
>> validation.log
echo "----------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT_adapt_alg_no_bl_squash/soln0.dat RESLT_adapt_alg_no_bl_squash/soln1.dat  \
    RESLT_adapt_alg_no_bl_squash/soln2.dat RESLT_adapt_alg_no_bl_squash/soln3.dat  \
    > result_adapt_alg_no_bl_squash.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adapt_alg_no_bl_squash.dat.gz \
    result_adapt_alg_no_bl_squash.dat 0.1 1.0e-8 >> validation.log
fi



# Validation for adapt alg collapsible channel with prescribed wall motion (BL squash)
#---------------------------------------------------------------------------

echo "Running validation for adaptive alg collapsible channel with prescribed  "
echo "wall motion with BL squashing in mesh"
rm -rf RESLT_adapt_alg_bl_squash
rm -rf RESLT
mkdir RESLT
# Do validation run
../collapsible_channel_adaptive_algebraic_bl_squash blabla > OUTPUT 
mv RESLT RESLT_adapt_alg_bl_squash
mv OUTPUT OUTPUT_adapt_alg_bl_squash
echo "done"
echo " " >> validation.log
echo "Adaptive alg collapsible channel with prescribed wall motion validation (BL squash)" \
>> validation.log
echo "----------------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT_adapt_alg_bl_squash/soln0.dat RESLT_adapt_alg_bl_squash/soln1.dat  \
    RESLT_adapt_alg_bl_squash/soln2.dat RESLT_adapt_alg_bl_squash/soln3.dat  \
    > result_adapt_alg_bl_squash.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adapt_alg_bl_squash.dat.gz \
    result_adapt_alg_bl_squash.dat 0.1 1.0e-8 >> validation.log
fi






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
