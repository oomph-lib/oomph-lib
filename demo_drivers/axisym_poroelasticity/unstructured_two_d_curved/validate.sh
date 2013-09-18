#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=1


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for unstructured two-d curved axisymmetric poroelasticity
#----------------------------------------------------------------------
cd Validation

if [ 0 -eq 0 ]; then

echo "Running unstructured two-d curved axisymmetric poroelasticity validation "
mkdir RESLT
../unstructured_two_d_curved --validation --dt 0.01 --n_steps 10 > OUTPUT_unstructured_two_d_curved
echo "done"
echo " " >> validation.log
echo "Unstructured two-d curved axisymmetric poroelasticity validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
touch unstructured_two_d_curved.dat
rm -f unstructured_two_d_curved.dat
cat RESLT/trace.dat | tail -n 1 | cut -d ' ' -f 4-6 >> unstructured_two_d_curved.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../../bin/fpdiff.py ../validata/unstructured_two_d_curved.dat.gz   \
  unstructured_two_d_curved.dat  >> validation.log
fi

fi

if [ 1 -eq 0 ]; then

echo "Running adaptive unstructured two-d curved axisymmetric poroelasticity validation "
mkdir RESLT
../adaptive_unstructured_two_d_curved --validation --dt 0.01 --n_steps 10 --max_adapt 0 > OUTPUT_adaptive_unstructured_two_d_curved
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured two-d curved axisymmetric poroelasticity validation" >> validation.log
echo "-------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
touch adaptive_unstructured_two_d_curved.dat
rm -f adaptive_unstructured_two_d_curved.dat
cat RESLT/trace.dat | tail -n 1 | cut -d ' ' -f 4-6 >> adaptive_unstructured_two_d_curved.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../../bin/fpdiff.py ../validata/unstructured_two_d_curved.dat.gz   \
  adaptive_unstructured_two_d_curved.dat  >> validation.log
fi

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
