#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation imposed boundary deformation using Lagrange multipliers
#-------------------------------------------------------------------

cd Validation

mkdir RESLT

echo "Running imposed boundary deformation using Lagrange multipliers "
../prescribed_displ_lagr_mult > OUTPUT_with_lagr_mult

echo "done"
echo " " >> validation.log
echo "Imposed boundary deformation using Lagrange multipliers" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/lagr2.dat RESLT/soln2.dat > result.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
fi

mv RESLT RESLT_with_lagr_mult

mkdir RESLT

echo "Running imposed boundary deformation without Lagrange multipliers "
../prescribed_displ_lagr_mult2 > OUTPUT_without_lagr_mult

echo "done"
echo " " >> validation.log
echo "Imposed boundary deformation without Lagrange multipliers" >> validation.log
echo "---------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln2.dat > result2.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result2.dat.gz \
    result2.dat  >> validation.log
fi

mv RESLT RESLT_without_lagr_mult


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

