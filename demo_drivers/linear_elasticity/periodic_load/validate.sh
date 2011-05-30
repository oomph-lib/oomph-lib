#! /bin/sh


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
         ./adapt.dat 0.1 2.0e-12 >> validation.log
fi


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
