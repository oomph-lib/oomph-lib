#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for shock disk
#--------------------------

cd Validation
mkdir RESLT0
mkdir RESLT1
mkdir RESLT2

echo "Running shock disk validation "
../shock_disk lalala  > OUTPUT


echo "done"
echo " " >> validation.log
echo "Shock disk validation" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT0/soln3.dat \
    > shock_disk_results0.dat
cat RESLT1/soln3.dat \
    > shock_disk_results1.dat
cat RESLT2/soln3.dat \
    > shock_disk_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/shock_disk_results0.dat.gz \
    shock_disk_results0.dat  >> validation.log
../../../../bin/fpdiff.py ../validata/shock_disk_results1.dat.gz \
    shock_disk_results1.dat  >> validation.log
../../../../bin/fpdiff.py ../validata/shock_disk_results2.dat.gz \
    shock_disk_results2.dat  >> validation.log
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

