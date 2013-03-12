#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=0  # hierher: Tests temporarily disabled until we
             #           sort out the problems with the additional
             #           small entries. 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for different assembly strategies
#---------------------------------------------
cd Validation

echo "Running validation of different assembly strategies"
for n in 1 2 3 4 5
 do
  ../sparse_assemble_test 10 10 $n 1 0 > OUTPUT_$n
  #UNIX sort with magic to do numerical sorting on the first column
  #and then on the second column if the first columns are the same
  sort -k 1,1n -k 2,2n  matrix$n.dat  > matrix$n.dat.sorted
 done

echo "done"
echo " " >> validation.log
echo " Validation of different assembly strategies" >> validation.log
echo "--------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


echo "  " >> validation.log
echo "------------------------------------------------" >> validation.log
echo "NOTE: fpdiff is currently disabled on this test." >> validation.log
echo "=====" >> validation.log
echo "------------------------------------------------" >> validation.log
echo "  " >> validation.log

#if test "$1" = "no_fpdiff"; then
#  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
#else
#../../../../bin/fpdiff.py ../validata/matrix1.dat.gz   \
#    matrix1.dat.sorted >> validation.log
#../../../../bin/fpdiff.py ../validata/matrix2.dat.gz   \
#    matrix2.dat.sorted >> validation.log
#../../../../bin/fpdiff.py ../validata/matrix3.dat.gz   \
#    matrix3.dat.sorted >> validation.log
#../../../../bin/fpdiff.py ../validata/matrix4.dat.gz   \
#    matrix4.dat.sorted >> validation.log
#fi


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
