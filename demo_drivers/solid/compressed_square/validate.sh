#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for compressed square
#---------------------------------

cd Validation



mkdir RESLT_norefine0
mkdir RESLT_norefine1
mkdir RESLT_norefine2
mkdir RESLT_norefine3
mkdir RESLT_norefine4
mkdir RESLT_norefine5
mkdir RESLT_norefine6
mkdir RESLT_norefine7
mkdir RESLT_norefine8
mkdir RESLT_norefine9
mkdir RESLT_norefine10
mkdir RESLT_norefine11
mkdir RESLT_norefine12
mkdir RESLT_norefine13
mkdir RESLT_norefine14



echo "Running compressed square "
../compressed_square > OUTPUT


echo "done"
echo " " >> validation.log
echo "compressed square" >> validation.log
echo "-----------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat \
RESLT_norefine0/soln1.dat \
RESLT_norefine1/soln1.dat \
RESLT_norefine2/soln1.dat \
RESLT_norefine3/soln1.dat \
RESLT_norefine4/soln1.dat \
RESLT_norefine5/soln1.dat \
RESLT_norefine6/soln1.dat \
RESLT_norefine7/soln1.dat \
RESLT_norefine8/soln1.dat \
RESLT_norefine9/soln1.dat \
RESLT_norefine10/soln1.dat \
RESLT_norefine11/soln1.dat \
RESLT_norefine12/soln1.dat \
RESLT_norefine13/soln1.dat \
RESLT_norefine14/soln1.dat \
    > result.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
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

