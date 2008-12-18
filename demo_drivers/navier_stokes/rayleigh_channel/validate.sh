#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running rayleigh_channel"
mkdir RESLT_CR RESLT_TH
../rayleigh_channel 0 0 > OUTPUT_rayleigh_channel1
echo "done (start from periodic solution)"
echo " " >> validation.log
echo "Rayleigh_channel validation with IC=periodic soln" >> validation.log
echo "-------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT_CR
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results1_CR.dat
cd ../RESLT_TH
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results1_TH.dat
cd ..
../rayleigh_channel 0 1  > OUTPUT_rayleigh_channel2
echo "done (start from impulsive start)"
echo " " >> validation.log
echo "Rayleigh_channel validation with impulsive start" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT_CR
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results2_CR.dat
cd ../RESLT_TH
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results2_TH.dat
cd ..

cat ./results1_CR.dat ./results1_TH.dat \
    ./results2_CR.dat ./results2_TH.dat \
    > results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         ./results.dat 0.1 1.0e-8 >> validation.log
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
