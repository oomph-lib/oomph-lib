#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Rayleigh channel
#--------------------------------
cd Validation

echo "Running rayleigh channel with/without ALE"
mkdir RESLT_CR RESLT_TH RESLT_CR_ALE RESLT_TH_ALE
../rayleigh_channel 0 0 > OUTPUT_rayleigh_channel
echo "done"
echo " " >> validation.log
echo "Rayleigh_channel with/without ALE" >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT_CR
cat  soln0.dat soln1.dat > ../results_CR.dat
cd ../RESLT_TH
cat  soln0.dat soln1.dat > ../results_TH.dat
cd ../RESLT_CR_ALE
cat  soln0.dat soln1.dat > ../results_CR_ALE.dat
cd ../RESLT_TH_ALE
cat  soln0.dat soln1.dat > ../results_TH_ALE.dat
cd ..

cat ./results_CR.dat ./results_TH.dat \
    ./results_CR_ALE.dat ./results_TH_ALE.dat > results.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/results.dat.gz  results.dat 0.1 1.0e-8 >> validation.log
fi


# Validation for Rayleigh driven cavity
#--------------------------------------

echo "Running rayleigh driven cavity with/without ALE"
mkdir RESLT2_CR RESLT2_TH RESLT2_CR_ALE RESLT2_TH_ALE
../rayleigh_circular_driven_cavity_adapt 0 0 > OUTPUT_rayleigh_driven_cavity
echo "done"
echo " " >> validation.log
echo "Rayleigh driven cavity with/without ALE" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT2_CR
cat  soln0.dat soln1.dat > ../results2_CR.dat
cd ../RESLT2_TH
cat  soln0.dat soln1.dat > ../results2_TH.dat
cd ../RESLT2_CR_ALE
cat  soln0.dat soln1.dat > ../results2_CR_ALE.dat
cd ../RESLT2_TH_ALE
cat  soln0.dat soln1.dat > ../results2_TH_ALE.dat
cd ..

cat ./results2_CR.dat ./results2_TH.dat \
    ./results2_CR_ALE.dat ./results2_TH_ALE.dat > results2.dat

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/results2.dat.gz  results2.dat 0.1 1.0e-8 >> validation.log
fi



# Append log to main validation log
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
