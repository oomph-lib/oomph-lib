#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for coarse vs fine base meshes demo
#------------------------------------------------
cd Validation

echo "Running coarse vs fine base meshes demo "
mkdir RESLT
echo "S" >> tmp_input.dat
../two_d_poisson_adapt 3 < tmp_input.dat > OUTPUT_two_d_poisson_adapt
rm tmp_input.dat
echo "done"
echo " " >> validation.log
echo "Coarse vs fine base meshes demo" >> validation.log
echo "-------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.cpp_style.dat  RESLT/soln1.cpp_style.dat  \
    > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Compare results: " >> validation.log
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
echo "Compare C and C++ output: " >> validation.log
../../../../bin/fpdiff.py RESLT/soln0.cpp_style.dat RESLT/soln0.c_style.dat\
    >> validation.log
../../../../bin/fpdiff.py RESLT/soln1.cpp_style.dat RESLT/soln1.c_style.dat\
    >> validation.log
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




