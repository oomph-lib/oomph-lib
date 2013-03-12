#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running 1D demo spectral poisson validation "
mkdir RESLT_1d
cd RESLT_1d
../../one_d_spectral > ../OUTPUT_one_d_spectral
cd ..
echo "done"
echo " " >> validation.log
echo "1D demo spectral poisson validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_1d/soln0.dat RESLT_1d/soln1.dat > 1d_spectral.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/1d_spectral.dat.gz   \
    1d_spectral.dat  >> validation.log
fi

echo "Running 2D demo spectral poisson validation "
mkdir RESLT_2d
../two_d_spectral_adapt > ./OUTPUT_two_d_spectral_adapt
echo "done"
echo " " >> validation.log
echo "2D demo spectral poisson validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_2d/soln0.dat > 2d_spectral.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_spectral.dat.gz   \
    2d_spectral.dat  >> validation.log
fi

echo "Running eighth sphere spectral poisson validation "
mkdir RESLT_3d
../three_d_spectral_adapt lala > OUTPUT_three_d_spectral_adapt
echo "done"
echo " " >> validation.log
echo "Eighth sphere spectral poisson validation" >> validation.log
echo "--------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_3d/soln0.dat RESLT_3d/soln1.dat > 3d_spectral.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/3d_spectral.dat.gz   \
    3d_spectral.dat  0.1 1.0e-8 >> validation.log
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
