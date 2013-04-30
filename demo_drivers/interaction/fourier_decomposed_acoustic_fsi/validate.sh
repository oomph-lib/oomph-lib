#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for acoustic fsi problem
#------------------------------------
cd Validation

echo "Running Fourier-decomposed acoustic fsi validation "
mkdir RESLT
../fourier_decomposed_acoustic_fsi --nstep 2 > OUTPUT_structured
echo "done"
echo " " >> validation.log
echo "Fourier-decomposed acoustic fsi validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/helmholtz_soln1.dat RESLT/elast_soln1.dat RESLT/trace.dat > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz  \
         result.dat >> validation.log
fi

mv RESLT RESLT_structured



# Validation for unstructured acoustic fsi problem
#-------------------------------------------------

echo "Running unstructured Fourier-decomposed acoustic fsi validation "
mkdir RESLT
../unstructured_fourier_decomposed_acoustic_fsi > OUTPUT_unstructured
echo "done"
echo " " >> validation.log
echo "Unstructured Fourier-decomposed acoustic fsi validation" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/trace.dat > unstructured_result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_result.dat.gz  \
         unstructured_result.dat >> validation.log
fi

mv RESLT RESLT_unstructured


# Validation for adaptive unstructured acoustic fsi problem
#-----------------------------------------------------------

echo "Running adaptive unstructured Fourier-decomposed acoustic fsi validation "
mkdir RESLT
../adaptive_unstructured_fourier_decomposed_acoustic_fsi > OUTPUT_adaptive_unstructured
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive Fourier-decomposed acoustic fsi validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/trace.dat > adaptive_unstructured_result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_result.dat.gz  \
         adaptive_unstructured_result.dat >> validation.log
fi

mv RESLT RESLT_adaptive_unstructured


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
