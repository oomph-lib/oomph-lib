#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for static disk compression
#---------------------------------------

cd Validation
mkdir RESLT_fold RESLT_pitch RESLT_hopf

echo "Running fold bifurcation validation "
cd RESLT_fold
../../fold > ../OUTPUT_fold
cd ..

echo "done"
echo " " >> validation.log
echo "Fold bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_fold/trace_mu.dat > fold.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/fold.dat.gz \
    fold.dat  >> validation.log
fi


echo "Running pitchfork bifurcation validation "
cd RESLT_pitch
../../pitchfork > ../OUTPUT_pitchfork
cd ..

echo "done"
echo " " >> validation.log
echo "Pitchfork bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_pitch/trace.dat > pitch.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/pitch.dat.gz pitch.dat  >> validation.log
fi

echo "Running hopf bifurcation validation "
cd RESLT_hopf
../../hopf > ../OUTPUT_hopf
cd ..

echo "done"
echo " " >> validation.log
echo "Hopf bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_hopf/trace_hopf.dat > hopf.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/hopf.dat.gz \
    hopf.dat  0.3 1.0e-14  >> validation.log
fi



# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../validation.log





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

