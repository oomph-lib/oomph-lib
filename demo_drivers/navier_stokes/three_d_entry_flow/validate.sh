#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for 3D entry flow Navier Stokes problem (quarter tube)
#-------------------------------------------------------------------
cd Validation

echo "Running 3D entry flow (quarter tube) Navier Stokes validation "
mkdir RESLT_TH
mkdir RESLT_CR

../three_d_entry_flow lalala > OUTPUT
echo "done"
echo " " >> validation.log
echo "3D entry flow (quarter tube) Navier Stokes validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln1.dat  RESLT_CR/soln1.dat > entry_flow_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/entry_flow_results.dat.gz  \
         entry_flow_results.dat >> validation.log
fi

mv RESLT_TH RESLT_TH_quarter
mv RESLT_CR RESLT_CR_quarter

# Validation for 3D entry flow (full tube) Navier Stokes problem
#--------------------------------------------------------------------

echo "Running 3D entry flow (full tube) Navier Stokes validation "
mkdir RESLT_TH

../full_tube > OUTPUT_full
echo "done"
echo " " >> validation.log
echo "3D entry flow (full tube) Navier Stokes validation" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln_Re50.dat  > full_tube.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/full_tube.dat.gz  \
         full_tube.dat >> validation.log
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
