#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Navier Stokes flow in a curved pipe problem
#---------------------------------------------------------
cd Validation

echo "Running curved pipe Navier Stokes validation "
mkdir RESLT_TH
mkdir RESLT_CR

../curved_pipe lalala > OUTPUT_curved_pipe
echo "done"
echo " " >> validation.log
echo "3D curved pipe Navier Stokes validation" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln_Re50.dat  RESLT_CR/soln_Re50.dat > curved_pipe.dat
mv RESLT_TH RESLT_TH_cyl
mv RESLT_CR RESLT_CR_cyl

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/curved_pipe.dat.gz  \
         curved_pipe.dat 0.1 1.0e-13 >> validation.log
fi

echo "Running helical pipe Navier Stokes validation "
mkdir RESLT_TH
mkdir RESLT_CR

../helical_pipe lalala > OUTPUT_helical_pipe
echo "done"
echo " " >> validation.log
echo "3D helical pipe Navier Stokes validation" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln_Re50.dat  RESLT_CR/soln_Re50.dat > helical_pipe.dat
mv RESLT_TH RESLT_TH_helix
mv RESLT_CR RESLT_CR_helix

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/helical_pipe.dat.gz  \
         helical_pipe.dat 0.1 1.0e-13 >> validation.log
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
