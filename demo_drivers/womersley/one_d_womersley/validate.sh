#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Womersley problems
#----------------------------------
cd Validation

echo "Running 1D impedance tube validation"
mkdir RESLT_impedance_tube
mkdir RESLT_impedance_tube_with_flux_control
../one_d_womersley -validation_run -outflow 1 > OUTPUT
../one_d_womersley -validation_run -outflow 2 >> OUTPUT
echo "done"
echo " " >> validation.log
echo "1D Womersley validation" >> validation.log
echo "-----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log



echo "1D impedance tube test: " >> validation.log
cat  RESLT_impedance_tube/soln0.dat RESLT_impedance_tube/soln1.dat RESLT_impedance_tube/soln9.dat RESLT_impedance_tube/soln10.dat RESLT_impedance_tube/trace.dat \
> 1_d_womersley.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/1_d_womersley.dat.gz  \
1_d_womersley.dat >> validation.log
fi


echo "1D Womersley with flux control test: " >> validation.log
cat  RESLT_impedance_tube_with_flux_control/soln0.dat RESLT_impedance_tube_with_flux_control/soln1.dat RESLT_impedance_tube_with_flux_control/soln9.dat RESLT_impedance_tube_with_flux_control/soln10.dat RESLT_impedance_tube_with_flux_control/trace.dat \
> 1_d_womersley_with_flux_control.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/1_d_womersley.dat.gz  \
1_d_womersley_with_flux_control.dat >> validation.log
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
