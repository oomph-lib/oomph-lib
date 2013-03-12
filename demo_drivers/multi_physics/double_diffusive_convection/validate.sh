#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running 2D double-diffusive convection validation "
mkdir RESLT
../dd_convection lalala > ./OUTPUT_dd
echo "done"
echo " " >> validation.log
echo "2D double-diffusive convection validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/soln5.dat > dd.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/dd.dat.gz \
    dd.dat  0.1 1.0e-14 >> validation.log
fi
mv RESLT RESLT_dd


echo "Running 2D double-diffusive convection (multimesh) validation "
mkdir RESLT
../multimesh_dd_convection lalala > ./OUTPUT_dd_multimesh
echo "done"
echo " " >> validation.log
echo "2D double-diffusive convection (multimesh) validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/soln5.dat > dd_multimesh.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/dd_multimesh.dat.gz \
    dd_multimesh.dat  0.1 1.0e-12 >> validation.log
fi
mv RESLT RESLT_dd_multimesh

echo "Running Refineable 2D double-diffusive convection (multimesh) validation "
mkdir RESLT_ref_multimesh
../multimesh_ref_dd_convection lalala > ./OUTPUT_ref_dd_multimesh
echo "done"
echo " " >> validation.log
echo "Refineable 2D double-diffusive convection (multimesh) validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_ref_multimesh/trace.dat RESLT_ref_multimesh/nst_soln5.dat \
    RESLT_ref_multimesh/temp_soln5.dat RESLT_ref_multimesh/conc_soln5.dat \
     > ref_dd_multimesh.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/ref_dd_multimesh.dat.gz \
    ref_dd_multimesh.dat 0.1 2.0e-12 >> validation.log
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
