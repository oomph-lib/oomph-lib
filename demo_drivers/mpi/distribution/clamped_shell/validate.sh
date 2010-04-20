#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for buckling of clamped shell with arclength continuation
#---------------------------------------------------------------------

mkdir RESLT
cd RESLT

echo "Running clamped_shell with arclength continuation validation  "
$MPI_RUN_COMMAND ../../clamped_shell_with_arclength_cont > ../OUTPUT_clamped_shell_with_arclength_cont

echo "done"
cd ..
echo " " >> validation.log
echo "Clamped shell buckling with arclength continuation validation" \
 >> validation.log
echo "(parallel continuation test)" >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/final_shape_on_proc0.dat RESLT/final_shape_on_proc1.dat RESLT/trace_on_proc0.dat RESLT/trace_disp_on_proc0.dat > shell_with_arclength_cont_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/shell_with_arclength_cont_results.dat.gz \
 shell_with_arclength_cont_results.dat 0.1 1.0e-9 >> validation.log
fi


mv RESLT RESLT_with_arclength




# Append output to global validation log file
#--------------------------------------------
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
