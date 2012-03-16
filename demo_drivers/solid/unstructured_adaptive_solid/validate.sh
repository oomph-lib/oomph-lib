#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for static fish deformation
#---------------------------------------

cd Validation
mkdir RESLT RESLT_pres_disp RESLT_pres_disp_incomp

echo "Running unstructured adaptive solid validation "
../unstructured_adaptive_solid > OUTPUT_adaptive_solid


echo "done"
echo " " >> validation.log
echo "Unstructured adaptive solid validation" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/s_energy.dat > res.dat
cat RESLT_pres_disp/s_energy.dat > res_pres_disp.dat
cat RESLT_pres_disp_incomp/s_energy.dat > res_pres_disp_incomp.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/res.dat.gz \
    res.dat 0.1 1.0e-12 >> validation.log
../../../../bin/fpdiff.py ../validata/res_pres_disp.dat.gz \
    res_pres_disp.dat 0.1 1.0e-7 >> validation.log
../../../../bin/fpdiff.py ../validata/res_pres_disp_incomp.dat.gz \
    res_pres_disp_incomp.dat 0.11 1.5e-7 >> validation.log
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

