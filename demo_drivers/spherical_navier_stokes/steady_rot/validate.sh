#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for spherical Couette flow
#--------------------------------
cd Validation

echo "Running Steady Spherical Rigid-Body Rotation validation "
mkdir RESLT
cd RESLT
../../steady_rot > ../OUTPUT_steady_rot
cd ..
echo "done"
echo " " >> validation.log
echo "Steady spherical rigid-body rotation validation" >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln_CR_8x8_10.dat RESLT/trace_CR.dat > sph_rot_CR.dat
cat  RESLT/soln_TH_8x8_10.dat RESLT/trace_TH.dat > sph_rot_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Crouzeix-Raviart elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/sph_rot_CR.dat.gz  \
         sph_rot_CR.dat 0.1 1.0e-12 >> validation.log
echo "Taylor-Hood elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/sph_rot_TH.dat.gz  \
         sph_rot_TH.dat 0.1 1.0e-12 >> validation.log

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
