#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for sphere scattering
#---------------------------------
cd Validation

echo "Running sphere scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../sphere_scattering > OUTPUT
echo "done"
echo " " >> validation.log
echo "Sphere validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  RESLT/soln1.dat RESLT/soln2.dat  RESLT/soln3.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
fi


mv RESLT RESLT_sphere_scattering
mv OUTPUT OUTPUT_sphere_scattering

# Validation for unstructured sphere scattering
#----------------------------------------------

echo "Running unstructured sphere scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../unstructured_sphere_scattering > OUTPUT
echo "done"
echo " " >> validation.log
echo "Unstructured sphere validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_results.dat.gz   \
    unstructured_results.dat  >> validation.log
fi


mv RESLT RESLT_unstructured_sphere_scattering
mv OUTPUT OUTPUT_unstructured_sphere_scattering




# Validation for adaptive unstructured sphere scattering
#-------------------------------------------------------

echo "Running adaptive unstructured sphere scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../unstructured_adaptive_sphere_scattering > OUTPUT
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured sphere validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_unstructured_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_results.dat.gz   \
    adaptive_unstructured_results.dat  >> validation.log
fi


mv RESLT RESLT_adaptive_unstructured_sphere_scattering
mv OUTPUT OUTPUT_unstructured_sphere_scattering




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
