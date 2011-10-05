#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for axisymmetric static cap
#---------------------------------------
cd Validation

echo "Running axisymmetric static cap validation "
mkdir RESLT_hijacked_external
mkdir RESLT_hijacked_internal
mkdir RESLT_elastic_hijacked_external
mkdir RESLT_elastic_hijacked_internal
#valgrind --leak-check=full -v 
../axi_static_cap > OUTPUT_axi_static_cap
echo "done"
echo " " >> validation.log
echo "Static axisymmetric cap validation" >> validation.log
echo "------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_hijacked_external/soln5.dat RESLT_hijacked_external/trace.dat \
     RESLT_hijacked_internal/soln5.dat RESLT_hijacked_internal/trace.dat \
     RESLT_elastic_hijacked_external/soln5.dat \
     RESLT_elastic_hijacked_external/trace.dat \
     RESLT_elastic_hijacked_internal/soln5.dat \
     RESLT_elastic_hijacked_internal/trace.dat \
     > axi_static_cap.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/axi_static_cap.dat.gz  \
   axi_static_cap.dat 0.1 1.0e-8 >> validation.log
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
