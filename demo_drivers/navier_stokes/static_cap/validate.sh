#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for 2D two layer static cap
#---------------------------------------
cd Validation

echo "Running 2D static two layer validation "
mkdir RESLT
../static_two_layer > OUTPUT_static_two_layer
echo "done"
echo " " >> validation.log
echo "Two-fluid static cap validation" >> validation.log
echo "-------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln5.dat RESLT/trace.dat > static_two_layer.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/static_two_layer.dat.gz  \
   static_two_layer.dat 0.1 1.0e-8 >> validation.log
fi

mv RESLT RESLT_static_two_layer




# Validation for 2D single layer static cap
#---------------------------------------
echo "Running 2D static single layer validation "
mkdir RESLT_hijacked_external
mkdir RESLT_hijacked_internal
mkdir RESLT_elastic_hijacked_external
mkdir RESLT_elastic_hijacked_internal
../static_single_layer > OUTPUT_static_single_layer
echo "done"
echo " " >> validation.log
echo "Static single layer validation" >> validation.log
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
     > static_single_layer.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/static_single_layer.dat.gz  \
   static_single_layer.dat 0.1 1.0e-8 >> validation.log
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
