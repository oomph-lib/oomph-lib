#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for buckling of clamped shell 
#-----------------------------------------

mkdir RESLT_clamped

echo "Running clamped shell "
../clamped_or_pinned_shell 0 > OUTPUT_clamped

echo "done"
echo " " >> validation.log
echo "Clamped shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_clamped/shell0.dat RESLT_clamped/shell4.dat RESLT_clamped/shell9.dat RESLT_clamped/trace.dat > clamped_shell.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/clamped_shell.dat.gz \
 clamped_shell.dat>> validation.log
fi




# Validation for buckling of pinned shell 
#-----------------------------------------

mkdir RESLT_pinned

echo "Running pinned shell "
../clamped_or_pinned_shell 1 > OUTPUT_clamped

echo "done"
echo " " >> validation.log
echo "Pinned shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_pinned/shell0.dat RESLT_pinned/shell4.dat RESLT_pinned/shell9.dat RESLT_pinned/trace.dat > pinned_shell.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pinned_shell.dat.gz \
 pinned_shell.dat>> validation.log
fi



# Validation for buckling of pinned periodic shell 
#-----------------------------------------

mkdir RESLT_pinned_periodic

echo "Running pinned periodic shell "
../clamped_or_pinned_shell 2 > OUTPUT_clamped

echo "done"
echo " " >> validation.log
echo "Pinned periodic shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_pinned_periodic/shell0.dat RESLT_pinned_periodic/shell4.dat RESLT_pinned_periodic/shell9.dat RESLT_pinned_periodic/trace.dat > pinned_periodic.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pinned_periodic.dat.gz \
 pinned_periodic.dat>> validation.log
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
