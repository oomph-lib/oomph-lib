#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=7


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation
cd Validation


# Validation for demo Stefan Boltzmann problem
#---------------------------------------------

echo "Running 2D Stefan Boltzmann validation"
mkdir RESLT
../stefan_boltzmann > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D Stefan Boltzmann validation" >> validation.log
echo "------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_sb.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_sb.dat.gz \
    result_sb.dat  >> validation.log
fi

mv RESLT RESLT_sb


# Validation for demo unsteady heat with pretend melting
#-------------------------------------------------------

echo "Running 2D pretend melt validation "
mkdir RESLT
../pretend_melt --validate > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady pretend melt validation " >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_pretend_melt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_pretend_melt.dat.gz \
    result_pretend_melt.dat  >> validation.log
fi

mv RESLT RESLT_pretend_melt


# Validation for demo unsteady heat with pretend melting (2)
#-----------------------------------------------------------

echo "Running 2D pretend melt validation (2)"
mkdir RESLT
../melt --disable_melting --validate > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady pretend melt validation (2)" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_pretend_melt2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_pretend_melt2.dat.gz \
    result_pretend_melt2.dat  >> validation.log
fi

mv RESLT RESLT_pretend_melt2




# Validation for demo unsteady heat with actual melting
#------------------------------------------------------

echo "Running 2D actual melt validation"
mkdir RESLT
../melt --validate > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady actual melt validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_melt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_melt.dat.gz \
    result_melt.dat  >> validation.log
fi

mv RESLT RESLT_melt



# Validation for spring contact
#------------------------------

echo "Running spring contact (single kink) validation "
mkdir RESLT
cd RESLT
../../spring_contact --single_kink > OUTPUT
cd ..
echo "done"
echo " " >> validation.log
echo "Spring contact (single kink) validation " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_spring_contact_single_kink.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_spring_contact_single_kink.dat.gz \
    result_spring_contact_single_kink.dat  >> validation.log
fi

mv RESLT RESLT_spring_contact_single_kink


# Validation for spring contact
#------------------------------

echo "Running spring contact (kuhn tucker) validation "
mkdir RESLT
cd RESLT
../../spring_contact --kuhn_tucker > OUTPUT
cd ..
echo "done"
echo " " >> validation.log
echo "Spring contact (kuhn tucker) validation " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_spring_contact_kuhn_tucker.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_spring_contact_kuhn_tucker.dat.gz \
    result_spring_contact_kuhn_tucker.dat  >> validation.log
fi

mv RESLT RESLT_spring_contact_kuhn_tucker



# Validation for spring contact
#------------------------------

echo "Running spring contact (old version) validation "
mkdir RESLT
cd RESLT
../../spring_contact --old_version > OUTPUT
cd ..
echo "done"
echo " " >> validation.log
echo "Spring contact (old version) validation " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_spring_contact_old_version.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_spring_contact_old_version.dat.gz \
    result_spring_contact_old_version.dat  >> validation.log
fi

mv RESLT RESLT_spring_contact_old_version




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
