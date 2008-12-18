#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation



# Validation for fish poisson without adaptation
#-----------------------------------------------

echo "Running fish poisson without adaptation validation "
mkdir RESLT
../fish_poisson_no_adapt > OUTPUT_fish_poisson_no_adapt
echo "done"
echo " " >> validation.log
echo "Fish poisson without adaptation validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > fish_poisson_no_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_no_adapt.dat.gz   \
    fish_poisson_no_adapt.dat  >> validation.log
fi
mv RESLT RESLT_fish_poisson_no_adapt




# Validation for fish poisson with adaptation
#--------------------------------------------

echo "Running fish poisson with adaptation validation "
mkdir RESLT
../fish_poisson_adapt > OUTPUT_fish_poisson_adapt
echo "done"
echo " " >> validation.log
echo "Fish poisson with adaptation validation" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln3.dat > fish_poisson_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_adapt.dat.gz   \
    fish_poisson_adapt.dat  >> validation.log
fi
mv RESLT RESLT_fish_poisson_adapt



# Validation for fish poisson with node updates
#------------------------------------------------

echo "Running fish poisson with node updates validation "
mkdir RESLT
../fish_poisson_node_update > OUTPUT_fish_poisson_node_update
echo "done"
echo " " >> validation.log
echo "Fish poisson with node update validation" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln2.dat > fish_poisson_node_update.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_node_update.dat.gz   \
    fish_poisson_node_update.dat  >> validation.log
fi
mv RESLT RESLT_fish_poisson_node_update





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




