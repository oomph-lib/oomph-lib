#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3

# Doc what we're using to run tests on two processors
echo " " 
echo "Running mpi tests with mpi run command: " $MPI_RUN_COMMAND
echo " " 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation
cp ../*partition.dat .


#----------------------------------------------------------------------

# Validation for fish poisson
#----------------------------

echo "Running fish poisson validation "
mkdir RESLT_select_refine
mkdir RESLT_incremental2
mkdir RESLT_fully_automatic

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../fish_poisson > OUTPUT_fish_poisson
echo "done"
echo " " >> validation.log
echo "Fish poisson validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_select_refine/soln0_on_proc0.dat RESLT_select_refine/soln0_on_proc1.dat \
    RESLT_select_refine/soln1_on_proc0.dat RESLT_select_refine/soln1_on_proc1.dat \
    RESLT_select_refine/soln2_on_proc0.dat RESLT_select_refine/soln2_on_proc1.dat \
    > fish_poisson_select_results.dat
cat RESLT_incremental2/soln1_on_proc0.dat RESLT_incremental2/soln1_on_proc1.dat \
    RESLT_incremental2/soln6_on_proc0.dat RESLT_incremental2/soln6_on_proc1.dat \
    RESLT_incremental2/soln16_on_proc0.dat RESLT_incremental2/soln16_on_proc1.dat \
    > fish_poisson_incremental_results.dat
cat RESLT_fully_automatic/soln0_on_proc0.dat RESLT_fully_automatic/soln0_on_proc1.dat \
    > fish_poisson_automatic_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/fish_poisson_select_results.dat.gz  \
         fish_poisson_select_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/fish_poisson_incremental_results.dat.gz  \
         fish_poisson_incremental_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/fish_poisson_automatic_results.dat.gz  \
         fish_poisson_automatic_results.dat >> validation.log
fi

mkdir RESLT_fish_poisson
mv RESLT_select_refine RESLT_fish_poisson
mv RESLT_incremental2 RESLT_fish_poisson
mv RESLT_fully_automatic RESLT_fish_poisson

#----------------------------------------------------------------------

# Append log to main validation log
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
