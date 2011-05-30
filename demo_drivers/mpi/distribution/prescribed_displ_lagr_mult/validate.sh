#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=4

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


#------------------------------------------------------

# Validation for prescribed displacement problem with Lagrange multipliers
#-------------------------------------------------------------------------
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

echo "Running prescribed displ. problem with Lagrange multipliers "
$MPI_RUN_COMMAND ../prescribed_displ_lagr_mult --validation > OUTPUT_prescribed_displ_lagr_mult
echo "done"
echo " " >> validation.log
echo "Prescribed displ. with Lagrange multipliers validation" >> validation.log
echo "------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0_on_proc0.dat RESLT/soln1_on_proc1.dat RESLT/soln2_on_proc0.dat \
    RESLT/lagr0_on_proc1.dat RESLT/lagr1_on_proc0.dat RESLT/lagr2_on_proc1.dat \
    > prescribed_displ_lagr_mult_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/prescribed_displ_lagr_mult_results.dat.gz  \
         prescribed_displ_lagr_mult_results.dat >> validation.log
fi


cat RESLT/soln3_on_proc0.dat RESLT/soln3_on_proc1.dat \
    RESLT/soln4_on_proc0.dat RESLT/soln4_on_proc1.dat \
    RESLT/lagr3_on_proc0.dat RESLT/lagr3_on_proc1.dat  \
    RESLT/lagr4_on_proc0.dat RESLT/lagr4_on_proc1.dat  \
    > prescribed_displ_lagr_mult_results_load_balanced.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/prescribed_displ_lagr_mult_results_load_balanced.dat.gz  \
         prescribed_displ_lagr_mult_results_load_balanced.dat >> validation.log
fi

mv RESLT RESLT_prescribed_displ_lagr_mult

#------------------------------------------------------


# Validation for testing resizing halo nodes
#-------------------------------------------
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

echo "Running test for resizing halo nodes "
$MPI_RUN_COMMAND ../resize_hanging_node_tester --validation > OUTPUT_resize_hanging_node_tester
echo "done"
echo " " >> validation.log
echo "Resizing halo nodes test" >> validation.log
echo "------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0_on_proc0.dat RESLT/soln1_on_proc1.dat RESLT/soln2_on_proc0.dat \
    RESLT/lagr0_on_proc1.dat RESLT/lagr1_on_proc0.dat RESLT/lagr2_on_proc1.dat \
    RESLT/nodes0_on_proc0.dat RESLT/nodes0_on_proc1.dat \
    RESLT/nodes1_on_proc0.dat RESLT/nodes1_on_proc1.dat \
    RESLT/nodes2_on_proc0.dat RESLT/nodes2_on_proc1.dat \
    > resize_test_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/resize_test_results.dat.gz  \
         resize_test_results.dat >> validation.log
fi


cat RESLT/soln3_on_proc0.dat RESLT/soln3_on_proc1.dat \
    RESLT/soln4_on_proc0.dat RESLT/soln4_on_proc1.dat \
    RESLT/lagr3_on_proc0.dat RESLT/lagr3_on_proc1.dat  \
    RESLT/lagr4_on_proc0.dat RESLT/lagr4_on_proc1.dat  \
    RESLT/nodes3_on_proc0.dat RESLT/nodes3_on_proc1.dat  \
    RESLT/nodes4_on_proc0.dat RESLT/nodes4_on_proc1.dat  \
    > resize_test_load_balanced.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/resize_test_load_balanced.dat.gz  \
         resize_test_load_balanced.dat >> validation.log
fi

mv RESLT RESLT_resize_test

#------------------------------------------------------

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
