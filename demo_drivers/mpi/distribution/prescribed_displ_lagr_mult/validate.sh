#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


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
mkdir RESLT_resize_test
mkdir RESLT

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

cat RESLT_resize_test/soln0_on_proc0.dat RESLT_resize_test/soln1_on_proc1.dat RESLT_resize_test/soln2_on_proc0.dat \
    RESLT_resize_test/lagr0_on_proc1.dat RESLT_resize_test/lagr1_on_proc0.dat RESLT_resize_test/lagr2_on_proc1.dat \
    RESLT_resize_test/nodes0_on_proc0.dat RESLT_resize_test/nodes0_on_proc1.dat \
    RESLT_resize_test/nodes1_on_proc0.dat RESLT_resize_test/nodes1_on_proc1.dat \
    RESLT_resize_test/nodes2_on_proc0.dat RESLT_resize_test/nodes2_on_proc1.dat \
    > resize_test_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/resize_test_results.dat.gz  \
         resize_test_results.dat >> validation.log
fi


cat RESLT_resize_test/soln3_on_proc0.dat RESLT_resize_test/soln3_on_proc1.dat \
    RESLT_resize_test/soln4_on_proc0.dat RESLT_resize_test/soln4_on_proc1.dat \
    RESLT_resize_test/lagr3_on_proc0.dat RESLT_resize_test/lagr3_on_proc1.dat  \
    RESLT_resize_test/lagr4_on_proc0.dat RESLT_resize_test/lagr4_on_proc1.dat  \
    RESLT_resize_test/nodes3_on_proc0.dat RESLT_resize_test/nodes3_on_proc1.dat  \
    RESLT_resize_test/nodes4_on_proc0.dat RESLT_resize_test/nodes4_on_proc1.dat  \
    > resize_test_load_balanced.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/resize_test_load_balanced.dat.gz  \
         resize_test_load_balanced.dat >> validation.log
fi

#------------------------------------------------------

# Append log to main validation log
cat validation.log >> ../../../../../validation.log

cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
