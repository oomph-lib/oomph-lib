#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

# Receive the mpirun command as the second argument
MPI_RUN_COMMAND="$2"

#Set the number of tests to be checked
NUM_TESTS=2

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

# Validation for 3D entry flow problem (TH & CR)
#-----------------------------------------------

echo "Running 3D entry flow (TH & CR) validation "
mkdir RESLT_TH
mkdir RESLT_CR
mkdir RESLT_MESH

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../three_d_entry_flow validate > OUTPUT_three_d_entry_flow
echo "done"
echo " " >> validation.log
echo "3D entry flow (TH & CR) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_TH/soln0_on_proc0.dat RESLT_TH/soln0_on_proc1.dat \
    RESLT_TH/soln3_on_proc0.dat RESLT_TH/soln3_on_proc1.dat \
    > three_d_entry_flow_TH_results.dat
cat RESLT_CR/soln1_on_proc0.dat RESLT_CR/soln1_on_proc1.dat \
    RESLT_CR/soln3_on_proc0.dat RESLT_CR/soln3_on_proc1.dat \
    > three_d_entry_flow_CR_results.dat

if test "$3" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/three_d_entry_flow_TH_results.dat.gz  \
         three_d_entry_flow_TH_results.dat >> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/three_d_entry_flow_CR_results.dat.gz  \
         three_d_entry_flow_CR_results.dat >> validation.log
fi

mkdir RESLT_three_d_entry_flow
mv RESLT_TH RESLT_three_d_entry_flow
mv RESLT_CR RESLT_three_d_entry_flow

#----------------------------------------------------------------------

# Append log to main validation log
cat validation.log >> $OOMPH_ROOT_DIR/validation.log

cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 10
