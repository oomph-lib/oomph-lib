#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


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

# Validation for adaptive driven cavity (TH & CR)
#------------------------------------------------

echo "Running adaptive rectangular driven cavity (TH & CR) validation "
mkdir RESLT_CR_MESH  RESLT RESLT_TH_MESH


# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../adaptive_driven_cavity validate > OUTPUT_adaptive_driven_cavity
echo "done"
echo " " >> validation.log
echo "Adaptive rectangular driven cavity (TH & CR) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > adaptive_cavity_TH_results.dat
cat RESLT/soln1_on_proc0.dat RESLT/soln1_on_proc1.dat \
    > adaptive_cavity_CR_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/adaptive_cavity_TH_results.dat.gz  \
         adaptive_cavity_TH_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/adaptive_cavity_CR_results.dat.gz  \
         adaptive_cavity_CR_results.dat >> validation.log
fi

mv RESLT RESLT_adaptive_driven_cavity
mv RESLT_TH_MESH  RESLT_TH_MESH_adaptive_driven_cavity
mv RESLT_CR_MESH  RESLT_CR_MESH_adaptive_driven_cavity

#----------------------------------------------------------------------

# Validation for adaptive driven cavity (TH & CR) with load balancing
#--------------------------------------------------------------------

echo "Running adaptive rectangular driven cavity validation with load balancing"
mkdir RESLT_LOAD_BALANCE 

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../adaptive_driven_cavity_load_balance --validate > OUTPUT_adaptive_driven_cavity_load_balance
echo "done"
echo " " >> validation.log
echo "Adaptive rectangular driven cavity validation with load balancing" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_LOAD_BALANCE/soln0_on_proc0.dat RESLT_LOAD_BALANCE/soln0_on_proc1.dat \
    RESLT_LOAD_BALANCE/soln9_on_proc0.dat RESLT_LOAD_BALANCE/soln9_on_proc1.dat \
    > adaptive_cavity_load_balance_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/adaptive_cavity_load_balance_results.dat.gz  \
         adaptive_cavity_load_balance_results.dat >> validation.log
fi


mv RESLT_LOAD_BALANCE RESLT_adaptive_driven_cavity_load_balance

#----------------------------------------------------------------------

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
