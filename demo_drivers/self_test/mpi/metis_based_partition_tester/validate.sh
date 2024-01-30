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


# Validation for metis based partitioning
#-------------------------------------------------------

echo "Running check metis based partitioning "
mkdir RESLT
mkdir RESLT_load_balance


# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../adaptive_driven_cavity >OUTPUT
echo "done"
echo " " >>validation.log
echo "Check metis based partitioning" >>validation.log
echo "------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT/passed_check.dat > distribution_results.dat


if test "$3" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/distribution_results.dat.gz \
                                    distribution_results.dat >>validation.log
fi


echo " " >>validation.log
echo "Check load balancing with metis based partitioning" >>validation.log
echo "--------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT_load_balance/passed_check.dat > load_balance_results.dat


if test "$3" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/load_balance_results.dat.gz \
                                    load_balance_results.dat >>validation.log
fi

#------------------------------------------------------

# Append log to main validation log
cat validation.log >>$OOMPH_ROOT_DIR/validation.log

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
