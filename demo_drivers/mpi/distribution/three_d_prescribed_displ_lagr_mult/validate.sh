#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1

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


#----------------------------------------------------------------------------
# Validation for 3d prescribed displacement problem with Lagrange multipliers
#----------------------------------------------------------------------------
mkdir RESLT_proc0
mkdir RESLT_proc1

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

echo "Running 3d prescribed displ. problem with Lagrange multipliers "
$MPI_RUN_COMMAND ../three_d_prescribed_displ_lagr_mult bla > OUTPUT
echo "done"
echo " " >> validation.log
echo "3D prescribed displ. with Lagrange multipliers validation" 
>> validation.log
echo "---------------------------------------------------------" 
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_proc0/soln1.dat  RESLT_proc0/lagr1.dat  \
    RESLT_proc1/soln1.dat  RESLT_proc1/lagr1.dat  \
    > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         results.dat >> validation.log
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
