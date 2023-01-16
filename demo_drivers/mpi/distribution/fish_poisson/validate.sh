#! /bin/sh

echo "Arguments $#: $@"

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

# Receive the mpirun command as the second argument
MPI_RUN_COMMAND="$2"

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
mkdir RESLT_select_mesh
mkdir RESLT_adapt_mesh
mkdir RESLT_incr_mesh
mkdir RESLT_incremental2
mkdir RESLT_fully_automatic

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../fish_poisson >OUTPUT_fish_poisson
echo "done"
echo " " >>validation.log
echo "Fish poisson validation" >>validation.log
echo "------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT_select_refine/soln0_on_proc0.dat RESLT_select_refine/soln0_on_proc1.dat \
    RESLT_select_refine/soln1_on_proc0.dat RESLT_select_refine/soln1_on_proc1.dat \
    RESLT_select_refine/soln2_on_proc0.dat RESLT_select_refine/soln2_on_proc1.dat \
    >fish_poisson_select_results.dat
cat RESLT_incremental2/soln1_on_proc0.dat RESLT_incremental2/soln1_on_proc1.dat \
    RESLT_incremental2/soln6_on_proc0.dat RESLT_incremental2/soln6_on_proc1.dat \
    RESLT_incremental2/soln16_on_proc0.dat RESLT_incremental2/soln16_on_proc1.dat \
    >fish_poisson_incremental_results.dat
cat RESLT_fully_automatic/soln0_on_proc0.dat RESLT_fully_automatic/soln0_on_proc1.dat \
    >fish_poisson_automatic_results.dat

if test "$3" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
    $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/fish_poisson_select_results.dat.gz \
        fish_poisson_select_results.dat >>validation.log
    $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/fish_poisson_incremental_results.dat.gz \
        fish_poisson_incremental_results.dat >>validation.log
    $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/fish_poisson_automatic_results.dat.gz \
        fish_poisson_automatic_results.dat >>validation.log
fi

mkdir RESLT_fish_poisson
mv RESLT_select_refine RESLT_fish_poisson
mv RESLT_incremental2 RESLT_fish_poisson
mv RESLT_fully_automatic RESLT_fish_poisson

#----------------------------------------------------------------------

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
