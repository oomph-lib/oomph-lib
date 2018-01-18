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
cp ../*partition.dat .

# Validation for FSI oscillating ring problem (algebraic node update)
#--------------------------------------------------------------------

echo "Running FSI oscillating ring (algebraic node update) validation "
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../fsi_osc_ring validate > OUTPUT_fsi_osc_ring
echo "done"
echo " " >> validation.log
echo "FSI oscillating ring (algebraic) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/new_soln0_on_proc0.dat RESLT/new_soln0_on_proc1.dat RESLT/new_soln1_on_proc0.dat RESLT/new_soln1_on_proc1.dat RESLT/new_wall_soln1_on_proc0.dat RESLT/new_wall_soln1_on_proc1.dat > fsi_osc_ring_external_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
# Note these limits are ridiculous but the solution is VERY rough
../../../../../bin/fpdiff.py ../validata/fsi_osc_ring_external_results.dat.gz  \
         fsi_osc_ring_external_results.dat 10.0 1.0e-3 >> validation.log
fi

mv RESLT RESLT_fsi_osc_ring

#---------------------------------------------------------------------

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
