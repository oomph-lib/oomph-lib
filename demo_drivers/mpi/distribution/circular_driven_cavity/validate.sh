#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


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

# Validation for circular driven cavity (TH and CR problems)
#-----------------------------------------------------------
mkdir RESLT
mkdir RESLT_MESH
mkdir RESLT_ERROR

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

echo "Running adaptive circular driven cavity (TH, CR problems) validation "
$MPI_RUN_COMMAND ../circular_driven_cavity > OUTPUT_circular_driven_cavity
echo "done"
echo " " >> validation.log
echo "Adaptive circular driven cavity validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln1_on_proc0.dat RESLT/soln1_on_proc1.dat > circular_cavity_TH1_results.dat
cat RESLT/soln2_on_proc0.dat RESLT/soln2_on_proc1.dat > circular_cavity_TH2_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/circular_cavity_TH1_results.dat.gz  \
         circular_cavity_TH1_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/circular_cavity_TH2_results.dat.gz  \
         circular_cavity_TH2_results.dat >> validation.log
fi

mv RESLT RESLT_circular_driven_cavity

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
