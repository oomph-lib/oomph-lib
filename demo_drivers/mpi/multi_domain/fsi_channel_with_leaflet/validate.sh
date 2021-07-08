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

# Validation for FSI channel with leaflet problem (algebraic node update)
#------------------------------------------------------------------------

echo "Running FSI channel with leaflet (algebraic node update) validation "
mkdir RESLT_FSI_LEAF

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../fsi_channel_with_leaflet validate > OUTPUT_fsi_channel_with_leaflet
echo "done"
echo " " >> validation.log
echo "FSI channel with leaflet (algebraic) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_FSI_LEAF/soln1_on_proc0.dat RESLT_FSI_LEAF/soln1_on_proc1.dat RESLT_FSI_LEAF/soln3_on_proc0.dat RESLT_FSI_LEAF/soln3_on_proc1.dat RESLT_FSI_LEAF/wall_soln2_on_proc0.dat RESLT_FSI_LEAF/wall_soln2_on_proc1.dat > fsi_channel_with_leaflet_external_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/fsi_channel_with_leaflet_external_results.dat.gz  \
         fsi_channel_with_leaflet_external_results.dat 0.4 1.0e-14 >> validation.log
fi

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
