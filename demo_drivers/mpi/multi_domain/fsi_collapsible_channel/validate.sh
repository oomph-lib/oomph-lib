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

# Validation for FSI collapsible channel with macro element node update
#----------------------------------------------------------------------

echo "Running FSI collapsible channel problem (macro element node update) "
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../fsi_collapsible_channel_macro_adapt validate > OUTPUT_fsi_collapsible_channel_macro_adapt
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel (macro element node update) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/new_soln0_on_proc0.dat RESLT/new_soln1_on_proc1.dat \
    RESLT/new_soln2_on_proc0.dat RESLT/new_soln3_on_proc1.dat \
    > fsi_collapsible_channel_macro_adapt_external_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/fsi_collapsible_channel_macro_adapt_external_results.dat.gz  \
         fsi_collapsible_channel_macro_adapt_external_results.dat 2.0 5.0e-5 >> validation.log
fi

mv RESLT RESLT_fsi_collapsible_channel_macro_adapt_external

#----------------------------------------------------------------------

# Validation for FSI collapsible channel with algebraic node update
#----------------------------------------------------------------------

echo "Running FSI collapsible channel problem (algebraic node update) "
mkdir RESLT_ALG

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../fsi_collapsible_channel_adapt validate > OUTPUT_fsi_collapsible_channel_adapt
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel (algebraic node update) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_ALG/new_soln0_on_proc0.dat RESLT_ALG/new_soln1_on_proc1.dat \
    RESLT_ALG/new_soln2_on_proc0.dat RESLT_ALG/new_soln3_on_proc1.dat \
    > fsi_collapsible_channel_adapt_external_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/fsi_collapsible_channel_adapt_external_results.dat.gz  \
         fsi_collapsible_channel_adapt_external_results.dat 2.0 5.0e-5 >> validation.log
fi

mv RESLT_ALG RESLT_fsi_collapsible_channel_adapt_external


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
