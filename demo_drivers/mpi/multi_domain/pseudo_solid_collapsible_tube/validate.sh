#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

if [ ! -f pseudo_solid_collapsible_tube ]; then
    echo " "
    echo "Pseudo-elastic collapsible tube validation "
    echo " "
    echo " SKIPPED BECAUSE WE DON'T HAVE TRILINOS OR HYPRE"
    echo " "
    echo " " >> validation.log
    echo "Pseudo-elastic collapsible tube validation" >> validation.log
    echo "------------------------------------------" >> validation.log
    echo " SKIPPED BECAUSE WE DON'T HAVE TRILINOS OR HYPRE" >> validation.log
    exit
fi


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

# Validation for pseudo-elastic collapsible tube
#-----------------------------------------------

echo "Running pseudo-elastic collapsible tube validation "
mkdir RESLT_proc0
mkdir RESLT_proc1


# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../pseudo_solid_collapsible_tube --use_iterative_solver --validate  > OUTPUT
echo "done"
echo " " >> validation.log
echo "Pseudo-elastic collapsible tube validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_proc0/solid_soln3.dat  RESLT_proc0/solid_soln4.dat \
     RESLT_proc0/fluid_soln3.dat  RESLT_proc0/fluid_soln4.dat \
     RESLT_proc1/solid_soln3.dat  RESLT_proc1/solid_soln4.dat \
     RESLT_proc1/fluid_soln3.dat  RESLT_proc1/fluid_soln4.dat \
     > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         results.dat 6.0 3.0e-5 >> validation.log
fi

# check iteration counts -- allow for 10% difference
grep iterations RESLT_proc1/OUTPUT.1 | awk '{print $5}' > iter_counts.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/iter_counts.dat.gz  \
         iter_counts.dat 10.0 1.0e-12 >> validation.log
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
