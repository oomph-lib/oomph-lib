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



# Validation for 3D cantilever deformation problem
#-------------------------------------------------

echo "Running 3D cantilever problem "
mkdir RESLT_MESH
mkdir RESLT_refine0
mkdir RESLT_refine1
mkdir RESLT_refine2
mkdir RESLT_refine3
mkdir RESLT_refine4
mkdir RESLT_refine5
mkdir RESLT_refine6
mkdir RESLT_refine7
mkdir RESLT_refine8
mkdir RESLT_refine9

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../three_d_cantilever_adapt > OUTPUT_three_d_cantilever_adapt
echo "done"
echo " " >> validation.log
echo "3D cantilever validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_refine0/soln1_on_proc0.dat RESLT_refine0/soln1_on_proc1.dat \
    RESLT_refine1/soln1_on_proc0.dat RESLT_refine1/soln1_on_proc1.dat \
    RESLT_refine2/soln1_on_proc0.dat RESLT_refine2/soln1_on_proc1.dat \
    RESLT_refine3/soln1_on_proc0.dat RESLT_refine3/soln1_on_proc1.dat \
    RESLT_refine4/soln1_on_proc0.dat RESLT_refine4/soln1_on_proc1.dat \
    RESLT_refine5/soln1_on_proc0.dat RESLT_refine5/soln1_on_proc1.dat \
    RESLT_refine6/soln1_on_proc0.dat RESLT_refine6/soln1_on_proc1.dat \
    RESLT_refine7/soln1_on_proc0.dat RESLT_refine7/soln1_on_proc1.dat \
    RESLT_refine8/soln1_on_proc0.dat RESLT_refine8/soln1_on_proc1.dat \
    RESLT_refine9/soln1_on_proc0.dat RESLT_refine9/soln1_on_proc1.dat \
    > three_d_cantilever_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/three_d_cantilever_results.dat.gz  \
         three_d_cantilever_results.dat 0.2 4.0e-7 >> validation.log
fi

mkdir RESLT_three_d_cantilever_adapt
mv RESLT_refine* RESLT_three_d_cantilever_adapt

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
