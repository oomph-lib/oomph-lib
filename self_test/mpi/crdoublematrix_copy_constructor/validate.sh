#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1

# Doc what we're using to run tests on variable processors
echo " " 
echo "Running mpi tests with mpi run command: " $MPI_VARIABLENP_RUN_COMMAND
echo "OOMPHNP = 1, 2, 3, 4"
echo " " 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validation for matrix addition test
#-----------------------------------------

echo "Running CRDoubleMatrix copy constructor test "

# one processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/1/g"`
$MPI_RUN_ON_NP_COMMAND ../crdoublematrix_copy_constructor

# two processors
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/2/g"`
$MPI_RUN_ON_NP_COMMAND ../crdoublematrix_copy_constructor

# three processors
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/3/g"`
$MPI_RUN_ON_NP_COMMAND ../crdoublematrix_copy_constructor

# four processors
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/4/g"`
$MPI_RUN_ON_NP_COMMAND ../crdoublematrix_copy_constructor


echo "done"
echo " " >> validation.log
echo "Distributed CRDoubleMatrix copy constructor test" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Cat all the outputs into one file.
cat \
t1m1_NP1R0 \
copy_of_t1m1_NP1R0 \
t1m1_NP2R0 \
copy_of_t1m1_NP2R0 \
t1m1_NP2R1 \
copy_of_t1m1_NP2R1 \
t1m1_NP3R0 \
copy_of_t1m1_NP3R0 \
t1m1_NP3R1 \
copy_of_t1m1_NP3R1 \
t1m1_NP3R2 \
copy_of_t1m1_NP3R2 \
t1m1_NP4R0 \
copy_of_t1m1_NP4R0 \
t1m1_NP4R1 \
copy_of_t1m1_NP4R1 \
t1m1_NP4R2 \
copy_of_t1m1_NP4R2 \
t1m1_NP4R3 \
copy_of_t1m1_NP4R3 \
t2m1_NP1R0 \
copy_of_t2m1_NP1R0 \
t2m1_NP2R0 \
copy_of_t2m1_NP2R0 \
t2m1_NP2R1 \
copy_of_t2m1_NP2R1 \
t2m1_NP3R0 \
copy_of_t2m1_NP3R0 \
t2m1_NP3R1 \
copy_of_t2m1_NP3R1 \
t2m1_NP3R2 \
copy_of_t2m1_NP3R2 \
t2m1_NP4R0 \
copy_of_t2m1_NP4R0 \
t2m1_NP4R1 \
copy_of_t2m1_NP4R1 \
t2m1_NP4R2 \
copy_of_t2m1_NP4R2 \
t2m1_NP4R3 \
copy_of_t2m1_NP4R3 \
t3m1_NP1R0 \
copy_of_t3m1_NP1R0 \
t3m1_NP2R0 \
copy_of_t3m1_NP2R0 \
t3m1_NP2R1 \
copy_of_t3m1_NP2R1 \
t3m1_NP3R0 \
copy_of_t3m1_NP3R0 \
t3m1_NP3R1 \
copy_of_t3m1_NP3R1 \
t3m1_NP3R2 \
copy_of_t3m1_NP3R2 \
t3m1_NP4R0 \
copy_of_t3m1_NP4R0 \
t3m1_NP4R1 \
copy_of_t3m1_NP4R1 \
t3m1_NP4R2 \
copy_of_t3m1_NP4R2 \
t3m1_NP4R3 \
copy_of_t3m1_NP4R3 \
> crdoublematrix_copy_constructor.dat

# Remove the individual outputs.
rm -rf \
t1m1_NP1R0 \
copy_of_t1m1_NP1R0 \
t1m1_NP2R0 \
copy_of_t1m1_NP2R0 \
t1m1_NP2R1 \
copy_of_t1m1_NP2R1 \
t1m1_NP3R0 \
copy_of_t1m1_NP3R0 \
t1m1_NP3R1 \
copy_of_t1m1_NP3R1 \
t1m1_NP3R2 \
copy_of_t1m1_NP3R2 \
t1m1_NP4R0 \
copy_of_t1m1_NP4R0 \
t1m1_NP4R1 \
copy_of_t1m1_NP4R1 \
t1m1_NP4R2 \
copy_of_t1m1_NP4R2 \
t1m1_NP4R3 \
copy_of_t1m1_NP4R3 \
t2m1_NP1R0 \
copy_of_t2m1_NP1R0 \
t2m1_NP2R0 \
copy_of_t2m1_NP2R0 \
t2m1_NP2R1 \
copy_of_t2m1_NP2R1 \
t2m1_NP3R0 \
copy_of_t2m1_NP3R0 \
t2m1_NP3R1 \
copy_of_t2m1_NP3R1 \
t2m1_NP3R2 \
copy_of_t2m1_NP3R2 \
t2m1_NP4R0 \
copy_of_t2m1_NP4R0 \
t2m1_NP4R1 \
copy_of_t2m1_NP4R1 \
t2m1_NP4R2 \
copy_of_t2m1_NP4R2 \
t2m1_NP4R3 \
copy_of_t2m1_NP4R3 \
t3m1_NP1R0 \
copy_of_t3m1_NP1R0 \
t3m1_NP2R0 \
copy_of_t3m1_NP2R0 \
t3m1_NP2R1 \
copy_of_t3m1_NP2R1 \
t3m1_NP3R0 \
copy_of_t3m1_NP3R0 \
t3m1_NP3R1 \
copy_of_t3m1_NP3R1 \
t3m1_NP3R2 \
copy_of_t3m1_NP3R2 \
t3m1_NP4R0 \
copy_of_t3m1_NP4R0 \
t3m1_NP4R1 \
copy_of_t3m1_NP4R1 \
t3m1_NP4R2 \
copy_of_t3m1_NP4R2 \
t3m1_NP4R3 \
copy_of_t3m1_NP4R3

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/crdoublematrix_copy_constructor.dat.gz  \
         crdoublematrix_copy_constructor.dat >> validation.log
fi

#----------------------------------------------------------------------

# Append log to main validation log
cat validation.log >> ../../../../validation.log

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
