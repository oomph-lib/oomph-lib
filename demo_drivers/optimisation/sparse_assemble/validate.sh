#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=0  # hierher: Tests temporarily disabled until we
             #           sort out the problems with the additional
             #           small entries. 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for different assembly strategies
#---------------------------------------------
cd Validation

echo "Running validation of different assembly strategies"
for n in 1 2 3 4 5
 do
  ../sparse_assemble_test 10 10 $n 1 0 > OUTPUT_$n
  #UNIX sort with magic to do numerical sorting on the first column
  #and then on the second column if the first columns are the same
  sort -k 1,1n -k 2,2n  matrix$n.dat  > matrix$n.dat.sorted
 done

echo "done"
echo " " >> validation.log
echo " Validation of different assembly strategies" >> validation.log
echo "--------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


echo "  " >> validation.log
echo "------------------------------------------------" >> validation.log
echo "NOTE: fpdiff is currently disabled on this test." >> validation.log
echo "=====" >> validation.log
echo "------------------------------------------------" >> validation.log
echo "  " >> validation.log

#if test "$1" = "no_fpdiff"; then
#  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
#else
#../../../../bin/fpdiff.py ../validata/matrix1.dat.gz   \
#    matrix1.dat.sorted >> validation.log
#../../../../bin/fpdiff.py ../validata/matrix2.dat.gz   \
#    matrix2.dat.sorted >> validation.log
#../../../../bin/fpdiff.py ../validata/matrix3.dat.gz   \
#    matrix3.dat.sorted >> validation.log
#../../../../bin/fpdiff.py ../validata/matrix4.dat.gz   \
#    matrix4.dat.sorted >> validation.log
#fi


# Append output to global validation log file
#--------------------------------------------
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
