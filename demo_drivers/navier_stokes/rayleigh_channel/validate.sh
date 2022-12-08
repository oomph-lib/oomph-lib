#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running rayleigh_channel"
mkdir RESLT_CR RESLT_TH
../rayleigh_channel 0 0 > OUTPUT_rayleigh_channel1
echo "done (start from periodic solution)"
echo " " >> validation.log
echo "Rayleigh_channel validation with IC=periodic soln" >> validation.log
echo "-------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT_CR
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results1_CR.dat
cd ../RESLT_TH
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results1_TH.dat
cd ..
../rayleigh_channel 0 1  > OUTPUT_rayleigh_channel2
echo "done (start from impulsive start)"
echo " " >> validation.log
echo "Rayleigh_channel validation with impulsive start" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT_CR
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results2_CR.dat
cd ../RESLT_TH
cat  soln0.dat soln1.dat soln2.dat soln3.dat soln4.dat soln5.dat \
    > ../results2_TH.dat
cd ..

cat ./results1_CR.dat ./results1_TH.dat \
    ./results2_CR.dat ./results2_TH.dat \
    > results.dat


if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results.dat.gz  \
         ./results.dat 0.1 1.0e-8 >> validation.log
fi



# Append log to main validation log
cat validation.log >> $OOMPH_ROOT_DIR/validation.log

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
