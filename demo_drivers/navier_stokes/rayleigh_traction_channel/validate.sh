#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


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

echo "Running rayleigh_traction_channel "
mkdir RESLT_CR RESLT_TH
../rayleigh_traction_channel 0 0 > OUTPUT_rayleigh_traction_channel1
echo "done (start from periodic solution)"
echo " " >> validation.log
echo "Rayleigh traction channel validation with IC=periodic soln" >> validation.log
echo "----------------------------------------------------------" >> validation.log
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
../rayleigh_traction_channel 0 1 > OUTPUT_rayleigh_traction_channel2
echo "done (start from impulsive start)"
echo " " >> validation.log
echo "Rayleigh traction channel validation with impulsive start" >> validation.log
echo "---------------------------------------------------------" >> validation.log
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


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz \
         ./results.dat 0.1 1.0e-8 >> validation.log
fi



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
