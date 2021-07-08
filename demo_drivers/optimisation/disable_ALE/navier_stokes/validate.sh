#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Rayleigh channel
#--------------------------------
cd Validation

echo "Running rayleigh channel with/without ALE"
mkdir RESLT_CR RESLT_TH RESLT_CR_ALE RESLT_TH_ALE
../rayleigh_channel 0 0 > OUTPUT_rayleigh_channel
echo "done"
echo " " >> validation.log
echo "Rayleigh_channel with/without ALE" >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT_CR
cat  soln0.dat soln1.dat > ../results_CR.dat
cd ../RESLT_TH
cat  soln0.dat soln1.dat > ../results_TH.dat
cd ../RESLT_CR_ALE
cat  soln0.dat soln1.dat > ../results_CR_ALE.dat
cd ../RESLT_TH_ALE
cat  soln0.dat soln1.dat > ../results_TH_ALE.dat
cd ..

cat ./results_CR.dat ./results_TH.dat \
    ./results_CR_ALE.dat ./results_TH_ALE.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/results.dat.gz  results.dat 0.1 1.0e-8 >> validation.log
fi


# Validation for Rayleigh driven cavity
#--------------------------------------

echo "Running rayleigh driven cavity with/without ALE"
mkdir RESLT2_CR RESLT2_TH RESLT2_CR_ALE RESLT2_TH_ALE
../rayleigh_circular_driven_cavity_adapt 0 0 > OUTPUT_rayleigh_driven_cavity
echo "done"
echo " " >> validation.log
echo "Rayleigh driven cavity with/without ALE" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT2_CR
cat  soln0.dat soln1.dat > ../results2_CR.dat
cd ../RESLT2_TH
cat  soln0.dat soln1.dat > ../results2_TH.dat
cd ../RESLT2_CR_ALE
cat  soln0.dat soln1.dat > ../results2_CR_ALE.dat
cd ../RESLT2_TH_ALE
cat  soln0.dat soln1.dat > ../results2_TH_ALE.dat
cd ..

cat ./results2_CR.dat ./results2_TH.dat \
    ./results2_CR_ALE.dat ./results2_TH_ALE.dat > results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/results2.dat.gz  results2.dat 0.1 1.0e-8 >> validation.log
fi



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
