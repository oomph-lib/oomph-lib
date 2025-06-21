#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for unstructured fluid
#----------------------------------
cd Validation

echo "Running 2D unstructured Navier Stokes validation "

# Get triangle files
cp ../fluid.fig.1.* .

mkdir RESLT_CR RESLT_TH

../unstructured_two_d_fluid >OUTPUT
echo "done"
echo " " >>validation.log
echo "2D unstructured Navier Stokes validation" >>validation.log
echo "----------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log

echo "Taylor Hood" >>validation.log
cat RESLT_TH/soln0.dat RESLT_TH/soln1.dat RESLT_TH/soln2.dat >results_TH.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results.dat.gz \
    results_TH.dat >>validation.log
fi

echo "Crouzeix Raviart" >>validation.log
cat RESLT_CR/soln0.dat RESLT_CR/soln1.dat RESLT_CR/soln2.dat >results_CR.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results_CR.dat.gz \
    results_CR.dat 0.2 1.0e-11 >>validation.log
fi

# Append log to main validation log
cat validation.log >>$OOMPH_ROOT_DIR/validation.log

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
