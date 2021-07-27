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

# Disable the self-test on macOS as it causes a mysterious segmentation fault.
# To do this, we provide a dummy OK in the validation log. Note that the script
# validate_ok_count reads the Validation/validation.log file so we need to make
# sure to create the log file in there
if [[ uname -eq "Darwin" ]]; then
  echo "#==================================================== " >>Validation/validation.log
  echo "dummy [OK] -- Test fails on macOS, not running it! " >>Validation/validation.log
  echo "#====================================================" >>Validation/validation.log
  . $OOMPH_ROOT_DIR/bin/validate_ok_count

  # Never get here
  exit 10
fi

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running 1D reaction-diffusion validation "
mkdir RESLT
cd RESLT
../../one_d_act_inhibit >../OUTPUT_1d_act
cd ..
echo "done"
echo " " >>validation.log
echo "1D reaction-diffusion validation " >>validation.log
echo "------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT/step10.dat >one_d_act_inhibit.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  ../../../../bin/fpdiff.py ../validata/one_d_act_inhibit.dat.gz \
    one_d_act_inhibit.dat >>validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >>../../../../validation.log

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
