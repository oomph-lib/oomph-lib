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

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running adaptive axisymmetric advection diffusion in pipe  validation "
mkdir RESLT
../pipe > OUTPUT_pipe
echo "done"
echo " " >> validation.log
echo "Adaptive axisymmetric advection diffusion in pipe validation " >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln3.dat > pipe.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pipe.dat.gz \
    pipe.dat  >> validation.log
fi

mv RESLT RESLT_pipe

echo "Running adaptive axisymmetric advection variable diffusion validation "
mkdir RESLT
../pipe_variable_diff > OUTPUT_pipe_variable_diff
echo "done"
echo " " >> validation.log
echo "Adaptive axisymmetric advection variable diffusion validation " >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln3.dat > pipe_var_diff.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pipe_var_diff.dat.gz \
    pipe_var_diff.dat  >> validation.log
fi

mv RESLT RESLT_pipe_var_diff


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
