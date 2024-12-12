#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# Set the number of tests to be checked
NUM_TESTS=5

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo driver
#---------------------------
cd Validation

echo "Running 2D vorticity smoother validation "
mkdir RESLT
../vorticity_smoother_validation > OUTPUT
echo "done"
echo " " >> validation.log
echo "2D vorticity smoother validation" >> validation.log
echo "--------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/vorticity_convergence.dat > convergence.dat
cat RESLT/soln0.dat > results0.dat
cat RESLT/soln1.dat > results1.dat
cat RESLT/soln2.dat > results2.dat
cat RESLT/soln3.dat > results3.dat

if test "$1" = "no_fpdiff"
then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/results0.dat.gz \
				     results0.dat 1.0 1.0e-10 >> validation.log
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/results1.dat.gz \
				     results1.dat 1.0 1.0e-08 >> validation.log
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/results2.dat.gz \
				     results2.dat 1.0 3.0e-08 >> validation.log
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/results3.dat.gz \
				     results3.dat 2.0 3.0e-07 >> validation.log
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/convergence.dat.gz \
				     convergence.dat >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log

cd ..

#######################################################################

# Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
