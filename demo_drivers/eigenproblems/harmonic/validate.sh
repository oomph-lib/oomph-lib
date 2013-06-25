#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running harmonic equation eigenproblem validation "
mkdir RESLT
cd RESLT
../../harmonic > ../OUTPUT_harmonic
cd ..
echo "done"
echo " " >> validation.log
echo "Harmonic eigenproblem validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/eigenvalues1.dat RESLT/soln1.dat > res_arpack.dat
cat RESLT/eigenvalues2.dat RESLT/soln2.dat > res_qz.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "ARPACK Test: " >> validation.log
../../../../bin/fpdiff.py  ../validata/harmonic_results.dat.gz   \
    res_arpack.dat  0.1 1.0e-13 >> validation.log
echo "QZ Test: " >> validation.log
../../../../bin/fpdiff.py  ../validata/harmonic_results.dat.gz   \
    res_qz.dat  0.1 1.0e-13 >> validation.log
fi


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
MAKEFLAGS=s
OOMPH_ROOT_DIR=$(make print-top_builddir)
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
