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

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running fish poisson validation "
mkdir RESLT_fully_automatic
mkdir RESLT_incremental
../fish_poisson > OUTPUT_fish_poisson
echo "done"
echo " " >> validation.log
echo "Fish poisson validation" >> validation.log
echo "-----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_fully_automatic/soln0.dat RESLT_incremental/soln3.dat \
 > fish_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_results.dat.gz   \
    fish_poisson_results.dat  >> validation.log
fi

echo "Running fish poisson validation without adaptation "
mkdir RESLT
../fish_poisson_no_adapt > OUTPUT_fish_poisson_no_adapt
echo "done"
echo " " >> validation.log
echo "Fish poisson without adaptation validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > fish_poisson_no_adapt_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_no_adapt_results.dat.gz  \
    fish_poisson_no_adapt_results.dat  >> validation.log
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
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
