#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for compressed square
#---------------------------------

cd Validation



mkdir RESLT0
mkdir RESLT1
mkdir RESLT2
mkdir RESLT3
mkdir RESLT4
mkdir RESLT5
mkdir RESLT6
mkdir RESLT7
mkdir RESLT8
mkdir RESLT9
mkdir RESLT10
mkdir RESLT11
mkdir RESLT12
mkdir RESLT13
mkdir RESLT14
mkdir RESLT15



echo "Running compressed square "
../compressed_square > OUTPUT


echo "done"
echo " " >> validation.log
echo "compressed square" >> validation.log
echo "-----------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat \
RESLT0/soln1.dat \
RESLT1/soln1.dat \
RESLT2/soln1.dat \
RESLT3/soln1.dat \
RESLT4/soln1.dat \
RESLT5/soln1.dat \
RESLT6/soln1.dat \
RESLT7/soln1.dat \
RESLT8/soln1.dat \
RESLT9/soln1.dat \
RESLT10/soln1.dat \
RESLT11/soln1.dat \
RESLT12/soln1.dat \
RESLT13/soln1.dat \
RESLT14/soln1.dat \
RESLT15/soln1.dat \
    > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  0.1 1.0e-7 >> validation.log
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
