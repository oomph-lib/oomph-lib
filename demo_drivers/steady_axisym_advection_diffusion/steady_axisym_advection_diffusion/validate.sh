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

# Validation for demo axisym steady advection diffusion with flux bc
#-------------------------------------------------------------------
cd Validation

echo "Running axisym steady advection diffusion with flux bc validation "
mkdir RESLT
../steady_axisym_advection_diffusion > OUTPUT
echo "done"
echo " "                                                           >> validation.log
echo "2D axisym steady advection diffusion with flux bc validation " >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " "                                                           >> validation.log
echo "Validation directory: "                                      >> validation.log
echo " "                                                           >> validation.log
echo "  " `pwd`                                                    >> validation.log
echo " "                                                           >> validation.log
cat RESLT/soln0.dat > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
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
