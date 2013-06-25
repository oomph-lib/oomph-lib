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


# Validation for demo unsteady heat with and without ALE
#-------------------------------------------------------
cd Validation

echo "Running 2D unsteady heat validation with and without ALE"
mkdir RESLT
mkdir RESLT_ALE
../two_d_unsteady_heat  > OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady heat validation with and without ALE" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT_ALE/soln0.dat RESLT_ALE/soln1.dat \
    > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
fi



# Validation for demo adaptive unsteady heat with and without ALE
#-----------------------------------------------------------------

echo "Running 2D adaptive unsteady heat validation with and without ALE"
mkdir RESLT_adapt
mkdir RESLT_adapt_ALE
../two_d_unsteady_heat_adapt  > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "2D adaptive unsteady heat validation with and without ALE" >> validation.log
echo "---------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_adapt/soln0.dat RESLT_adapt/soln1.dat \
    RESLT_adapt_ALE/soln0.dat RESLT_adapt_ALE/soln1.dat \
    > result_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/result_adapt.dat.gz \
    result_adapt.dat  >> validation.log
fi






# Append output to global validation log file
#--------------------------------------------
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
