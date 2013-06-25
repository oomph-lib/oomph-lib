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


# Validation for non-FSI Turek problem (algebraic node update)
#-------------------------------------------------------------
cd Validation

echo "Running non-FSI Turek problem (algebraic node update)"
mkdir RESLT

# Do validation run
../turek_flag_non_fsi_alg blabla > OUTPUT_alg
echo "done"
echo " " >> validation.log
echo "Non-FSI Turek problem (algebraic node update)" \
>> validation.log
echo "---------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat \
    > result_alg.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_alg.dat.gz \
    result_alg.dat >> validation.log
fi


mv RESLT RESLT_alg


# Validation for non-FSI Turek problem (Domain/MacroElement-based node update)
#-----------------------------------------------------------------------------

echo "Running non-FSI Turek problem (Domain-based node update)"
mkdir RESLT

# Do validation run
../turek_flag_non_fsi_macro blabla > OUTPUT_macro
echo "done"
echo " " >> validation.log
echo "Non-FSI Turek problem (Domain/MacroElement-based node update)" \
>> validation.log
echo "-------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat \
    > result_macro.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_macro.dat.gz \
    result_macro.dat >> validation.log
fi


mv RESLT RESLT_macro


# Append log to main validation log
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
