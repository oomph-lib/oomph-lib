#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for buckling of clamped shell 
#-----------------------------------------

mkdir RESLT_clamped

echo "Running clamped shell "
../clamped_or_pinned_shell 0 > OUTPUT_clamped

echo "done"
echo " " >> validation.log
echo "Clamped shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_clamped/shell0.dat RESLT_clamped/shell4.dat RESLT_clamped/shell9.dat RESLT_clamped/trace.dat > clamped_shell.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/clamped_shell.dat.gz \
 clamped_shell.dat>> validation.log
fi




# Validation for buckling of pinned shell 
#-----------------------------------------

mkdir RESLT_pinned

echo "Running pinned shell "
../clamped_or_pinned_shell 1 > OUTPUT_clamped

echo "done"
echo " " >> validation.log
echo "Pinned shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_pinned/shell0.dat RESLT_pinned/shell4.dat RESLT_pinned/shell9.dat RESLT_pinned/trace.dat > pinned_shell.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pinned_shell.dat.gz \
 pinned_shell.dat>> validation.log
fi



# Validation for buckling of pinned periodic shell 
#-----------------------------------------

mkdir RESLT_pinned_periodic

echo "Running pinned periodic shell "
../clamped_or_pinned_shell 2 > OUTPUT_clamped

echo "done"
echo " " >> validation.log
echo "Pinned periodic shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_pinned_periodic/shell0.dat RESLT_pinned_periodic/shell4.dat RESLT_pinned_periodic/shell9.dat RESLT_pinned_periodic/trace.dat > pinned_periodic.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pinned_periodic.dat.gz \
 pinned_periodic.dat>> validation.log
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
