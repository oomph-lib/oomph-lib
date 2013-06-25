#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for steady buckling ring with displacement control
#--------------------------------------------------------------

cd Validation
mkdir RESLT_global
mkdir RESLT_no_global

echo "Running steady ring validation "
../steady_ring > OUTPUT_steady_ring

echo "done"
echo " " >> validation.log
echo "Steady ring validation" >> validation.log
echo "----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_global/ring0.dat \
    RESLT_global/ring5.dat \
    RESLT_global/ring12.dat \
    RESLT_global/ring20.dat \
    RESLT_global/trace.dat\
    > ring_results.dat

cat RESLT_no_global/ring0.dat \
    RESLT_no_global/ring5.dat \
    RESLT_no_global/ring12.dat \
    RESLT_no_global/ring20.dat \
    RESLT_no_global/trace.dat\
    > ring_results2.dat
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../../bin/fpdiff.py ../validata/ring_results.dat.gz \
   ring_results.dat >> validation.log
  ../../../../bin/fpdiff.py ../validata/ring_results.dat.gz \
   ring_results2.dat >> validation.log
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
