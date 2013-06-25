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

# Validation for small amplitude oscillating ring
#------------------------------------------------

cd Validation
mkdir RESLT

echo "Running small amplitude oscillating ring validation "
../lin_unsteady_ring 1 0 0 > OUTPUT_lin_unsteady_ring

echo "done"
echo " " >> validation.log
echo "Small amplitude oscillating ring validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/ring1.dat \
    RESLT/ring3.dat \
    RESLT/ring5.dat \
    RESLT/ring9.dat \
    RESLT/trace_ring.dat \
    > lin_unsteady_results.dat
mv RESLT RESLT_lin_unsteady_ring

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../../bin/fpdiff.py ../validata/lin_unsteady_results.dat.gz \
   lin_unsteady_results.dat >> validation.log
fi

#######################################################################

# Validation for large amplitude oscillating ring
#------------------------------------------------

mkdir RESLT

echo "Running large amplitude oscillating ring validation "
../unsteady_ring 1 > OUTPUT_unsteady_ring

echo "done"
echo " " >> validation.log
echo "Large amplitude oscillating ring validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/ring1.dat \
    RESLT/ring3.dat \
    RESLT/ring5.dat \
    RESLT/ring9.dat \
    RESLT/trace_ring.dat \
    >  unsteady_results.dat
mv RESLT RESLT_unsteady_ring


mkdir RESLT

echo "Running restarted large amplitude oscillating ring validation "
../unsteady_ring 1 RESLT_unsteady_ring/ring_restart5.dat \
 > OUTPUT_restarted_unsteady_ring

echo "done"
echo " " >> validation.log
echo "Restarted large amplitude oscillating ring validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/ring1.dat \
    RESLT/ring3.dat \
    RESLT/ring5.dat \
    RESLT/ring9.dat \
    RESLT/trace_ring.dat \
    >  restarted_unsteady_results.dat
mv RESLT RESLT_restarted_unsteady_ring



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../../bin/fpdiff.py ../validata/unsteady_results.dat.gz \
   unsteady_results.dat >> validation.log
  ../../../../bin/fpdiff.py ../validata/restarted_unsteady_results.dat.gz \
   restarted_unsteady_results.dat >> validation.log
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
