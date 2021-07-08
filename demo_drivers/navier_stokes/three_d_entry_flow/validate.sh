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

# Validation for 3D entry flow Navier Stokes problem (quarter tube)
#-------------------------------------------------------------------
cd Validation

echo "Running 3D entry flow (quarter tube) Navier Stokes validation "
mkdir RESLT_TH
mkdir RESLT_CR

../three_d_entry_flow lalala > OUTPUT
echo "done"
echo " " >> validation.log
echo "3D entry flow (quarter tube) Navier Stokes validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln1.dat  RESLT_CR/soln1.dat > entry_flow_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/entry_flow_results.dat.gz  \
         entry_flow_results.dat >> validation.log
fi

mv RESLT_TH RESLT_TH_quarter
mv RESLT_CR RESLT_CR_quarter

# Validation for 3D entry flow (full tube) Navier Stokes problem
#--------------------------------------------------------------------

echo "Running 3D entry flow (full tube) Navier Stokes validation "
mkdir RESLT_TH

../full_tube > OUTPUT_full
echo "done"
echo " " >> validation.log
echo "3D entry flow (full tube) Navier Stokes validation" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln_Re50.dat  > full_tube.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/full_tube.dat.gz  \
         full_tube.dat >> validation.log
fi

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
