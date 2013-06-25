#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running 2D double-diffusive convection validation "
mkdir RESLT
../dd_convection lalala > ./OUTPUT_dd
echo "done"
echo " " >> validation.log
echo "2D double-diffusive convection validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/soln5.dat > dd.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/dd.dat.gz \
    dd.dat  0.1 1.0e-14 >> validation.log
fi
mv RESLT RESLT_dd


echo "Running 2D double-diffusive convection (multimesh) validation "
mkdir RESLT
../multimesh_dd_convection lalala > ./OUTPUT_dd_multimesh
echo "done"
echo " " >> validation.log
echo "2D double-diffusive convection (multimesh) validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/soln5.dat > dd_multimesh.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/dd_multimesh.dat.gz \
    dd_multimesh.dat  0.1 1.0e-12 >> validation.log
fi
mv RESLT RESLT_dd_multimesh

echo "Running Refineable 2D double-diffusive convection (multimesh) validation "
mkdir RESLT_ref_multimesh
../multimesh_ref_dd_convection lalala > ./OUTPUT_ref_dd_multimesh
echo "done"
echo " " >> validation.log
echo "Refineable 2D double-diffusive convection (multimesh) validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_ref_multimesh/trace.dat RESLT_ref_multimesh/nst_soln5.dat \
    RESLT_ref_multimesh/temp_soln5.dat RESLT_ref_multimesh/conc_soln5.dat \
     > ref_dd_multimesh.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/ref_dd_multimesh.dat.gz \
    ref_dd_multimesh.dat 0.1 2.0e-12 >> validation.log
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
