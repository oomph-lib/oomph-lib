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

cd Validation

# Validation for buckling of clamped shell with arclength continuation
#---------------------------------------------------------------------

mkdir RESLT
cd RESLT

echo "Running clamped_shell with arclength continuation validation "
../../clamped_shell_with_arclength_cont >../OUTPUT_clamped_shell_with_arclength_cont

echo "done"
cd ..
echo " " >>validation.log
echo "Clamped shell buckling with arclength continuation validation" \
  >>validation.log
echo "---------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT/final_shape.dat RESLT/trace.dat RESLT/trace_disp.dat >shell_with_arclength_cont_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  ../../../../bin/fpdiff.py ../validata/shell_with_arclength_cont_results.dat.gz \
    shell_with_arclength_cont_results.dat 0.1 1.0e-9 >>validation.log
fi

mv RESLT RESLT_with_arclength

# Validation for buckling of clamped shell with displacement control
#-------------------------------------------------------------------

mkdir RESLT
cd RESLT

echo "Running clamped_shell with displacement control validation "
../../clamped_shell >../OUTPUT_clamped_shell

echo "done"
cd ..
echo " " >>validation.log
echo "Clamped shell buckling validation" >>validation.log
echo "---------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT/final_shape.dat RESLT/trace.dat >shell.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  ../../../../bin/fpdiff.py ../validata/shell.dat.gz \
    shell.dat 0.1 1.0e-9 >>validation.log
fi

mv RESLT RESLT_shell

# Append output to global validation log file
#--------------------------------------------
cat validation.log >>../../../../validation.log

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
