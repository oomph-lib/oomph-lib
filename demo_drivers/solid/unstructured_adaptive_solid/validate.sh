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

# Validation for static fish deformation
#---------------------------------------

cd Validation
mkdir RESLT RESLT_pres_disp RESLT_pres_disp_incomp

echo "Running unstructured adaptive solid validation "
../unstructured_adaptive_solid > OUTPUT_adaptive_solid


echo "done"
echo " " >> validation.log
echo "Unstructured adaptive solid validation" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/s_energy.dat > res.dat
cat RESLT_pres_disp/s_energy.dat > res_pres_disp.dat
cat RESLT_pres_disp_incomp/s_energy.dat > res_pres_disp_incomp.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/res.dat.gz \
    res.dat 0.1 1.0e-7 >> validation.log
../../../../bin/fpdiff.py ../validata/res_pres_disp.dat.gz \
    res_pres_disp.dat 0.1 1.0e-7 >> validation.log
../../../../bin/fpdiff.py ../validata/res_pres_disp_incomp.dat.gz \
    res_pres_disp_incomp.dat 0.11 1.5e-7 >> validation.log
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
