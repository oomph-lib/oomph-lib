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

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running Marangoni Convection validation with periodic boundaries"
mkdir RESLT
../marangoni_convection lalala > OUTPUT_marangoni_convection
echo "done"
echo " " >> validation.log
echo "Marangoni Convection validation with periodic boundaries" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > mar_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/mar_results.dat.gz   \
    mar_results.dat 0.1 1.0e-14  >> validation.log
fi

mv RESLT RESLT_mar


echo "Running Marangoni Convection with rigid boundaries and contact angle"
mkdir RESLT
../marangoni_convection_box lalala > OUTPUT_marangoni_convection_box
echo "done"
echo " " >> validation.log
echo "Marangoni Convection validation with rigid boundaries" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln3_75.dat RESLT/soln5_75.dat > mar_box_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/mar_box_results.dat.gz   \
    mar_box_results.dat 0.1 1.0e-14  >> validation.log
fi

mv RESLT RESLT_mar_box

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
