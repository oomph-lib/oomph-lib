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

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running Two-layer soluble surfactant validation "
mkdir RESLT
../two_layer_soluble_surfactant lalala > ./OUTPUT_two_layer_sol_surf
echo "done"
echo " " >> validation.log
echo "Two-layer soluble surfactant validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/int5.dat RESLT/soln5.dat > 2layer_sol_surf.dat
mv RESLT RESLT_2layer

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2layer_sol_surf.dat.gz \
    2layer_sol_surf.dat  0.1 1.0e-14 >> validation.log
fi


echo "Running Refineable two-layer soluble surfactant validation "
mkdir RESLT
../refineable_two_layer_soluble_surfactant lalala > ./OUTPUT_ref_two_layer_sol_surf
echo "done"
echo " " >> validation.log
echo "Refineable two-layer soluble surfactant validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/int5.dat RESLT/soln5.dat > ref_2layer_sol_surf.dat
mv RESLT RESLT_ref_2layer

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/ref_2layer_sol_surf.dat.gz \
    ref_2layer_sol_surf.dat  0.25 1.0e-14 >> validation.log
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
