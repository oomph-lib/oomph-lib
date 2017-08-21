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

# Validation for rectangular driven cavity
#-----------------------------------------
cd Validation

echo "Running rectangular driven cavity validation with subtracted singularity"
mkdir RESLT
../driven_cavity > OUTPUT_driven_cavity
echo "done"
echo " " >> validation.log
echo "Rectangular driven cavity with subtracted singularity validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

#omit the possibly non-portable lines with infs in 
awk '{if ((NR!=2)&&(NR!=920)){print $0}}' RESLT/soln0.dat \
> driven_cavity_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/../bin/fpdiff.py ../validata/driven_cavity_results.dat.gz  \
         driven_cavity_results.dat >> validation.log
fi

mv RESLT RESLT_driven_cavity





echo "Running circular couette validation with subtracted pseudo-singularity"
mkdir RESLT
../circular_couette --re 10.0 > OUTPUT_circular_couette
echo "done"
echo " " >> validation.log
echo "Circular couette validation with subtracted pseudo-singularity validation" >> validation.log
echo "-------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  > circular_couette_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/../bin/fpdiff.py ../validata/circular_couette_results.dat.gz  \
         circular_couette_results.dat >> validation.log
fi

mv RESLT RESLT_circular_couette





# Append log to main validation log
cat validation.log >> $OOMPH_ROOT_DIR/../validation.log

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
