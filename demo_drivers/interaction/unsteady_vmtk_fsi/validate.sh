#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for vmtk fsi
#------------------------
cd Validation

cp ../fluid.1.ele .
cp ../fluid.1.face .
cp ../fluid.1.node .
cp ../fluid_quadratic_nodes.dat .
cp ../solid.1.ele .
cp ../solid.1.face .
cp ../solid.1.node .
cp ../solid_quadratic_nodes.dat .
cp ../boundary_enumeration.dat .

mkdir RESLT

../unsteady_vmtk_fsi bla > OUTPUT
echo "done"
echo " " >> validation.log
echo "Unsteady vmtk FSI (quadratic FSI surface) validation" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0.dat  RESLT/fluid_soln2.dat   \
    RESLT/solid_soln0.dat  RESLT/solid_soln2.dat   \
    > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         results.dat 0.1 1.0e-11 >> validation.log
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
