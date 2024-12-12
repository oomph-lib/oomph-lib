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

cd Validation


# Validate mesh generation from inline triangle (internal boundaries)
#---------------------------------------------------------------------------

echo "Running inline triangle mesh generation test (internal boundaries)"
mkdir RESLT
../mesh_from_inline_triangle_internal_boundaries --validation > OUTPUT

echo "done"
echo " " >> validation.log
echo "triangle inline mesh generation test (internal boundaries)" >> validation.log
echo "-----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  \
    > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat 2 1.0e-14 >> validation.log
fi

mv RESLT RESLT_internal_boundaries


# Validate mesh generation from inline triangle_extra (internal boundaries and definition of
# regions and holes)
#---------------------------------------------------------------------------

echo "Running inline triangle mesh generation test (internal boundaries -- definition of holes and regions --)"
mkdir RESLT
../mesh_from_inline_triangle_internal_boundaries_extra --validation > OUTPUT

echo "done"
echo " " >> validation.log
echo "triangle inline mesh generation test (internal boundaries -- definition of holes and regions --)" >> validation.log
echo "-----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  \
    > results_extra.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_extra.dat.gz   \
    results_extra.dat  2 1.0e-14 >> validation.log
fi

mv RESLT RESLT_internal_boundaries_extra


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
