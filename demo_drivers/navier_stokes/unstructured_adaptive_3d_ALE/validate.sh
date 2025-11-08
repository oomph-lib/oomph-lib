#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

# Bypass self-tests if executable hasn't been built (almost
# certainly because there's no cgal
if [ ! -f uns_adapt_3d ]; then
    echo " "
    echo "Skipping unstructured adaptive 3D Navier Stokes ALE validation,"
    echo "presumably because driver code wasn't compiled because we don't"
    echo "have CGAL."
    echo " "
    echo " " >validation.log
    echo "Skipping unstructured adaptive 3D Navier Stokes ALE validation," >>validation.log
    echo "presumably because driver code wasn't compiled because we don't" >>validation.log
    echo "have CGAL." >>validation.log
    echo " " >>validation.log
    cat validation.log >>$OOMPH_ROOT_DIR/validation.log
    exit
fi

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for unstructured fluid
#----------------------------------
cd Validation

echo "Running unstructured adaptive 3D Navier Stokes ALE validation "
mkdir RESLT

../uns_adapt_3d --validation >OUTPUT
echo "done"
echo " " >>validation.log
echo "Unstructured adaptive 3D Navier Stokes ALE validation" >>validation.log
echo "------------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT/diss.dat >results.dat

if test "$2" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
    $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results.dat.gz \
        results.dat 2.3 1.0e-14 >>validation.log
fi

# Append log to main validation log
cat validation.log >>$OOMPH_ROOT_DIR/validation.log

cd ..

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 10
