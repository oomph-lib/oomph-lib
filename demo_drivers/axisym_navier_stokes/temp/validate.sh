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

# Validation for axisymmetric static cap
#---------------------------------------
cd Validation

echo "Running axisymmetric static cap validation "
mkdir RESLT_hijacked_external
mkdir RESLT_hijacked_internal
mkdir RESLT_elastic_hijacked_external
mkdir RESLT_elastic_hijacked_internal
#valgrind --leak-check=full -v 
../axi_static_cap > OUTPUT_axi_static_cap
echo "done"
echo " " >> validation.log
echo "Static axisymmetric cap validation" >> validation.log
echo "------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_hijacked_external/soln5.dat RESLT_hijacked_external/trace.dat \
     RESLT_hijacked_internal/soln5.dat RESLT_hijacked_internal/trace.dat \
     RESLT_elastic_hijacked_external/soln5.dat \
     RESLT_elastic_hijacked_external/trace.dat \
     RESLT_elastic_hijacked_internal/soln5.dat \
     RESLT_elastic_hijacked_internal/trace.dat \
     RESLT/soln0.vtu > axi_static_cap.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/axi_static_cap.dat.gz  \
   axi_static_cap.dat 0.1 1.0e-8 >> validation.log
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
