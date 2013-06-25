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

# Validation for 2D two layer static cap
#---------------------------------------
cd Validation

echo "Running 2D static two layer validation "
mkdir RESLT RESLT_elastic
../static_two_layer > OUTPUT_static_two_layer
echo "done"
echo " " >> validation.log
echo "Two-fluid static cap validation" >> validation.log
echo "-------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln5.dat RESLT/trace.dat > static_two_layer.dat
cat  RESLT_elastic/soln5.dat RESLT_elastic/trace.dat >> static_two_layer.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/static_two_layer.dat.gz  \
   static_two_layer.dat 0.1 1.0e-8 >> validation.log
fi

mv RESLT RESLT_static_two_layer
mv RESLT_elastic RESLT_elastic_static_two_layer





# Validation for 2D single layer static cap
#---------------------------------------
echo "Running 2D static single layer validation "
mkdir RESLT_hijacked_external
mkdir RESLT_hijacked_internal
mkdir RESLT_elastic_hijacked_external
mkdir RESLT_elastic_hijacked_internal
../static_single_layer > OUTPUT_static_single_layer
echo "done"
echo " " >> validation.log
echo "Static single layer validation" >> validation.log
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
     > static_single_layer.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/static_single_layer.dat.gz  \
   static_single_layer.dat 0.1 1.0e-8 >> validation.log
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
