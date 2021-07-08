#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=8

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for 2D two layer static cap
#---------------------------------------
cd Validation

echo "Running 3D static cap validation "
mkdir RESLT_TH_internal RESLT_TH_external RESLT_CR_internal RESLT_CR_external
mkdir RESLT_TH_internal_elastic RESLT_TH_external_elastic RESLT_CR_internal_elastic RESLT_CR_external_elastic
../3d_static_cap > OUTPUT_3d_static_cap
../3d_static_cap_elastic > OUTPUT_3d_static_cap_elastic
echo "done"
echo " " >> validation.log
echo "3D static cap validation" >> validation.log
echo "-------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH_internal/soln0.dat RESLT_TH_internal/soln1.dat > TH_int.dat 
cat  RESLT_TH_external/soln0.dat RESLT_TH_external/soln1.dat > TH_ext.dat 
cat  RESLT_CR_internal/soln0.dat RESLT_CR_internal/soln1.dat > CR_int.dat 
cat  RESLT_CR_external/soln0.dat RESLT_CR_external/soln1.dat > CR_ext.dat 
cat  RESLT_TH_internal_elastic/soln0.dat RESLT_TH_internal_elastic/soln1.dat \
    > TH_int_elastic.dat 
cat  RESLT_TH_external_elastic/soln0.dat RESLT_TH_external_elastic/soln1.dat \
    > TH_ext_elastic.dat 
cat  RESLT_CR_internal_elastic/soln0.dat RESLT_CR_internal_elastic/soln1.dat \
    > CR_int_elastic.dat 
cat  RESLT_CR_external_elastic/soln0.dat RESLT_CR_external_elastic/soln1.dat \
    > CR_ext_elastic.dat 


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Spine Tests" >> validation.log
echo >> validation.log
echo "Taylor Hood (hijacked internal pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/TH_int.dat.gz  \
   TH_int.dat 0.1 1.0e-5 >> validation.log
echo "Taylor Hood (hijacked external pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/TH_ext.dat.gz  \
   TH_ext.dat 0.1 1.1e-5 >> validation.log
echo "Crouzeix Raviart (hijacked internal pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/CR_int.dat.gz  \
   CR_int.dat 0.1 1.0e-5 >> validation.log
echo "Crouzeix Raviart (hijacked external pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/CR_ext.dat.gz  \
   CR_ext.dat 0.1 1.0e-5 >> validation.log

echo >> validation.log
echo "PseudoSolidMeshUpdate Tests" >> validation.log
echo >> validation.log
echo "Taylor Hood (hijacked internal pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/TH_int_elastic.dat.gz  \
   TH_int_elastic.dat 0.1 5.0e-7 >> validation.log
echo "Taylor Hood (hijacked external pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/TH_ext_elastic.dat.gz  \
   TH_ext_elastic.dat 0.1 5.0e-7 >> validation.log
echo "Crouzeix Raviart (hijacked internal pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/CR_int_elastic.dat.gz  \
   CR_int_elastic.dat 0.1 5.0e-7 >> validation.log
echo "Crouzeix Raviart (hijacked external pressure)">> validation.log
../../../../bin/fpdiff.py  ../validata/CR_ext_elastic.dat.gz  \
   CR_ext_elastic.dat 0.1 5.0e-7 >> validation.log

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
