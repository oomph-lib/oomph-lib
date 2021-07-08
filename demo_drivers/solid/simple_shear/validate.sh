#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=6

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for static fish deformation
#---------------------------------------

cd Validation
mkdir RESLT RESLT_ref
mkdir RESLT_pres RESLT_pres_ref
mkdir RESLT_cont_pres RESLT_cont_pres_ref

echo "Running static simple shear deformation validation "
../simple_shear > OUTPUT_simple_shear


echo "done"
echo " " >> validation.log
echo "Static simple shear deformation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln1.dat RESLT/stress1.dat > res.dat
cat RESLT_pres/soln1.dat RESLT_pres/stress1.dat > res_pres.dat
cat RESLT_cont_pres/soln1.dat RESLT_cont_pres/stress1.dat > res_cont_pres.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/res.dat.gz \
    res.dat 0.1 2.0e-10 >> validation.log
../../../../bin/fpdiff.py ../validata/res_pres.dat.gz \
    res_pres.dat 0.1 3.5e-10 >> validation.log
../../../../bin/fpdiff.py ../validata/res_cont_pres.dat.gz \
    res_cont_pres.dat 0.1 5.0e-8 >> validation.log
fi



echo "Running refineable static simple shear deformation validation "
../refineable_simple_shear > OUTPUT_refineable_simple_shear


echo "done"
echo " " >> validation.log
echo "Refineable static simple shear deformation validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_ref/soln1.dat RESLT_ref/stress1.dat > res_ref.dat
cat RESLT_pres_ref/soln1.dat RESLT_pres_ref/stress1.dat > res_pres_ref.dat
cat RESLT_cont_pres_ref/soln1.dat RESLT_cont_pres_ref/stress1.dat \
 > res_cont_pres_ref.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/res_ref.dat.gz \
    res_ref.dat 0.1 4.0e-10 >> validation.log
../../../../bin/fpdiff.py ../validata/res_pres_ref.dat.gz \
    res_pres_ref.dat 0.1 4.0e-10 >> validation.log
../../../../bin/fpdiff.py ../validata/res_cont_pres_ref.dat.gz \
    res_cont_pres_ref.dat 0.1 5.0e-8 >> validation.log
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
