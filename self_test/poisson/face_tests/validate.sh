#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4



# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation



# Validation for 2D TPoisson faces
#--------------------------------------
cd Validation

echo "Running 2D TPoisson faces validation "
mkdir RESLT_T2d
cd RESLT_T2d
../../t_faces_2d > OUTPUT_t_faces_2d
cd ..
echo "done"
echo " " >> validation.log
echo " 2D TPoisson faces validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_T2d/linear_errors_0.dat  RESLT_T2d/linear_errors_1.dat RESLT_T2d/quadratic_errors_0.dat RESLT_T2d/quadratic_errors_1.dat RESLT_T2d/cubic_errors_0.dat RESLT_T2d/cubic_errors_1.dat > ./T2d.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/T2d.dat.gz  \
        T2d.dat >> validation.log
fi



# Validation for 2D QPoisson faces
#--------------------------------------

echo "Running 2D QPoisson faces validation "
mkdir RESLT_Q2d
cd RESLT_Q2d
../../q_faces_2d > OUTPUT_q_faces_2d
cd ..
echo "done"
echo " " >> validation.log
echo " 2D QPoisson faces validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_Q2d/linear_errors.dat  RESLT_Q2d/quadratic_errors.dat RESLT_Q2d/cubic_errors.dat > ./Q2d.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/Q2d.dat.gz  \
         Q2d.dat >> validation.log
fi


# Validation for 3D QPoisson faces
#--------------------------------------

echo "Running 3D QPoisson faces validation "
mkdir RESLT_Q3d
cd RESLT_Q3d
../../q_faces_3d  > OUTPUT_q_faces_3d
cd ..
echo "done"
echo " " >> validation.log
echo " 3D QPoisson faces validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_Q3d/linear_errors.dat  RESLT_Q3d/quadratic_errors.dat RESLT_Q3d/cubic_errors.dat > ./Q3d.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/Q3d.dat.gz  \
         Q3d.dat >> validation.log
fi



# Validation for 3D TPoisson faces
#--------------------------------------

echo "Running 3D TPoisson faces validation "
mkdir RESLT_T3d
cd RESLT_T3d
../../t_faces_3d  > OUTPUT_t_faces_3d
cd ..
echo "done"
echo " " >> validation.log
echo " 3D TPoisson faces validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_T3d/quadratic_normals3_0.dat RESLT_T3d/quadratic_normals0_1.dat RESLT_T3d/quadratic_normals1_2.dat RESLT_T3d/quadratic_normals2_3.dat > ./T3d.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/T3d.dat.gz  \
         T3d.dat >> validation.log
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
