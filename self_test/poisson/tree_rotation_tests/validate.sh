#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=78

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for 2D Rotations
#--------------------------------------
cd Validation

echo "Running 2D Tree rotation validation "
mkdir RESLT1 RESLT2 RESLT3
../tree_2d > OUTPUT_2d_rotation_validation
echo "done"
echo " " >> validation.log
echo " 2D Tree rotation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo " LINEAR ELEMENTS: " >> validation.log
echo >> validation.log
echo " 90 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT1/soln0.dat RESLT1/soln1.dat >> validation.log
echo " 180 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT1/soln0.dat RESLT1/soln2.dat >> validation.log
echo " 270 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT1/soln0.dat RESLT1/soln3.dat >> validation.log

echo " QUADRATIC ELEMENTS:" >> validation.log
echo >> validation.log
echo " 90 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT2/soln0.dat RESLT2/soln1.dat >> validation.log
echo " 180 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT2/soln0.dat RESLT2/soln2.dat >> validation.log
echo " 270 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT2/soln0.dat RESLT2/soln3.dat >> validation.log

echo " CUBIC SPECTRAL ELEMENTS: " >> validation.log
echo >> validation.log
echo " 90 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT3/soln0.dat RESLT3/soln1.dat >> validation.log
echo " 180 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT3/soln0.dat RESLT3/soln2.dat >> validation.log
echo " 270 degree rotation " >> validation.log
../../../../bin/fpdiff.py RESLT3/soln0.dat RESLT3/soln3.dat >> validation.log

fi

mv RESLT1 RESLT_2d_linear
mv RESLT2 RESLT_2d_quad
mv RESLT3 RESLT_2d_spec_cubic

# Validation for 3D rotations
#--------------------------------

echo "Running 3D Tree rotation validation "
mkdir RESLT1 RESLT2 RESLT3
../tree_3d > OUTPUT_3d_rotation_validation
echo "done"
echo " " >> validation.log
echo " 3D Tree rotation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo " LINEAR ELEMENTS: " >> validation.log
echo >> validation.log
for INDEX in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
  do
  echo "Rotation Number " $INDEX >> validation.log
  ../../../../bin/fpdiff.py RESLT1/soln0.dat RESLT1/soln${INDEX}.dat >> validation.log
done

echo " QUADRATIC ELEMENTS: " >> validation.log
echo >> validation.log
for INDEX in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
  do
  echo "Rotation Number " $INDEX >> validation.log
  ../../../../bin/fpdiff.py RESLT2/soln0.dat RESLT2/soln${INDEX}.dat >> validation.log
done

echo " CUBIC SPECTRAL ELEMENTS: " >> validation.log
echo >> validation.log
for INDEX in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
  do
  echo "Rotation Number " $INDEX >> validation.log
  ../../../../bin/fpdiff.py RESLT3/soln0.dat RESLT3/soln${INDEX}.dat >> validation.log
done

fi


mv RESLT1 RESLT_3d_linear
mv RESLT2 RESLT_3d_quad
mv RESLT3 RESLT_3d_spec_cubic

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
