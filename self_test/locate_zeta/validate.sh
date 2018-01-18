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


cd Validation

echo "Running 2D locate_zeta validation "
mkdir RESLT
../locate_zeta_tester > OUTPUT_2D
echo "done"
echo " " >> validation.log
echo "2D locate zeta validation" >> validation.log
echo "-------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
cat RESLT/success_cgal.dat RESLT/success_ref_bin.dat RESLT/success_non_ref_bin.dat > results_2d.dat
../../../bin/fpdiff.py ../validata/results_2d.dat.gz results_2d.dat >> validation.log
fi

mv RESLT RESLT_2D


echo "Running 2D triangle locate_zeta validation "
mkdir RESLT
../locate_zeta_tester_triangle > OUTPUT_triangle
echo "done"
echo " " >> validation.log
echo "2D triangle locate zeta validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
cat RESLT/success_cgal.dat RESLT/success_ref_bin.dat RESLT/success_non_ref_bin.dat > results_triangle.dat
../../../bin/fpdiff.py ../validata/results_triangle.dat.gz results_triangle.dat >> validation.log
fi

mv RESLT RESLT_triangle


echo "Running 3D locate_zeta validation "
mkdir RESLT
../locate_zeta_tester_3d > OUTPUT_3D
echo "done"
echo " " >> validation.log
echo "3D locate zeta validation" >> validation.log
echo "-------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
cat RESLT/success_cgal.dat RESLT/success_ref_bin.dat RESLT/success_non_ref_bin.dat > results_3d.dat
../../../bin/fpdiff.py ../validata/results_3d.dat.gz results_3d.dat >> validation.log
fi

mv RESLT RESLT_3D


echo "Running 3D tetgen locate_zeta validation "
mkdir RESLT
cp ../cube_hole.1.node ../cube_hole.1.ele ../cube_hole.1.face .
../locate_zeta_tester_tetgen cube_hole.1.node cube_hole.1.ele cube_hole.1.face > OUTPUT_tetgen

echo "done"
echo " " >> validation.log
echo "3D tetgen locate zeta validation" >> validation.log
echo "--------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
cat RESLT/success_cgal.dat RESLT/success_ref_bin.dat RESLT/success_non_ref_bin.dat > results_tetgen.dat
../../../bin/fpdiff.py ../validata/results_tetgen.dat.gz results_tetgen.dat >> validation.log
fi

mv RESLT RESLT_tetgen

# Append log to main validation log
cat validation.log >> ../../../validation.log

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
