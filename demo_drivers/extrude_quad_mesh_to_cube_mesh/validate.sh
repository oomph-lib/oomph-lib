#!/bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# Set the number of tests to be checked
NUM_TESTS=5

# Setup validation directory
#---------------------------
touch Validation
rm -rf Validation
mkdir Validation

# Alias for the location of the fpdiff python script; note that the oomph
# root directory is one level higher because we descended into the Validation
# directory after getting the relative address for the root directory
alias fpdiff="../$OOMPH_ROOT_DIR/bin/fpdiff.py"

# Validation for extrusion of triangle-generated mesh
#-----------------------------------------------------
cd Validation
mkdir RESLT
cp -r ../triangle_meshes ./

echo "Running extrude_triangle_generated_mesh"
../extrude_triangle_generated_mesh > OUTPUT
echo "done "
echo " " >> validation.log
echo "Mesh extrusion validation" >> validation.log
echo "-------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  validata_folder="../validata/extrude_triangle_generated_mesh/"
  fpdiff $validata_folder/soln0.dat.gz ./RESLT/soln0.dat 0.1 1.0e-8 >> validation.log
  fpdiff $validata_folder/soln1.dat.gz ./RESLT/soln1.dat 0.1 1.0e-8 >> validation.log
fi

mv RESLT RESLT_triangle

# Validation for extrusion of mesh with macro-element representation
#-------------------------------------------------------------------
mkdir RESLT

echo "Running extrude_with_macro_element_representation"
../extrude_with_macro_element_representation > OUTPUT
echo "done "
echo " " >> validation.log
echo "Mesh extrusion with macro-element validation" >> validation.log
echo "--------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  validata_folder="../validata/extrude_with_macro_element_representation/"
  fpdiff $validata_folder/soln0.dat.gz ./RESLT/soln0.dat 0.1 1.0e-8 >> validation.log
  fpdiff $validata_folder/soln1.dat.gz ./RESLT/soln1.dat 0.1 1.0e-8 >> validation.log
  fpdiff $validata_folder/soln2.dat.gz ./RESLT/soln2.dat 0.1 1.0e-8 >> validation.log
fi

mv RESLT RESLT_macro_element

# Append log to main validation log
cat validation.log >> ../$OOMPH_ROOT_DIR/validation.log

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
