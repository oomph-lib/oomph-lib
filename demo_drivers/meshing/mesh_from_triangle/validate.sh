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


# Validation for mesh generation from triangle with box mesh
#------------------------------------------------------------

echo "Running triangle mesh generation Poisson validation with box mesh"
mkdir RESLT
../mesh_from_triangle_poisson ../box_hole.1.node ../box_hole.1.ele \
../box_hole.1.poly > OUTPUT_box_poisson  

echo "done"
echo " " >> validation.log
echo "triangle mesh generation Poisson validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat  \
    > box_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/box_poisson_results.dat.gz   \
    box_poisson_results.dat  >> validation.log
fi

mv RESLT RESLT_box_poisson


# Validation for mesh generation from triangle with box mesh (ADAPT)
#------------------------------------------------------------

echo "Running triangle mesh generation Poisson validation with box mesh (ADAPT)"
mkdir RESLT
../mesh_from_triangle_poisson_adapt ../box_hole_adapt.1.node ../box_hole_adapt.1.ele \
../box_hole_adapt.1.poly > OUTPUT_box_poisson_adapt  

echo "done"
echo " " >> validation.log
echo "triangle mesh generation Poisson validation (ADAPT)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm0.dat RESLT/norm1.dat  \
    > box_poisson_adapt_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/box_poisson_adapt_results.dat.gz   \
    box_poisson_adapt_results.dat  >> validation.log
fi

mv RESLT RESLT_box_poisson_adapt


# Validation for mesh generation from triangle with Navier Stokes
#-----------------------------------------------------------------

echo "Running triangle mesh generation Navier Stokes validation with box mesh"
mkdir RESLT
../mesh_from_triangle_navier_stokes ../flow_past_box.1.node \
../flow_past_box.1.ele ../flow_past_box.1.poly \
 > OUTPUT_box_navier_stokes  

echo "done"
echo " " >> validation.log
echo "triangle mesh generation Navier Stokes validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat \
    > box_navier_stokes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/box_navier_stokes_results.dat.gz   \
    box_navier_stokes_results.dat  >> validation.log
fi

mv RESLT RESLT_box_navier_stokes

# Validation for mesh generation from triangle with Navier Stokes (ADAPT)
#-----------------------------------------------------------------

echo "Running triangle mesh generation Navier Stokes validation with box mesh (ADAPT)"
mkdir RESLT
../mesh_from_triangle_navier_stokes_adapt ../flow_past_box_adapt.1.node \
../flow_past_box_adapt.1.ele ../flow_past_box_adapt.1.poly \
 > OUTPUT_box_navier_stokes_adapt  

echo "done"
echo " " >> validation.log
echo "triangle mesh generation Navier Stokes validation (ADAPT)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm0.dat RESLT/norm1.dat \
    > box_navier_stokes_adapt_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/box_navier_stokes_adapt_results.dat.gz   \
    box_navier_stokes_adapt_results.dat  >> validation.log
fi

mv RESLT RESLT_box_navier_stokes_adapt



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
