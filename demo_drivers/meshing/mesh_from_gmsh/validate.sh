#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for Airy cantilever
#-------------------------------

cd Validation

echo "Running a 3D Gmsh version of Airy cantilever validation "
mkdir RESLT
../mesh_from_gmsh_solid ../cantilever.msh > OUTPUT


echo "done"
echo " " >> validation.log
echo "3D Gmsh Airy cantilever validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > result.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz result.dat >> validation.log
fi

mv RESLT RESLT_mesh_gmsh



# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log
