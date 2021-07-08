#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


# Bypass self-tests if executable hasn't been built (almost
# certainly because there's no cgal
if [ ! -f planar_facet_bounded_mesh_from_gmsh ]; then
    echo " "
    echo "Skipping four adaptive tet-mesh tests, presumably because"
    echo "driver codes weren't compiled because we don't have CGAL."
    echo " "
    echo " " > validation.log
    echo "Skipping four adaptive tet-mesh tests, presumably because" >> validation.log
    echo "driver codes weren't compiled because we don't have CGAL." >> validation.log
    echo " " >> validation.log
    cat validation.log >> ../../../validation.log
    exit
fi


#Set the number of tests to be checked
NUM_TESTS=4

# Have we specified command to run gmsh?
#---------------------------------------
run_gmsh=1

# Run this twice just in case Makefile.am re-generates itself and
# creates some output in the process
make -s --no-print-directory spit_out_gmsh_command
gmsh_command=`make -s --no-print-directory spit_out_gmsh_command`
if [ "$gmsh_command" = "" ]; then 
    NUM_TESTS=2    
    echo " "
    echo "Haven't specified command to run gmsh on command line."
    echo "Therefore not running gmsh-based self-tests!"
    echo " "
    run_gmsh=0
else
    echo " "
    echo "Specified command to run gmsh on command line:"
    echo " "
    echo "     "$gmsh_command
    echo " "
    rm -f .gmsh_version_number_dont_create_this_file
    `$gmsh_command --version 2>  .gmsh_version_number_dont_create_this_file`
    version=`cat .gmsh_version_number_dont_create_this_file`
    rm -f .gmsh_version_number_dont_create_this_file
    echo "gmsh version: "$version
    if [ "$version" != "3.0.6" ]; then
        echo "Wrong version number; I require 3.0.6 (for now) to be on the"
        echo "safe side. To avoid problems with self-tests, I'm not running"
        echo "them. Please investigate if your version of gmsh"
        echo "is compatible..."
        NUM_TESTS=2    
        run_gmsh=0
    else
        echo "Version of gmsh is correct."
    fi
fi

# Suppress (very costly!) bulk output -- re-enable if there are any problems...
suppress_bulk_output_flag=""
suppress_bulk_output_flag="--suppress_bulk_output"

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validation for mesh generation from tetgen
#-------------------------------------------

echo "Running adaptive tetgen mesh generation for curved boundaries"
mkdir RESLT
../curved_facet_bounded_mesh_from_tetgen $suppress_bulk_output_flag > OUTPUT_tetgen_curved
echo "done"
echo " " >> validation.log
echo "Adaptive tetgen mesh generation for curved boundaries validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/volumes2.dat RESLT/norm2.dat \
    > tetgen_curved_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/tetgen_curved_results.dat.gz \
    tetgen_curved_results.dat  >> validation.log
fi

mv RESLT RESLT_tetgen_curved


# Validation for mesh generation from gmsh
#-------------------------------------------
if [ $run_gmsh -eq 1 ]; then

echo "Running adaptive gmsh mesh generation for curved boundaries"
mkdir RESLT
../curved_facet_bounded_mesh_from_gmsh --gmsh_command_line $gmsh_command $suppress_bulk_output_flag > OUTPUT_gmsh_curved
echo "done"
echo " " >> validation.log
echo "Adaptive gmsh mesh generation for curved boundaries validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/volumes2.dat RESLT/norm2.dat \
    > gmsh_curved_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/gmsh_curved_results.dat.gz \
    gmsh_curved_results.dat  >> validation.log
fi

mv RESLT RESLT_gmsh_curved

fi


###### hierher
#if [ 1 -eq 0 ]; then


# Validation for mesh generation from tetgen
#-------------------------------------------

echo "Running adaptive tetgen mesh generation for planar boundaries"
mkdir RESLT
../planar_facet_bounded_mesh_from_tetgen $suppress_bulk_output_flag > OUTPUT_tetgen_planar
echo "done"
echo " " >> validation.log
echo "Adaptive tetgen mesh generation for planar boundaries validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/volumes2.dat RESLT/norm2.dat \
    > tetgen_planar_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/tetgen_planar_results.dat.gz \
    tetgen_planar_results.dat  >> validation.log
fi

mv RESLT RESLT_tetgen_planar


# Validation for mesh generation from gmsh
#-------------------------------------------
if [ $run_gmsh -eq 1 ]; then

echo "Running adaptive gmsh mesh generation for planar boundaries"
mkdir RESLT
../planar_facet_bounded_mesh_from_gmsh --gmsh_command_line $gmsh_command $suppress_bulk_output_flag > OUTPUT_gmsh_planar
echo "done"
echo " " >> validation.log
echo "Adaptive gmsh mesh generation for planar boundaries validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/volumes2.dat RESLT/norm2.dat \
    > gmsh_planar_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/gmsh_planar_results.dat.gz \
    gmsh_planar_results.dat  >> validation.log
fi

mv RESLT RESLT_gmsh_planar
fi

#### hierher
#fi


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
