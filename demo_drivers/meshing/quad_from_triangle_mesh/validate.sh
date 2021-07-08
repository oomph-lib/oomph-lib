#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


# Set the number of tests to be checked
NUM_TESTS=11


# Set up validation directory and jump into it
#---------------------------------------------------------------------
touch Validation
rm -rf Validation
mkdir Validation
cd Validation

#-------------------------------------
# Validation for 2D moving block
#-------------------------------------

# Get triangle files
cp ../box_hole.1.* .

#---------------------------------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D moving block with Crouzeix-Raviart elements"

# Run the driver code
../adaptive_moving_block_navier_stokes > OUTPUT_2D_moving_block

# Update the validation log
echo "done"
echo " " >> validation.log
echo "2D moving block problem with Crouzeix-Raviart elements" >> validation.log
echo "------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/soln0.dat > adaptive_moving_block_navier_stokes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_moving_block_navier_stokes_results.dat.gz   \
    adaptive_moving_block_navier_stokes_results.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_moving_block


#-------------------------------------
# Validation for 2D unstructured solid
#-------------------------------------

# Get triangle files
cp ../solid.1.* .

#---------------------------------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D unstructured solid with pure displacement formulation"

# Run the driver code
../unstructured_two_d_solid --case 0 > OUTPUT_2D_pure_disp

# Update the validation log
echo "done"
echo " " >> validation.log
echo "2D unstructured solid with pure displacement formulation" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/trace.dat > unstructured_two_d_solid_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_two_d_solid_results0.dat.gz   \
    unstructured_two_d_solid_results0.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_pure_disp

#---------------------------------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D unstructured solid with pressure/displacement formulation"

# Run the driver code
../unstructured_two_d_solid --case 1 > OUTPUT_2D_pres_disp

# Update the validation log
echo "done"
echo " " >> validation.log
echo "2D unstructured solid with pressure/displacement formulation" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/trace.dat > unstructured_two_d_solid_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_two_d_solid_results1.dat.gz   \
    unstructured_two_d_solid_results1.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_pres_disp

#---------------------------------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D unstructured solid with pressure/displacement formulation (incompressible)"

# Run the driver code
../unstructured_two_d_solid --case 2 > OUTPUT_2D_pres_disp_incomp

# Update the validation log
echo "done"
echo " " >> validation.log
echo "2D unstructured solid with pressure/displacement formulation (incompressible)" >> validation.log
echo "-----------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/trace.dat > unstructured_two_d_solid_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_two_d_solid_results2.dat.gz   \
    unstructured_two_d_solid_results2.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_pres_disp_incomp


#----------------------------------------------
# Validation for 2D adaptive unstructured solid
#----------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D adaptive unstructured solid with pure displacement formulation"

# Run the driver code
../adaptive_unstructured_two_d_solid --case 0 > OUTPUT_2D_adaptive_pure_disp

# Update the validation log
echo "done"
echo " " >> validation.log
echo "Adaptive 2D unstructured solid with pure displacement formulation" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/trace.dat > adaptive_unstructured_two_d_solid_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_two_d_solid_results0.dat.gz   \
    adaptive_unstructured_two_d_solid_results0.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_adaptive_pure_disp

#---------------------------------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D adaptive unstructured solid with pressure/displacement formulation"

# Run the driver code
../adaptive_unstructured_two_d_solid --case 1 > OUTPUT_2D_adaptive_pres_disp

# Update the validation log
echo "done"
echo " " >> validation.log
echo "Adaptive 2D unstructured solid with pressure/displacement formulation" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/trace.dat > adaptive_unstructured_two_d_solid_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_two_d_solid_results1.dat.gz   \
   adaptive_unstructured_two_d_solid_results1.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_adaptive_pres_disp

#---------------------------------------------------------------------

# Make the result directory
mkdir RESLT

# Notify user
echo "Running validation for 2D adaptive unstructured solid with pressure/displacement formulation (incompressible)"

# Run the driver code
../adaptive_unstructured_two_d_solid --case 2 > OUTPUT_2D_adaptive_pres_disp_incomp

# Update the validation log
echo "done"
echo " " >> validation.log
echo "Adaptive 2D unstructured solid with pressure/displacement formulation (incompressible)" >> validation.log
echo "-----------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# Move the results to be checked
cat RESLT/trace.dat > adaptive_unstructured_two_d_solid_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_two_d_solid_results2.dat.gz   \
    adaptive_unstructured_two_d_solid_results2.dat 2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_adaptive_pres_disp_incomp


#-------------------------------------------------
# Validation for unstructured adaptive scattering
#-------------------------------------------------

echo "Running adaptive unstructured scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../adaptive_unstructured_scattering_quad --case 0 > OUTPUT_unstructured_adapt_0
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_unstructured_scattering_quad_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results0.dat.gz   \
    adaptive_unstructured_scattering_quad_results0.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_0


#---------------------------------------------------------------------


echo "Running adaptive unstructured scattering validation. first order abc"
mkdir RESLT
../adaptive_unstructured_scattering_quad --case 1 > OUTPUT_unstructured_adapt_1
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (first order ABC)" >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_unstructured_scattering_quad_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results1.dat.gz   \
    adaptive_unstructured_scattering_quad_results1.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_1


#---------------------------------------------------------------------


echo "Running unstructured adaptive scattering validation. second order abc"
mkdir RESLT
../adaptive_unstructured_scattering_quad --case 2 > OUTPUT_unstructured_adapt_2
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (second order ABC)" >> validation.log
echo "--------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_unstructured_scattering_quad_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results2.dat.gz   \
    adaptive_unstructured_scattering_quad_results2.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_2


#---------------------------------------------------------------------


echo "Running unstructured adaptive scattering validation. third order abc"
mkdir RESLT
../adaptive_unstructured_scattering_quad --case 3 > OUTPUT_adapt_3
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (third order ABC)" >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_unstructured_scattering_quad_results3.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results3.dat.gz   \
    adaptive_unstructured_scattering_quad_results3.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_3


#---------------------------------------------------------------------


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log

cd ..

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
