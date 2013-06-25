#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for sphere scattering
#---------------------------------
cd Validation

echo "Running sphere scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../sphere_scattering > OUTPUT
echo "done"
echo " " >> validation.log
echo "Sphere validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  RESLT/soln1.dat RESLT/soln2.dat  RESLT/soln3.dat \
    RESLT/total_power0.dat  RESLT/total_power1.dat RESLT/total_power2.dat  \
    RESLT/total_power3.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
fi


mv RESLT RESLT_sphere_scattering
mv OUTPUT OUTPUT_sphere_scattering

# Validation for unstructured sphere scattering
#----------------------------------------------

echo "Running unstructured sphere scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../unstructured_sphere_scattering > OUTPUT
echo "done"
echo " " >> validation.log
echo "Unstructured sphere validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_results.dat.gz   \
    unstructured_results.dat  >> validation.log
fi


mv RESLT RESLT_unstructured_sphere_scattering
mv OUTPUT OUTPUT_unstructured_sphere_scattering




# Validation for adaptive unstructured sphere scattering
#-------------------------------------------------------

echo "Running adaptive unstructured sphere scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../unstructured_adaptive_sphere_scattering > OUTPUT
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured sphere validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_unstructured_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_results.dat.gz   \
    adaptive_unstructured_results.dat  >> validation.log
fi


mv RESLT RESLT_adaptive_unstructured_sphere_scattering
mv OUTPUT OUTPUT_unstructured_sphere_scattering




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
