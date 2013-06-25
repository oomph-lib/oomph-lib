#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=7


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running Boussinesq Convection validation with FD for off-diagonals"
mkdir RESLT
../boussinesq_convection lalala > OUTPUT_boussinesq_convection
echo "done"
echo " " >> validation.log
echo "Boussinesq Convection validation with FD for off-diagonals" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > bous_convection_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    bous_convection_results.dat 0.1 1.0e-8  >> validation.log
fi

mv RESLT RESLT_non_refineable



echo "Running Boussinesq Convection validation with FD for entire Jacobian"
mkdir RESLT
../boussinesq_convection_fd lalala > OUTPUT_boussinesq_convection_FD
echo "done"
echo " " >> validation.log
echo "Boussinesq Convection validation with FD for entire Jacobian" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > bous_convection_fd_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    bous_convection_fd_results.dat 0.1 5.0e-7  >> validation.log
fi

mv RESLT RESLT_non_refineable_fd


echo "Running Boussinesq Convection validation with analytic Jacobian"
mkdir RESLT
../boussinesq_convection_analytic lalala > OUTPUT_boussinesq_convection_analytic
echo "done"
echo " " >> validation.log
echo "Boussinesq Convection validation with analytic Jacobian" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > bous_convection_anal_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    bous_convection_anal_results.dat 0.1 5.0e-7  >> validation.log
fi

mv RESLT RESLT_non_refineable_anal


echo "Running Refineable Boussinesq Convection validation "
mkdir RESLT
../refineable_b_convection > OUTPUT_refineable_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq Convection validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat > ref_bous_convection_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/ref_results.dat.gz   \
   ref_bous_convection_results.dat  0.1 1.0e-7 >> validation.log
fi

mv RESLT RESLT_refineable


# Validation for Boussinesq convection problem using multi-domain method
#-----------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method) "
mkdir RESLT
../multi_domain_boussinesq_convection validate > OUTPUT_multi_domain_b_convection
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0.dat RESLT/fluid_soln5.dat \
    RESLT/temperature_soln0.dat RESLT/temperature_soln5.dat \
    > multi_domain_b_convection.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/multi_domain_b_convection.dat.gz  \
         multi_domain_b_convection.dat 0.1 1.5e-7 >> validation.log
fi

mv RESLT RESLT_multi_domain_boussinesq_convection


# Validation for refineable Boussinesq convection problem, multi-domain
#----------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (multi-domain method) "
mkdir RESLT
../multi_domain_ref_b_convection > OUTPUT_multi_domain_ref_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (multi-domain) validation" >> validation.log
echo "----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0.dat RESLT/fluid_soln1.dat \
    RESLT/temperature_soln0.dat RESLT/temperature_soln1.dat \
    > multi_domain_ref_b_convection.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/multi_domain_ref_b_convection.dat.gz  \
         multi_domain_ref_b_convection.dat 0.1 1.0e-7 >> validation.log
fi

mv RESLT RESLT_multi_domain_ref_b_convection



# Validation for refineable Boussinesq convection problem, multi-domain, FD for external data
#--------------------------------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (multi-domain method; FD for external data) "
mkdir RESLT
../multi_domain_ref_b_convection_fd_for_external_data  > OUTPUT_multi_domain_ref_b_convection_fd_for_external_data 
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (multi-domain; FD for external data) validation" >> validation.log
echo "--------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0.dat RESLT/fluid_soln1.dat \
    RESLT/temperature_soln0.dat RESLT/temperature_soln1.dat \
    > multi_domain_ref_b_convection_fd_for_external_data.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/multi_domain_ref_b_convection.dat.gz  \
         multi_domain_ref_b_convection_fd_for_external_data.dat 0.1 1.0e-7 >> validation.log
fi

mv RESLT RESLT_multi_domain_ref_b_convection_fd_for_external_data 


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
