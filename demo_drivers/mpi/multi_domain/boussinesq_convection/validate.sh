#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=5

# Doc what we're using to run tests on two processors
echo " " 
echo "Running mpi tests with mpi run command: " $MPI_RUN_COMMAND
echo " " 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation
cp ../*partition*.dat .

# Validation for Boussinesq convection problem using multi-domain method
#-----------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method, analytic, connected partitioning) "
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 15

$MPI_RUN_COMMAND ../multi_domain_boussinesq_convection validate 2 > OUTPUT_multi_domain_boussinesq_convection_2
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain, analytic, connected partitioning) validation" >> validation.log
echo "---------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0_on_proc0.dat RESLT/fluid_soln1_on_proc1.dat RESLT/fluid_soln2_on_proc0.dat \
    RESLT/fluid_soln3_on_proc1.dat RESLT/fluid_soln4_on_proc0.dat RESLT/fluid_soln5_on_proc1.dat \
    RESLT/temperature_soln0_on_proc0.dat RESLT/temperature_soln1_on_proc1.dat RESLT/temperature_soln2_on_proc0.dat \
    RESLT/temperature_soln3_on_proc1.dat RESLT/temperature_soln4_on_proc0.dat RESLT/temperature_soln5_on_proc1.dat \
    > multi_domain_boussinesq_convection_analytic_2_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multi_domain_boussinesq_convection_analytic_2_results.dat.gz  \
         multi_domain_boussinesq_convection_analytic_2_results.dat 0.1 1.5e-7 >> validation.log
fi

mkdir RESLT_multi_domain_boussinesq_convection_analytic_2
mv RESLT/*dat RESLT_multi_domain_boussinesq_convection_analytic_2
rm -rf RESLT

# Validation for Boussinesq convection problem using multi-domain method
#-----------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method, analytic, one mesh per processor) "
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 15

ls -l RESLT

$MPI_RUN_COMMAND ../multi_domain_boussinesq_convection validate 1 > OUTPUT_multi_domain_boussinesq_convection
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain, analytic, one mesh per processor) validation" >> validation.log
echo "---------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0_on_proc0.dat RESLT/fluid_soln1_on_proc1.dat RESLT/fluid_soln2_on_proc0.dat \
    RESLT/fluid_soln3_on_proc1.dat RESLT/fluid_soln4_on_proc0.dat RESLT/fluid_soln5_on_proc1.dat \
    RESLT/temperature_soln0_on_proc0.dat RESLT/temperature_soln1_on_proc1.dat RESLT/temperature_soln2_on_proc0.dat \
    RESLT/temperature_soln3_on_proc1.dat RESLT/temperature_soln4_on_proc0.dat RESLT/temperature_soln5_on_proc1.dat \
    > multi_domain_boussinesq_convection_analytic_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multi_domain_boussinesq_convection_analytic_results.dat.gz  \
         multi_domain_boussinesq_convection_analytic_results.dat 0.1 1.5e-7 >> validation.log
fi

mv RESLT RESLT_multi_domain_boussinesq_convection_analytic

# Validation for Boussinesq convection problem using multi-domain method (FD)
#----------------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method, FD) "
mkdir RESLT_FD

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 15

$MPI_RUN_COMMAND ../multi_domain_boussinesq_convection_fd validate 1 > OUTPUT_multi_domain_boussinesq_convection_fd
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain, FD) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_FD/fluid_soln0_on_proc0.dat RESLT_FD/fluid_soln1_on_proc1.dat \
    RESLT_FD/fluid_soln2_on_proc0.dat RESLT_FD/fluid_soln3_on_proc1.dat \
    RESLT_FD/fluid_soln4_on_proc0.dat RESLT_FD/fluid_soln5_on_proc1.dat \
    RESLT_FD/temperature_soln0_on_proc0.dat RESLT_FD/temperature_soln1_on_proc1.dat \
    RESLT_FD/temperature_soln2_on_proc0.dat RESLT_FD/temperature_soln3_on_proc1.dat \
    RESLT_FD/temperature_soln4_on_proc0.dat RESLT_FD/temperature_soln5_on_proc1.dat \
    > multi_domain_boussinesq_convection_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multi_domain_boussinesq_convection_results.dat.gz  \
         multi_domain_boussinesq_convection_results.dat 0.1 5e-7 >> validation.log
fi

mv RESLT_FD RESLT_multi_domain_boussinesq_convection_fd


#-----------------------------------------------------------------------
# Validation for refineable Boussinesq convection problem, single domain
#-----------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (dist. single domain) "
mkdir RESLT_SINGLE


# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND ../refineable_b_convection validate > OUTPUT_refineable_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (single domain) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_SINGLE/soln0_on_proc0.dat RESLT_SINGLE/soln0_on_proc1.dat \
    RESLT_SINGLE/soln1_on_proc0.dat RESLT_SINGLE/soln1_on_proc1.dat \
    > refineable_b_convection_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/refineable_b_convection_results.dat.gz  \
         refineable_b_convection_results.dat 0.1 1.0e-7 >> validation.log
fi

mv RESLT_SINGLE RESLT_refineable_b_convection

#----------------------------------------------------------------------

# Validation for refineable Boussinesq convection problem, multi-domain
#----------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (multi-domain method) "
mkdir RESLT_MULTI

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_COMMAND  ../multi_domain_ref_b_convection validate > OUTPUT_multi_domain_ref_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (multi-domain) validation" >> validation.log
echo "---------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_MULTI/fluid_soln0_on_proc0.dat RESLT_MULTI/fluid_soln0_on_proc1.dat \
    RESLT_MULTI/fluid_soln1_on_proc0.dat RESLT_MULTI/fluid_soln1_on_proc1.dat \
    RESLT_MULTI/temperature_soln0_on_proc0.dat RESLT_MULTI/temperature_soln0_on_proc1.dat \
    RESLT_MULTI/temperature_soln1_on_proc0.dat RESLT_MULTI/temperature_soln1_on_proc1.dat \
    RESLT_MULTI/temperature_soln2_on_proc0.dat RESLT_MULTI/temperature_soln2_on_proc1.dat \
    RESLT_MULTI/temperature_soln3_on_proc0.dat RESLT_MULTI/temperature_soln3_on_proc1.dat \
    RESLT_MULTI/temperature_soln4_on_proc0.dat RESLT_MULTI/temperature_soln4_on_proc1.dat \
    RESLT_MULTI/temperature_soln5_on_proc0.dat RESLT_MULTI/temperature_soln5_on_proc1.dat \
    > multi_domain_ref_b_convection_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multi_domain_ref_b_convection_results.dat.gz  \
         multi_domain_ref_b_convection_results.dat 0.1 1.0e-7 >> validation.log
fi

mv RESLT_MULTI RESLT_multi_domain_ref_b_convection


# Append log to main validation log
cat validation.log >> ../../../../../validation.log

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
