#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=16


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for adaptive scattering
#-----------------------------------
cd Validation

echo "Running adaptive scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../adaptive_scattering --case 0 > OUTPUT_adapt_0
echo "done"
echo " " >> validation.log
echo "Adpative cattering validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > adaptive_scattering_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_scattering_results0.dat.gz   \
    adaptive_scattering_results0.dat  >> validation.log
fi
mv RESLT RESLT_adapt_0




echo "Running adaptive scattering validation. first order abc"
mkdir RESLT
../adaptive_scattering --case 1 > OUTPUT_adapt_1
echo "done"
echo " " >> validation.log
echo "Adaptive scattering validation (first order ABC)" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > adaptive_scattering_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_scattering_results1.dat.gz   \
    adaptive_scattering_results1.dat  >> validation.log
fi
mv RESLT RESLT_adapt_1





echo "Running adaptive scattering validation. second order abc"
mkdir RESLT
../adaptive_scattering --case 2 > OUTPUT_adapt_2
echo "done"
echo " " >> validation.log
echo "Adapative scattering validation (second order ABC)" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > adaptive_scattering_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_scattering_results2.dat.gz   \
    adaptive_scattering_results2.dat  >> validation.log
fi
mv RESLT RESLT_adapt_2




echo "Running adaptive scattering validation. third order abc"
mkdir RESLT
../adaptive_scattering --case 3 > OUTPUT_adapt_3
echo "done"
echo " " >> validation.log
echo "Adaptive scattering validation (third order ABC)" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > adaptive_scattering_results3.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_scattering_results3.dat.gz   \
    adaptive_scattering_results3.dat  >> validation.log
fi
mv RESLT RESLT_adapt_3


# Validation for scattering
#--------------------------

echo "Running  scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../scattering --case 0 > OUTPUT_0
echo "done"
echo " " >> validation.log
echo "Scattering  (Dirichlet-to-Neumann BC)" >> validation.log
echo "-------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > scattering_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/scattering_results0.dat.gz   \
    scattering_results0.dat  >> validation.log
fi
mv RESLT RESLT_0




echo "Running  scattering validation. first order abc"
mkdir RESLT
../scattering --case 1 > OUTPUT_1
echo "done"
echo " " >> validation.log
echo "Scattering validation (first order ABC)" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > scattering_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/scattering_results1.dat.gz   \
    scattering_results1.dat  >> validation.log
fi
mv RESLT RESLT_1





echo "Running  scattering validation. second order abc"
mkdir RESLT
../scattering --case 2 > OUTPUT_2
echo "done"
echo " " >> validation.log
echo "Scattering validation (second order ABC)" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > scattering_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/scattering_results2.dat.gz   \
    scattering_results2.dat  >> validation.log
fi
mv RESLT RESLT_2




echo "Running  scattering validation. third order abc"
mkdir RESLT
../scattering --case 3 > OUTPUT_3
echo "done"
echo " " >> validation.log
echo "Scattering validation (third order ABC)" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > scattering_results3.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/scattering_results3.dat.gz   \
    scattering_results3.dat  >> validation.log
fi
mv RESLT RESLT_3





#########################################################################


# Validation for unstructured adaptive scattering
#-------------------------------------------------

echo "Running adaptive unstructured scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../adaptive_unstructured_scattering --case 0 > OUTPUT_unstructured_adapt_0
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_adaptive_scattering_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results0.dat.gz   \
    unstructured_adaptive_scattering_results0.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_0



echo "Running adaptive unstructured scattering validation. first order abc"
mkdir RESLT
../adaptive_unstructured_scattering --case 1 > OUTPUT_unstructured_adapt_1
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (first order ABC)" >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_adaptive_scattering_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results1.dat.gz   \
    unstructured_adaptive_scattering_results1.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_1





echo "Running unstructured adaptive scattering validation. second order abc"
mkdir RESLT
../adaptive_unstructured_scattering --case 2 > OUTPUT_unstructured_adapt_2
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (second order ABC)" >> validation.log
echo "--------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_adaptive_scattering_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results2.dat.gz   \
    unstructured_adaptive_scattering_results2.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_2




echo "Running unstructured adaptive scattering validation. third order abc"
mkdir RESLT
../adaptive_unstructured_scattering --case 3 > OUTPUT_adapt_3
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive scattering validation (third order ABC)" >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_adaptive_scattering_results3.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_scattering_results3.dat.gz   \
    unstructured_adaptive_scattering_results3.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_adapt_3


# Validation for unstructured scattering
#---------------------------------------

echo "Running unstructured scattering validation. Dirichlet-to-Neumann BC"
mkdir RESLT
../unstructured_scattering --case 0 > OUTPUT_unstructured_0
echo "done"
echo " " >> validation.log
echo "Unstructured scattering validation (Dirichlet-to-Neumann BC)" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_scattering_results0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_scattering_results0.dat.gz   \
    unstructured_scattering_results0.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_0



echo "Running unstructured scattering validation. first order abc"
mkdir RESLT
../unstructured_scattering --case 1 > OUTPUT_unstructured_1
echo "done"
echo " " >> validation.log
echo "Unstructured scattering validation (first order ABC)" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_scattering_results1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_scattering_results1.dat.gz   \
    unstructured_scattering_results1.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_1





echo "Running unstructured scattering validation. second order abc"
mkdir RESLT
../unstructured_scattering --case 2 > OUTPUT_unstructured_2
echo "done"
echo " " >> validation.log
echo "Unstructured scattering validation (second order ABC)" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_scattering_results2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_scattering_results2.dat.gz   \
    unstructured_scattering_results2.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_2




echo "Running unstructured scattering validation. third order abc"
mkdir RESLT
../unstructured_scattering --case 3 > OUTPUT_unstructured_3
echo "done"
echo " " >> validation.log
echo "Unstructured scattering validation (third order ABC)" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > unstructured_scattering_results3.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_scattering_results3.dat.gz   \
    unstructured_scattering_results3.dat  >> validation.log
fi
mv RESLT RESLT_unstructured_3





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
