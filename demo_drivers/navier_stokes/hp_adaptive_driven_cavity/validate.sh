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

cd Validation

# Validation for hp-adaptive rectangular driven cavity
#-----------------------------------------------------

echo "Running hp-adaptive rectangular driven cavity validation "
mkdir RESLT
../hp_adaptive_driven_cavity validate > OUTPUT
echo "done"
echo " " >> validation.log
echo "hp-adaptive rectangular driven cavity validation" >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         results.dat >> validation.log
fi

# Validation for hp-adaptive circular driven cavity
#--------------------------------------------------

echo "Running hp-adaptive circular driven cavity validation "
mkdir RESLT_circ_bdry0 RESLT_circ_bdry1
../circular_driven_cavity_hp_adapt > OUTPUT_circ
echo "done"
echo " " >> validation.log
echo "hp-adaptive circular driven cavity validation" >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat  RESLT_circ_bdry0/soln0.dat RESLT_circ_bdry0/soln1.dat > results_circ_bdry0.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_circ_bdry0.dat.gz  \
         results_circ_bdry0.dat >> validation.log
fi

cat  RESLT_circ_bdry1/soln0.dat RESLT_circ_bdry1/soln1.dat > results_circ_bdry1.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_circ_bdry1.dat.gz  \
         results_circ_bdry1.dat >> validation.log
fi

# Append log to main validation log
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
