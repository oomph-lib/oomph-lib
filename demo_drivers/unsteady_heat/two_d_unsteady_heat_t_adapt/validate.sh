#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for adaptive unsteady heat
#--------------------------------------
echo "Running adaptive 2D unsteady heat validation "
mkdir RESLT
../two_d_unsteady_heat_t_adapt  > OUTPUT_for_restart
echo "done run for restart"
echo " " >> validation.log
echo "2D adaptive unsteady heat validation " >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat \
    > result.dat
mv RESLT RESLT_for_restart

mkdir RESLT
../two_d_unsteady_heat_t_adapt RESLT_for_restart/restart3.dat > OUTPUT_restarted
echo "done restarted run"
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat \
    >> result.dat
mv RESLT RESLT_restarted



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
fi

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
