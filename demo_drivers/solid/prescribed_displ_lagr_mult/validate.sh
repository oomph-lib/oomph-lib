#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation imposed boundary deformation using Lagrange multipliers
#-------------------------------------------------------------------

cd Validation

mkdir RESLT

echo "Running imposed boundary deformation using Lagrange multipliers "
../prescribed_displ_lagr_mult > OUTPUT_with_lagr_mult

echo "done"
echo " " >> validation.log
echo "Imposed boundary deformation using Lagrange multipliers" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/lagr2.dat RESLT/soln2.dat > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result.dat  >> validation.log
fi

mv RESLT RESLT_with_lagr_mult

mkdir RESLT

echo "Running imposed boundary deformation without Lagrange multipliers "
../prescribed_displ_lagr_mult2 > OUTPUT_without_lagr_mult

echo "done"
echo " " >> validation.log
echo "Imposed boundary deformation without Lagrange multipliers" >> validation.log
echo "---------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln2.dat > result2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result2.dat.gz \
    result2.dat  >> validation.log
fi

mv RESLT RESLT_without_lagr_mult



mkdir RESLT

echo "Running precond imposed boundary deformation without Lagrange multipliers "
../prescribed_displ_lagr_mult_precond --block_upper_for_elastic_block > OUTPUT_without_lagr_mult_precond

echo "done"
echo " " >> validation.log
echo "Preconditioned imposed boundary deformation without Lagrange multipliers" >> validation.log
echo "------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

max_iter=`grep iterations OUTPUT_without_lagr_mult_precond | awk 'BEGIN{max=0}{if ($NF>max) max=$NF}END{print max}'`
echo " " >> validation.log
echo "Max number of GMRES iterations: "$max_iter>> validation.log
echo " " >> validation.log
if [ $max_iter -le 40 ]; then
      echo "   [OK] -- for number of iterations below 40" >> validation.log
else
      echo "   [FAILED] -- number of iterations above 40" >> validation.log
fi

cat RESLT/soln2.dat > result_precond.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_precond.dat.gz \
    result_precond.dat  >> validation.log
fi

mv RESLT RESLT_precond


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
