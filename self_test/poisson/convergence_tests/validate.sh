#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4



# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation



# Validation for 2D TPoisson convergence
#--------------------------------------
cd Validation

echo "Running 2D TPoisson convergence validation "
mkdir RESLT
../t_convergence_2d lala > OUTPUT_t_convergence_2d
echo "done"
echo " " >> validation.log
echo " 2D TPoisson convergence validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat RESLT/convergence.dat > t_convergence_2d_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/t_convergence_2d_results.dat.gz  \
         t_convergence_2d_results.dat >> validation.log
fi

mv RESLT RESLT_t_convergence_2d



# Validation for 2D QPoisson convergence
#--------------------------------------

echo "Running 2D QPoisson convergence validation "
mkdir RESLT
../q_convergence_2d lala > OUTPUT_q_convergence_2d
echo "done"
echo " " >> validation.log
echo " 2D QPoisson convergence validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat RESLT/convergence.dat > q_convergence_2d_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/q_convergence_2d_results.dat.gz  \
         q_convergence_2d_results.dat >> validation.log
fi

mv RESLT RESLT_q_convergence_2d


# Validation for 3D QPoisson convergence
#--------------------------------------

echo "Running 3D QPoisson convergence validation "
mkdir RESLT
../q_convergence_3d lala > OUTPUT_q_convergence_3d
echo "done"
echo " " >> validation.log
echo " 3D QPoisson convergence validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat RESLT/convergence.dat > q_convergence_3d_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/q_convergence_3d_results.dat.gz  \
         q_convergence_3d_results.dat >> validation.log
fi

mv RESLT RESLT_q_convergence_3d


# Validation for 3D TPoisson convergence
#--------------------------------------

echo "Running 3D TPoisson convergence validation "
mkdir RESLT
../t_convergence_3d lala > OUTPUT_t_convergence_3d
echo "done"
echo " " >> validation.log
echo " 3D TPoisson convergence validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat RESLT/soln3.dat RESLT/convergence.dat > t_convergence_3d_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/t_convergence_3d_results.dat.gz  \
         t_convergence_3d_results.dat >> validation.log
fi

mv RESLT RESLT_t_convergence_3d




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
