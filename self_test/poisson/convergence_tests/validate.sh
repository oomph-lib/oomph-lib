#! /bin/sh


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


if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
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


if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
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


if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
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


if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
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
OK_COUNT=`grep -c 'OK' Validation/validation.log`
if  [ $OK_COUNT -eq $NUM_TESTS ]; then
 echo " "
 echo "======================================================================"
 echo " " 
 echo "All tests in" 
 echo " " 
 echo "    `pwd`    "
 echo " "
 echo "passed successfully."
 echo " "
 echo "======================================================================"
 echo " " 
else
  if [ $OK_COUNT -lt $NUM_TESTS ]; then
   echo " "
   echo "======================================================================"
   echo " " 
   echo "Only $OK_COUNT of $NUM_TESTS test(s) passed; see"
   echo " " 
   echo "    `pwd`/Validation/validation.log"
   echo " " 
   echo "for details" 
   echo " " 
   echo "======================================================================"
   echo " "
  else 
   echo " "
   echo "======================================================================"
   echo " " 
   echo "More OKs than tests! Need to update NUM_TESTS in"
   echo " " 
   echo "    `pwd`/validate.sh"
   echo " "
   echo "======================================================================"
   echo " "
  fi
fi
