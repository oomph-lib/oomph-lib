#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validation for mesh generation from tetgen with cube hole mesh
#----------------------------------------------------------------

echo "Running tetgen mesh generation Poisson validation with cube_hole"
mkdir RESLT
../mesh_from_tetgen_poisson ../cube_hole.1.node ../cube_hole.1.ele \
../cube_hole.1.face > OUTPUT_tetgen_poisson
echo "done"
echo " " >> validation.log
echo "tetgen mesh generation Poisson validation" >> validation.log
echo "-----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln1.dat RESLT/soln3.dat \
    > cube_hole_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/cube_hole_poisson_results.dat.gz \
    cube_hole_poisson_results.dat  >> validation.log
fi

mv RESLT RESLT_poisson





# Validation for Navier Stokes from tetgen
#-----------------------------------------

echo "Running tetgen Navier Stokes validation"
mkdir RESLT
../mesh_from_tetgen_navier_stokes ../cube_hole.1.node ../cube_hole.1.ele \
../cube_hole.1.face > OUTPUT_tetgen_navier_stokes
echo "done"
echo " " >> validation.log
echo "tetgen mesh generation Navier Stokes validation" >> validation.log
echo "-----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > cube_hole_navier_stokes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/cube_hole_navier_stokes_results.dat.gz \
    cube_hole_navier_stokes_results.dat >> validation.log
fi

mv RESLT RESLT_navier_stokes


# Append output to global validation log file
#--------------------------------------------
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
 exit 0
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
   exit 1
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
  exit 2
  fi
fi




# Never get here
exit 10
