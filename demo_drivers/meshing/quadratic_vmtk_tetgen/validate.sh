#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=6


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

ln -s ../short_coarse_mesh_files


# Validation for quadratic mesh smoothing
#----------------------------------------


mkdir RESLT_fluid
mkdir RESLT_solid
mkdir RESLT_fluid_with_poisson
mkdir RESLT_solid_with_poisson
mkdir RESLT_solid_with_linear_elasticity
mkdir RESLT_fluid_with_linear_elasticity

echo "Running quadratic mesh smoothing"
../snap_mesh bla > OUTPUT
echo "done"


cp RESLT_fluid/mesh_after_smooth.dat fluid_mesh_after_smooth.dat 
cp RESLT_solid/mesh_after_smooth.dat solid_mesh_after_smooth.dat 

cp RESLT_fluid_with_poisson/mesh_after_smooth.dat \
         fluid_with_poisson_mesh_after_smooth.dat
cp RESLT_solid_with_poisson/mesh_after_smooth.dat \
         solid_with_poisson_mesh_after_smooth.dat

cp RESLT_fluid_with_linear_elasticity/mesh_after_smooth.dat \
         fluid_with_linear_elasticity_mesh_after_smooth.dat
cp RESLT_solid_with_linear_elasticity/mesh_after_smooth.dat \
         solid_with_linear_elasticity_mesh_after_smooth.dat




echo " " >> validation.log
echo "Quadratic mesh smoothing for fluid (lin elast)" >> validation.log
echo "----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fluid_with_linear_elasticity_mesh_after_smooth.dat.gz \
      fluid_with_linear_elasticity_mesh_after_smooth.dat  >> validation.log
fi



echo " " >> validation.log
echo "Quadratic mesh smoothing for solid (lin elast)" >> validation.log
echo "----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/solid_with_linear_elasticity_mesh_after_smooth.dat.gz \
      solid_with_linear_elasticity_mesh_after_smooth.dat  >> validation.log
fi




echo " " >> validation.log
echo "Quadratic mesh smoothing for fluid (Poisson)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fluid_with_poisson_mesh_after_smooth.dat.gz \
      fluid_with_poisson_mesh_after_smooth.dat  >> validation.log
fi



echo " " >> validation.log
echo "Quadratic mesh smoothing for solid (Poisson)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/solid_with_poisson_mesh_after_smooth.dat.gz \
      solid_with_poisson_mesh_after_smooth.dat  >> validation.log
fi





echo " " >> validation.log
echo "Quadratic mesh smoothing for fluid (nonlin)" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fluid_mesh_after_smooth.dat.gz \
      fluid_mesh_after_smooth.dat  >> validation.log
fi


echo " " >> validation.log
echo "Quadratic mesh smoothing for solid (nonlin)" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/solid_mesh_after_smooth.dat.gz \
      solid_mesh_after_smooth.dat  >> validation.log
fi



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




