#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for vmtk fsi
#------------------------
cd Validation

echo "Running quadratic vmtk geometry test"

# Get the mesh geometry
ln -s ../fluid_iliac_short_fine.1.ele  fluid.1.ele  
ln -s ../fluid_iliac_short_fine.1.face fluid.1.face
ln -s ../fluid_iliac_short_fine.1.node fluid.1.node
ln -s ../solid_iliac_short_fine.1.ele  solid.1.ele 
ln -s ../solid_iliac_short_fine.1.face solid.1.face
ln -s ../solid_iliac_short_fine.1.node solid.1.node
ln -s ../short_fine_iliac_boundary_enumeration.dat boundary_enumeration.dat
ln -s ../short_fine_iliac_quadratic_fsi_boundary.dat quadratic_fsi_boundary.dat
ln -s ../short_fine_iliac_quadratic_outer_solid_boundary.dat quadratic_outer_solid_boundary.dat


mkdir RESLT

../unsteady_vmtk_fsi quadratic_test > OUTPUT_quadratic_test
echo "done"
echo " " >> validation.log
echo "vmtk FSI validation" >> validation.log
echo "-------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_nodes0.dat   \
    RESLT/solid_nodes0.dat   \
    > quadratic_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/quadratic_results.dat.gz  \
         quadratic_results.dat 0.1 1.0e-11 >> validation.log
fi

mv RESLT RESLT_quadratic_test
rm  fluid.1.ele  
rm  fluid.1.face
rm  fluid.1.node
rm  solid.1.ele 
rm  solid.1.face
rm  solid.1.node
rm  boundary_enumeration.dat
rm  quadratic_fsi_boundary.dat
rm  quadratic_outer_solid_boundary.dat



# Get the mesh geometry
ln -s ../fluid_iliac_short_coarse.1.ele  fluid.1.ele  
ln -s ../fluid_iliac_short_coarse.1.face fluid.1.face
ln -s ../fluid_iliac_short_coarse.1.node fluid.1.node
ln -s ../solid_iliac_short_coarse.1.ele  solid.1.ele 
ln -s ../solid_iliac_short_coarse.1.face solid.1.face
ln -s ../solid_iliac_short_coarse.1.node solid.1.node
ln -s ../short_coarse_iliac_boundary_enumeration.dat boundary_enumeration.dat
ln -s ../short_coarse_iliac_quadratic_fsi_boundary.dat quadratic_fsi_boundary.dat
ln -s ../short_coarse_iliac_quadratic_outer_solid_boundary.dat quadratic_outer_solid_boundary.dat


mkdir RESLT

../unsteady_vmtk_fsi coarse_test > OUTPUT_coarse_test
echo "done"
echo " " >> validation.log
echo "vmtk FSI validation" >> validation.log
echo "-------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/fluid_soln0.dat  RESLT/fluid_soln3.dat   \
    RESLT/solid_soln0.dat  RESLT/solid_soln3.dat   \
    > coarse_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/coarse_results.dat.gz  \
         coarse_results.dat 0.1 1.0e-11 >> validation.log
fi

mv RESLT RESLT_coarse_test
rm -f fluid.1.ele  
rm -f fluid.1.face
rm -f fluid.1.node
rm -f solid.1.ele 
rm -f solid.1.face
rm -f solid.1.node
rm -f boundary_enumeration.dat
rm -f quadratic_fsi_boundary.dat
rm -f quadratic_outer_solid_boundary.dat

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
