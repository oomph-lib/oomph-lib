#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validation for mesh generation from triangle with box mesh
#------------------------------------------------------------

echo "Running triangle mesh generation Poisson validation with box mesh"
mkdir RESLT
../mesh_from_triangle_poisson ../box_hole.1.node ../box_hole.1.ele \
../box_hole.1.poly > OUTPUT_box_poisson  

echo "done"
echo " " >> validation.log
echo "triangle mesh generation Poisson validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat  \
    > box_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/box_poisson_results.dat.gz   \
    box_poisson_results.dat  >> validation.log
fi

mv RESLT RESLT_box_poisson




# Validation for mesh generation from triangle with Navier Stokes
#-----------------------------------------------------------------

echo "Running triangle mesh generation Navier Stokes validation with box mesh"
mkdir RESLT
../mesh_from_triangle_navier_stokes ../flow_past_box.1.node \
../flow_past_box.1.ele ../flow_past_box.1.poly \
 > OUTPUT_box_navier_stokes  

echo "done"
echo " " >> validation.log
echo "triangle mesh generation Navier Stokes validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat \
    > box_navier_stokes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/box_navier_stokes_results.dat.gz   \
    box_navier_stokes_results.dat  >> validation.log
fi

mv RESLT RESLT_box_navier_stokes



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
