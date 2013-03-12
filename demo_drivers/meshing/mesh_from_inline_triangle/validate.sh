#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validate mesh generation from inline triangle (curvilinear outer boundary)
#---------------------------------------------------------------------------

echo "Running inline triangle mesh generation test (curvilinear outer boundary)"
mkdir RESLT
../mesh_from_inline_triangle --validation > OUTPUT

echo "done"
echo " " >> validation.log
echo "triangle inline mesh generation test (curvilinear outer boundary)" >> validation.log
echo "-----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  \
    > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
fi

mv RESLT RESLT_curvilinear_outer_boundary



# Validate mesh generation from inline triangle (polygonal outer boundary)
#---------------------------------------------------------------------------

echo "Running inline triangle mesh generation test (polygonal outer boundary)"
mkdir RESLT
../mesh_from_inline_triangle_polygon --validation > OUTPUT

echo "done"
echo " " >> validation.log
echo "triangle inline mesh generation test (polygonal outer boundary)" >> validation.log
echo "--------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  \
    > results_polygon.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_polygon.dat.gz   \
    results_polygon.dat  >> validation.log
fi

mv RESLT RESLT_polygonal_outer_boundary




# Validate mesh generation from inline triangle (polygonal outer boundary)
#---------------------------------------------------------------------------

echo "Running inline triangle mesh generation test (no adaptation)"
mkdir RESLT
../mesh_from_inline_triangle_no_adapt --validation > OUTPUT

echo "done"
echo " " >> validation.log
echo "triangle inline mesh generation test (no adapt)" >> validation.log
echo "-----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  \
    > results_no_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_no_adapt.dat.gz   \
    results_no_adapt.dat  >> validation.log
fi

mv RESLT RESLT_no_adapt


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
