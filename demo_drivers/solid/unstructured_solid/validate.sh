#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation


# Validation for 2D unstructured solid
#-------------------------------------


# Get triangle files
cp ../*fig.1.* .

mkdir RESLT RESLT_pres_disp RESLT_pres_disp_incomp

echo "Running 2D unstructured solid "
../unstructured_two_d_solid > OUTPUT_2D

echo "done"
echo " " >> validation.log
echo "2D unstructured solid" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_two_d.dat
cat RESLT_pres_disp/soln2.dat > 2d_pres_disp.dat
cat RESLT_pres_disp_incomp/soln2.dat > 2d_pres_disp_inc.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_two_d.dat.gz \
    result_two_d.dat  >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_pres_disp.dat.gz \
    2d_pres_disp.dat  >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_pres_disp_inc.dat.gz \
    2d_pres_disp_inc.dat  >> validation.log
fi

mv RESLT RESLT_2d
mv RESLT_pres_disp RESLT_2d_pres_disp
mv RESLT_pres_disp_incomp RESLT_2d_pres_disp_incomp

# Validation for 3D unstructured solid
#-------------------------------------

# Get tetgen files
cp ../cube_hole.* .

mkdir RESLT

echo "Running 3D unstructured solid "
../unstructured_three_d_solid > OUTPUT_3D

echo "done"
echo " " >> validation.log
echo "3D unstructured solid" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_three_d.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_three_d.dat.gz \
    result_three_d.dat  >> validation.log
fi

mv RESLT RESLT_3d


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log




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

