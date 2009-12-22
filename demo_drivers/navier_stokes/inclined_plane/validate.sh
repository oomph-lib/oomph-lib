#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for single layer free surface Navier Stokes problem
#---------------------------------------------------------------
cd Validation

echo "Running inclined plane free surface Navier--Stokes validation "
mkdir RESLT_TH RESLT_CR
cd RESLT_TH
../../inclined_plane > ../OUTPUT_inclined_plane
echo "Taylor Hood done"
cd ..
cd RESLT_CR
../../inclined_plane2 > ../OUTPUT_inclined_plane_CR
echo "Crouzeix Raviart done"
cd ..
echo " " >> validation.log
echo "Inclined plane free surface Navier--Stokes validation" >> validation.log
echo "------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_TH/spine_output.dat    RESLT_TH/spine_step*_2.dat  > TH_spine.dat
cat RESLT_TH/elastic_output.dat  RESLT_TH/elastic_step*_2.dat  > TH_elastic.dat
cat RESLT_CR/spine_output.dat    RESLT_CR/spine_step*_2.dat  > CR_spine.dat
cat RESLT_CR/elastic_output.dat  RESLT_CR/elastic_step*_2.dat  > CR_elastic.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Taylor Hood Elements: Spine " >> validation.log
../../../../bin/fpdiff.py ../validata/TH_spine.dat.gz  \
         TH_spine.dat 0.1 1.0e-12 >> validation.log
echo "Taylor Hood Elements: Pseudo-Elastic " >> validation.log
../../../../bin/fpdiff.py ../validata/TH_elastic.dat.gz  \
         TH_elastic.dat 0.1 1.0e-12 >> validation.log

echo "Crouzeix Raviart Elements: Spine " >> validation.log
../../../../bin/fpdiff.py ../validata/CR_spine.dat.gz  \
         CR_spine.dat 0.1 1.0e-12 >> validation.log

echo "Crouzeix Raviart Elements: Pseudo-Elastic " >> validation.log
../../../../bin/fpdiff.py ../validata/CR_elastic.dat.gz  \
         CR_elastic.dat 0.1 1.0e-12 >> validation.log
fi


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
