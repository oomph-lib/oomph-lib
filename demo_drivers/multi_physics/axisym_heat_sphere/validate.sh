#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation 
#----------------------------
cd Validation

echo "Running Axisymmetric heat transfer from sphere"
mkdir RESLT
cd RESLT
../../axisym_heat_sphere > ../OUTPUT_axisym_heat_sphere
cd ..
echo "done"
echo " " >> validation.log
echo "Axisymmetric heat transfer from sphere" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln_Re10.dat RESLT/trace.dat > heat_sphere.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/heat_sphere.dat.gz   \
    heat_sphere.dat 0.1 1.0e-14  >> validation.log
fi


echo "Running Axisymmetric heat transfer from sphere (multi domain)"
mkdir RESLT_multi
cd RESLT_multi
../../multi_domain_axisym_heat_sphere > ../OUTPUT_multi_domain_axisym_heat_sphere
cd ..
echo "done"
echo " " >> validation.log
echo "Axisymmetric heat transfer from sphere" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_multi/vel_soln_Re10.dat RESLT_multi/conc_soln_Re10.dat > heat_sphere_multi.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/heat_sphere_multi.dat.gz   \
    heat_sphere_multi.dat 0.1 1.0e-14  >> validation.log
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
