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

echo "Running Axisymmetric convection from sphere"
mkdir RESLT
cd RESLT
../../b_convection_sphere > ../OUTPUT_b_convection_sphere
cd ..
echo "done"
echo " " >> validation.log
echo "Axisymmetric convection from sphere" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln_Re1_Ra50.dat RESLT/trace.dat > b_conv_sph.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/b_conv_sph.dat.gz   \
    b_conv_sph.dat 0.1 1.0e-14  >> validation.log
fi


echo "Running Axisymmetric convection from sphere (multi domain)"
mkdir RESLT_multi
cd RESLT_multi
../../multi_domain_b_convect_sph > ../OUTPUT_multi_domain_b_convect_sph
cd ..
echo "done"
echo " " >> validation.log
echo "Axisymmetric convection from sphere (multi domain)" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_multi/vel_soln_Re1_Ra50.dat RESLT_multi/conc_soln_Re1_Ra50.dat > b_conv_sph_multi.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/b_conv_sph_multi.dat.gz   \
    b_conv_sph_multi.dat 0.1 1.0e-14  >> validation.log
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
