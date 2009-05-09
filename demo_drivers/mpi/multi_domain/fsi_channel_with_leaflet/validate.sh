#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=1

# Doc what we're using to run tests on two processors
echo " " 
echo "Running mpi tests with mpi run command: " $MPI_RUN_COMMAND
echo " " 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation
cp ../*partition.dat .

# Validation for FSI channel with leaflet problem (algebraic node update)
#------------------------------------------------------------------------

echo "Running FSI channel with leaflet (algebraic node update) validation "
mkdir RESLT_FSI_LEAF_EXT
$MPI_RUN_COMMAND ../fsi_channel_with_leaflet_external validate > OUTPUT_fsi_channel_with_leaflet_external
echo "done"
echo " " >> validation.log
echo "FSI channel with leaflet (algebraic) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_FSI_LEAF_EXT/soln1_on_proc0.dat RESLT_FSI_LEAF_EXT/soln1_on_proc1.dat RESLT_FSI_LEAF_EXT/soln3_on_proc0.dat RESLT_FSI_LEAF_EXT/soln3_on_proc1.dat RESLT_FSI_LEAF_EXT/wall_soln2_on_proc0.dat RESLT_FSI_LEAF_EXT/wall_soln2_on_proc1.dat > fsi_channel_with_leaflet_external_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/fsi_channel_with_leaflet_external_results.dat.gz  \
         fsi_channel_with_leaflet_external_results.dat >> validation.log
fi

#---------------------------------------------------------------------

# Append log to main validation log
cat validation.log >> ../../../../../validation.log

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
