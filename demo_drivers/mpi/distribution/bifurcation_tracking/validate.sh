#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for distributed pitchfork bifurcation detection
#---------------------------------------------------------------------

mkdir RESLT RESLT_track_pitch
cd RESLT

echo "Running pitchfork bifurcation detection (parallel)  "
$MPI_RUN_COMMAND ../../pitchfork > ../OUTPUT_pitchfork

echo "done"
cd ..
echo " " >> validation.log
echo "Pitchfork bifurcation detection (parallel)  validation" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace0.dat RESLT/trace1.dat > pitch.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/pitch.dat.gz \
 pitch.dat 0.1 1.0e-13 >> validation.log
fi


mv RESLT RESLT_pitchfork


echo "Running pitchfork tracking validation (parallel)"
cd RESLT_track_pitch
$MPI_RUN_COMMAND ../../track_pitch > ../OUTPUT_track_pitch
cd ..

echo "done"
echo " " >> validation.log
echo "Pitchfork bifurcation tracking validation (parallel)" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_track_pitch/trace_pitch0.dat > track_pitch.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/track_pitch.dat.gz \
    track_pitch.dat  0.1 5.0e-8  >> validation.log
fi


# Append output to global validation log file
#--------------------------------------------
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
