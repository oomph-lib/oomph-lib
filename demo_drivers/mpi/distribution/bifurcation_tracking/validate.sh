#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

# Receive the mpirun command as the second argument
MPI_RUN_COMMAND="$2"

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

if test "$3" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/pitch.dat.gz \
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

if test "$3" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/track_pitch.dat.gz \
    track_pitch.dat  0.1 1e-06  >> validation.log
fi


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> $OOMPH_ROOT_DIR/validation.log

cd ..


#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 10
