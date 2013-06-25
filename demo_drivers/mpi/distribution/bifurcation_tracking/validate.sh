#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

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
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
