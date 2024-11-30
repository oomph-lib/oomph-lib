#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=7

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for fold tracking 
#---------------------------------------

cd Validation
mkdir RESLT_fold RESLT_pitch RESLT_hopf RESLT_adaptive_pitch \
      RESLT_adaptive_hopf RESLT_track_pitch RESLT_periodic_orbit

echo "Running fold bifurcation validation "
cd RESLT_fold
../../fold > ../OUTPUT_fold
cd ..

echo "done"
echo " " >> validation.log
echo "Fold bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_fold/trace_mu.dat > fold.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/fold.dat.gz \
    fold.dat  >> validation.log
fi


echo "Running pitchfork bifurcation validation "
cd RESLT_pitch
../../pitchfork > ../OUTPUT_pitchfork
cd ..

echo "done"
echo " " >> validation.log
echo "Pitchfork bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_pitch/trace.dat > pitch.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/pitch.dat.gz pitch.dat >> validation.log
fi

echo "Running hopf bifurcation validation "
cd RESLT_hopf
../../hopf > ../OUTPUT_hopf
cd ..

echo "done"
echo " " >> validation.log
echo "Hopf bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_hopf/trace_hopf.dat > hopf.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/hopf.dat.gz \
    hopf.dat  0.6 1.0e-14  >> validation.log
fi

echo "Running pitchfork tracking validation "
cd RESLT_track_pitch
../../track_pitch > ../OUTPUT_track_pitch
cd ..

echo "done"
echo " " >> validation.log
echo "Pitchfork bifurcation tracking validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_track_pitch/trace_pitch.dat > track_pitch.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/track_pitch.dat.gz \
    track_pitch.dat  0.1 3.0e-10  >> validation.log
fi

echo "Running adaptive pitchfork bifurcation validation "
cd RESLT_adaptive_pitch
../../adaptive_pitchfork > ../OUTPUT_adaptive_pitchfork
cd ..

echo "done"
echo " " >> validation.log
echo "Adaptive pitchfork bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_adaptive_pitch/trace.dat RESLT_adaptive_pitch/bif_soln.dat \
    > adaptive_pitch.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/adaptive_pitch.dat.gz adaptive_pitch.dat \
 0.1 5.0e-7 >> validation.log
fi



echo "Running adaptive hopf bifurcation validation "
cd RESLT_adaptive_hopf
ln -sf ../../adapt_hopf_eigen.dat eigen.dat
../../adaptive_hopf > ../OUTPUT_adaptive_hopf
cd ..

echo "done"
echo " " >> validation.log
echo "Adaptive hopf bifurcation validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_adaptive_hopf/trace.dat RESLT_adaptive_hopf/bif_soln.dat \
    > adaptive_hopf.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/adaptive_hopf.dat.gz adaptive_hopf.dat \
 0.1 5.0e-8 >> validation.log
fi

echo "Running periodic_orbit validation "
cd RESLT_periodic_orbit
../../periodic_orbit > ../OUTPUT_periodic_orbit
cd ..

echo "done"
echo " " >> validation.log
echo "Periodic orbit validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_periodic_orbit/trace.dat RESLT_periodic_orbit/orbit_trace.dat \
  RESLT_periodic_orbit/first_orbit.dat > periodic_orbit.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/periodic_orbit.dat.gz periodic_orbit.dat \
 0.1 1.0e-12 >> validation.log
fi


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../validation.log


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
