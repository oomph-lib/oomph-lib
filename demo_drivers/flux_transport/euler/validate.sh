#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running 1D Euler validation "
mkdir RESLT_1D
cd RESLT_1D
../../one_d_euler > ../OUTPUT_1d_euler
cd ..
echo "done"
echo " " >> validation.log
echo "1D Euler validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_1D/sod_100_time0.2.dat > 1d_sod.dat
cat RESLT_1D/lax_100_time0.2.dat > 1d_lax.dat

echo "Sod problem " >> validation.log
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/1d_sod.dat.gz \
    1d_sod.dat  0.1 3.0e-12 >> validation.log
fi

echo "Lax problem " >> validation.log
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/1d_lax.dat.gz \
    1d_lax.dat  0.1 3.0e-12 >> validation.log
fi

echo "Running 2D Euler validation "
mkdir RESLT_2D
cd RESLT_2D
../../two_d_euler > ../OUTPUT_2d_euler
cd ..
echo "done"
echo " " >> validation.log
echo "2D Euler validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_2D/trace_disc.dat RESLT_2D/disc_256_time0.05.dat > 2d_euler.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_euler.dat.gz \
    2d_euler.dat  0.1 1.0e-12 >> validation.log
fi


echo "Running 2D Euler (Couette) validation "
mkdir RESLT_couette
cd RESLT_couette
../../couette > ../OUTPUT_2d_couette
cd ..
echo "done"
echo " " >> validation.log
echo "2D Euler (Couette) validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_couette/trace_disc.dat RESLT_couette/disc_64_time0.2.dat > couette.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/couette.dat.gz \
    couette.dat  0.1 1.0e-12 >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log

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
