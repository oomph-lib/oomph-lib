#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running 1D discontinous/continuous advection validation "
mkdir RESLT_1D
cd RESLT_1D
../../one_d_advection > ../OUTPUT_1d_adv
cd ..
echo "done"
echo " " >> validation.log
echo "1D continuous convection validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_1D/trace_cont.dat RESLT_1D/cont_8_time1.dat > 1d_cont_adv.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/1d_cont_adv.dat.gz \
    1d_cont_adv.dat  0.1 1.0e-14 >> validation.log
fi

echo " " >> validation.log
echo "1D discontinuous convection validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_1D/trace_disc.dat RESLT_1D/disc_8_time1.dat > 1d_disc_adv.dat



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/1d_disc_adv.dat.gz \
    1d_disc_adv.dat 0.1 1.0e-14 >> validation.log
fi


echo "Running 2D discontinous advection validation "
mkdir RESLT_2D
cd RESLT_2D
../../two_d_advection > ../OUTPUT_2d_adv
cd ..
echo "done"
echo " " >> validation.log
echo "2D discontinuous convection validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_2D/trace_disc.dat RESLT_2D/disc_64_time0.5.dat > 2d_disc_adv.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_disc_adv.dat.gz \
    2d_disc_adv.dat  0.1 1.0e-14 >> validation.log
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
