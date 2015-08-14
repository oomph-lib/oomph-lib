#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo advection diffusion
#----------------------------------------
cd Validation

echo "Running axisymmetric insoluble surfactant Rayleigh--Plateau validation "
mkdir RESLT
../rayleigh_instability_insoluble_surfactant lalala > ./OUTPUT_ray_surf
echo "done"
echo " " >> validation.log
echo "Axisymmetric insoluble surfactant Rayleigh--Plateau  validation " >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/int5.dat  > rayleigh_surf.dat
mv RESLT RESLT_axi

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/rayleigh_surf.dat.gz \
    rayleigh_surf.dat 0.1 1.0e-14  >> validation.log
fi


echo "Running 3D insoluble surfactant Rayleigh--Plateau validation "
mkdir RESLT
../3d_rayleigh_instability_surfactant lalala > ./OUTPUT_3d_ray_surf
echo "done"
echo " " >> validation.log
echo "3D insoluble surfactant Rayleigh--Plateau  validation " >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/surface5.dat  > 3d_rayleigh_surf.dat
mv RESLT RESLT_3D

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/3d_rayleigh_surf.dat.gz \
    3d_rayleigh_surf.dat 0.1 2.0e-7  >> validation.log
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
