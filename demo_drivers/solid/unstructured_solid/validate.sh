#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation


# Validation for 2D unstructured solid
#-------------------------------------


# Get triangle files
cp ../*fig.1.* .

mkdir RESLT RESLT_pres_disp RESLT_pres_disp_incomp

echo "Running 2D unstructured solid "
../unstructured_two_d_solid > OUTPUT_2D

echo "done"
echo " " >> validation.log
echo "2D unstructured solid" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_two_d.dat
cat RESLT_pres_disp/soln2.dat > 2d_pres_disp.dat
cat RESLT_pres_disp_incomp/soln2.dat > 2d_pres_disp_inc.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_two_d.dat.gz \
    result_two_d.dat  >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_pres_disp.dat.gz \
    2d_pres_disp.dat  >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_pres_disp_inc.dat.gz \
    2d_pres_disp_inc.dat  >> validation.log
fi

mv RESLT RESLT_2d
mv RESLT_pres_disp RESLT_2d_pres_disp
mv RESLT_pres_disp_incomp RESLT_2d_pres_disp_incomp

# Validation for 3D unstructured solid
#-------------------------------------

# Get tetgen files
cp ../cube_hole.* .

mkdir RESLT

echo "Running 3D unstructured solid "
../unstructured_three_d_solid > OUTPUT_3D

echo "done"
echo " " >> validation.log
echo "3D unstructured solid" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_three_d.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_three_d.dat.gz \
    result_three_d.dat  >> validation.log
fi

mv RESLT RESLT_3d


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log




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
