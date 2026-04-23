#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation
cd Validation

#######################################################################

# Validation for tensioned string (orig version; parametrisation by arclength)
#-----------------------------------------------------------------------------

mkdir RESLT

echo "Running tensioned string validation "
../tensioned_string  > OUTPUT_tensioned_string

echo "done"
echo " " >> validation.log
echo "Tensioned string validation" >> validation.log
echo "---------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/beam1.dat \
    RESLT/beam3.dat \
    RESLT/beam5.dat \
    RESLT/beam9.dat \
    RESLT/trace_beam.dat \
    > beam_results.dat

mv RESLT RESLT_orig

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/beam_results.dat.gz \
                            beam_results.dat >> validation.log
fi


# Validation for tensioned string (new version but stick with parametrisation by arclength)
#------------------------------------------------------------------------------------------

mkdir RESLT

echo "Running tensioned string validation (new but should give same results as orig)"
../tensioned_string2 --stretch_ratio 1.0 > OUTPUT_tensioned_string2

echo "done"
echo " " >> validation.log
echo "Tensioned string validation (new)" >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/beam1.dat \
    RESLT/beam3.dat \
    RESLT/beam5.dat \
    RESLT/beam9.dat \
    RESLT/trace_beam.dat \
    > beam_results_new.dat

mv RESLT RESLT_new

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
   $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/beam_results.dat.gz \
   beam_results_new.dat >> validation.log
fi



# Validation for tensioned string (new version; large sigma0; parametrisation by arclength)
#------------------------------------------------------------------------------------------

mkdir RESLT

echo "Running tensioned string validation (new; large sigma0; unit stretch ratio)"
../tensioned_string2 --sigma0 10.0 --stretch_ratio 1.0 > OUTPUT_tensioned_string_large_prestress_stretch_ratio_1pt0

echo "done"
echo " " >> validation.log
echo "Tensioned string validation (new; large sigma0; unit stretch ratio)" >> validation.log
echo "-------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/beam1.dat \
    RESLT/beam3.dat \
    RESLT/beam5.dat \
    RESLT/beam9.dat \
    RESLT/trace_beam.dat \
    > beam_results_new_large_prestress_stretch_ratio_1pt0.dat

mv RESLT RESLT_new_large_prestress_stretch_ratio_1pt0

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/beam_results_new_large_prestress_stretch_ratio_1pt0.dat.gz \
   beam_results_new_large_prestress_stretch_ratio_1pt0.dat >> validation.log
fi

# Validation for tensioned string (new version; large sigma0; parametrisation by non-arclength)
#----------------------------------------------------------------------------------------------

mkdir RESLT

echo "Running tensioned string validation (new; large sigma0; non-unit stretch ratio)"
../tensioned_string2 --sigma0 10.0 --stretch_ratio 0.9 > OUTPUT_tensioned_string_large_prestress_stretch_ratio_0.9

echo "done"
echo " " >> validation.log
echo "Tensioned string validation (new; large sigma0; non0unit stretch ratio)" >> validation.log
echo "-----------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/beam1.dat \
    RESLT/beam3.dat \
    RESLT/beam5.dat \
    RESLT/beam9.dat \
    RESLT/trace_beam.dat \
    > beam_results_new_large_prestress_stretch_ratio_0pt9.dat

mv RESLT RESLT_new_large_prestress_stretch_ratio_0pt9

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/beam_results_new_large_prestress_stretch_ratio_0pt9.dat.gz \
   beam_results_new_large_prestress_stretch_ratio_0pt9.dat >> validation.log
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
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 10
