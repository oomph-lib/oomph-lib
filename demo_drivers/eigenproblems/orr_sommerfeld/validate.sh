#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for orr sommerfeld
#------------------------------
cd Validation

if [ -f ../../orr_sommerfeld ]; then


echo "Running orr_sommerfeld validation "
mkdir RESLT
cd RESLT
../../orr_sommerfeld > ../OUTPUT_orr_sommerfeld
cd ..
echo "done"
echo " " >> validation.log
echo "Orr-Sommerfeld validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/neutral.dat > orr_sommerfeld_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/orr_sommerfeld_results.dat.gz   \
    orr_sommerfeld_results.dat  0.3  1.0e-14 >> validation.log
fi


else

echo ""
echo "Not running orr_sommerfeld test; needs trilinos"
echo ""
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log

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
