#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2



# Compile code manually outside automake framework
#-------------------------------------------------
make -f makefile.sample my_demo_code



# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation



# Validation with automake/autoconf
#----------------------------------
cd Validation

echo "Running linking test: with automake/autoconf"
mkdir RESLT
cd RESLT
../../demo_code > ../OUTPUT_with_auto
cd ..
echo "done"
echo " " >> validation.log
echo "linking test with automake/autoconf" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat > results_with.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results_with.dat  >> validation.log
fi

mv RESLT RESLT_with_auto



# Validation without automake/autoconf
#----------------------------------

echo "Running linking test: without automake/autoconf"
mkdir RESLT
cd RESLT
../../my_demo_code > ../OUTPUT_without_auto
cd ..
echo "done"
echo " " >> validation.log
echo "linking test without automake/autoconf" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat > results_without.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results_without.dat  >> validation.log
fi

mv RESLT RESLT_without_auto





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
