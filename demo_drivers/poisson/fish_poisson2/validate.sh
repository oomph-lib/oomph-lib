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

cd Validation



# Validation for fish poisson without adaptation
#-----------------------------------------------

echo "Running fish poisson without adaptation validation "
mkdir RESLT
../fish_poisson_no_adapt > OUTPUT_fish_poisson_no_adapt
echo "done"
echo " " >> validation.log
echo "Fish poisson without adaptation validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > fish_poisson_no_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_no_adapt.dat.gz   \
    fish_poisson_no_adapt.dat  >> validation.log
fi
mv RESLT RESLT_fish_poisson_no_adapt




# Validation for fish poisson with adaptation
#--------------------------------------------

echo "Running fish poisson with adaptation validation "
mkdir RESLT
../fish_poisson_adapt > OUTPUT_fish_poisson_adapt
echo "done"
echo " " >> validation.log
echo "Fish poisson with adaptation validation" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln3.dat > fish_poisson_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_adapt.dat.gz   \
    fish_poisson_adapt.dat  >> validation.log
fi
mv RESLT RESLT_fish_poisson_adapt



# Validation for fish poisson with node updates
#------------------------------------------------

echo "Running fish poisson with node updates validation "
mkdir RESLT
../fish_poisson_node_update > OUTPUT_fish_poisson_node_update
echo "done"
echo " " >> validation.log
echo "Fish poisson with node update validation" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln2.dat > fish_poisson_node_update.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fish_poisson_node_update.dat.gz   \
    fish_poisson_node_update.dat  >> validation.log
fi
mv RESLT RESLT_fish_poisson_node_update





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
