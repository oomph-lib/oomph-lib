#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=13

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for BlockSelector object
#-----------------------------------------
cd Validation

echo "Running BlockSelector validation "

TESTPROGRAM="block_selector_test"
## First untar the validata
cd ../validata && tar -xzf OUTPUT.tar.gz
cd ../Validation

###########################################################################
../$TESTPROGRAM --test_default_constructor
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: Default constructor test" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_default_constructor_test  \
        OUTFILE_default_constructor_test >> validation.log
fi

###########################################################################
../$TESTPROGRAM --test_constructor_with_param
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: Constructor with parameters" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_constructor_with_param_wanted  \
        OUTFILE_constructor_with_param_wanted >> validation.log

    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_constructor_with_param_not_wanted  \
        OUTFILE_constructor_with_param_not_wanted >> validation.log

    ./.././validate_helper.sh ./../validata/OUTFILE_constructor_with_param_replace \
        OUTFILE_constructor_with_param_replace >> validation.log 
fi


###########################################################################
../$TESTPROGRAM --test_select_block_function
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: select_block function" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_select_block_wanted  \
        OUTFILE_select_block_wanted >> validation.log

    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_select_block_not_wanted  \
        OUTFILE_select_block_not_wanted >> validation.log

    ./.././validate_helper.sh ./../validata/OUTFILE_select_block_replace \
        OUTFILE_select_block_replace >> validation.log 
fi

###########################################################################
../$TESTPROGRAM --test_want_block_function
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: want_block function" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_want_block  \
        OUTFILE_want_block >> validation.log
fi
###########################################################################
../$TESTPROGRAM --test_do_not_want_block_function
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: do_not_want_block function" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_do_not_want_block  \
        OUTFILE_do_not_want_block >> validation.log
fi

###########################################################################
../$TESTPROGRAM --test_do_not_want_block_function_replace > do_not_want_block_warning 2>&1
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: do_not_want_block function with replacement" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_do_not_want_block_replace  \
        OUTFILE_do_not_want_block_replace >> validation.log
fi

## Now check for the warning.
REPLACEWARNING=$(grep "Oomph-lib WARNING" do_not_want_block_warning)
WARNSTR="WARNING"

if test "${REPLACEWARNING#*$WARNSTR}" != "$REPLACEWARNING"
then
  echo " [OK] found the warning." >> validation.log
else
  echo " [FAILED] cannot find warning." >> validation.log
fi

###########################################################################
../$TESTPROGRAM --test_set_row_index_function
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: set_row_index function" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_set_row_index  \
        OUTFILE_set_row_index >> validation.log
fi

###########################################################################
../$TESTPROGRAM --test_set_column_index_function
echo "done"
echo " " >> validation.log
echo "BlockSelector validation: set_column_index function" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/OUTFILE_set_column_index  \
        OUTFILE_set_column_index >> validation.log
fi


##########################################################################
## All test done, remove the extracted files in validata
cd ../validata
rm -rf OUTFILE_*
cd ../Validation

###########################################################################

# Append log to main validation log
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
