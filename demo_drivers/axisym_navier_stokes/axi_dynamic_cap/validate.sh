#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=5


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for linear stability run with restart
#--------------------------------------------------
validate(){
    # Expects to be called as
    # validate executable input actual_output_filename

    # Expected output file name is computed automatically.
  
  mkdir Validation/RESLT

  var="./$1 $2 > Validation/OUTPUT"
  echo $var
  eval $var
  echo "done"
  LOG="Validation/validation.log"
  echo " " >> $LOG 
  echo "Validation run" >> $LOG
  echo "---------------------------------------------" >> $LOG
  echo " " >> $LOG
  echo "Validation directory: " >> $LOG
  echo " " >> $LOG
  echo "  " `pwd` >> $LOG
  echo " " >> $LOG

  ## Compute expected file name 
  TEMP=$(basename -- $2)
  FILE=${1%% *}"_"${TEMP%.*}"_results.dat"
  echo Validation/$FILE
  cat  Validation/RESLT/$3 > Validation/$FILE
  rm Validation/RESLT -r
  
  if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> $LOG
  else
    ../../../bin/fpdiff.py validata/$FILE.gz  \
           Validation/$FILE 0.1 2e-7 >> $LOG
  fi
  
  # Append log to main validation log
  cat $LOG >> ../../../validation.log
}

# Utility scripts
validate "create_parameter_files --folder Validation/RESLT --overwrite --parameters" validata/parameters-with-restart.dat parameters.dat 

# Base state scripts
validate steady_run validata/parameters.dat trace.dat 
validate steady_run validata/parameters-with-restart.dat trace.dat 
validate unsteady_run validata/parameters.dat trace.dat 
validate unsteady_run validata/parameters-with-restart.dat trace.dat 

./run_tests

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
