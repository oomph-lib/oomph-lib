#! /bin/sh

#------------------------------------------------------------------
# Script to doc the number of successes/failures of all self tests
#------------------------------------------------------------------



#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' ../../validation.log`
FAIL_COUNT=`grep -c 'FAILED' ../../validation.log`

cp ../../validation.log .

echo " "
echo "======================================================================"
echo " " 
if  [ $FAIL_COUNT -eq 0 ]; then
  echo "All $OK_COUNT test[s] passed successfully. "
  echo " "
  echo "See " 
  echo " "
  echo "    `pwd`/validation.log    "
  echo " "
  echo "for details."
  echo " "
  echo "======================================================================"
  echo " " 
  exit 0
else 
  echo "Only $OK_COUNT test[s] passed successfully "
  echo "while $FAIL_COUNT test[s] failed! "
  echo " "
  echo "See " 
  echo " "
  echo "    `pwd`/validation.log    "
  echo " "
  echo "for details."
  echo " "
  ../../bin/find_failed_directories.bash
  echo " "
  echo "======================================================================"
  echo " " 
  exit 1
fi


