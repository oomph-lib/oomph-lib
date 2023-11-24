#! /bin/sh

#------------------------------------------------------------------
# Script to doc the number of successes/failures of all self tests
#------------------------------------------------------------------



#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' ../../validation.log`
FAIL_COUNT=`grep -c 'FAIL' ../../validation.log`

cp ../../validation.log .

echo " "
echo "======================================================================"
echo " " 
if  [ $FAIL_COUNT -eq 0 ]; then
  echo "All $OK_COUNT compiled test[s] passed successfully. "
  echo " "
  echo "See " 
  echo " "
  echo "    `pwd`/validation.log    "
  echo " "
  echo "for details."
  echo " "
  return_flag=0
  cd ../..; bin/find_compilation_failures_for_demo_drivers.bash | tee .failed_compilation.txt
  failures=`grep "Check these out" .failed_compilation.txt | wc -w`
  rm -f .failed_compilation.txt
  if [ $failures -ne 0 ]; then
      return_flag=1
  fi
  exit $return_flag
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
  cd ../..; bin/find_compilation_failures_for_demo_drivers.bash
  exit 1
fi


