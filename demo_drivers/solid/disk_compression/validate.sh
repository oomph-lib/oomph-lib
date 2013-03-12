#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=6

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for static disk compression
#---------------------------------------

cd Validation
mkdir RESLT_static_disk_compression0
ln -s RESLT_static_disk_compression0 RESLT0

mkdir RESLT_static_disk_compression1
ln -s RESLT_static_disk_compression1 RESLT1

mkdir RESLT_static_disk_compression2
ln -s RESLT_static_disk_compression2 RESLT2

mkdir RESLT_static_disk_compression3
ln -s RESLT_static_disk_compression3 RESLT3

mkdir RESLT_static_disk_compression4
ln -s RESLT_static_disk_compression4 RESLT4

mkdir RESLT_static_disk_compression5
ln -s RESLT_static_disk_compression5 RESLT5


echo "Running static disk compression validation "
../disk_compression la > OUTPUT_static_disk_compression


echo "done"
echo " " >> validation.log
echo "Static disk compression validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_static_disk_compression0/soln2.dat \
    RESLT_static_disk_compression0/trace.dat \
    > static_disk_compression_results0.dat

cat RESLT_static_disk_compression1/soln2.dat \
    RESLT_static_disk_compression1/trace.dat \
    > static_disk_compression_results1.dat

cat RESLT_static_disk_compression2/soln2.dat \
    RESLT_static_disk_compression2/trace.dat \
    > static_disk_compression_results2.dat

cat RESLT_static_disk_compression3/soln2.dat \
    RESLT_static_disk_compression3/trace.dat \
    > static_disk_compression_results3.dat

cat RESLT_static_disk_compression4/soln2.dat \
    RESLT_static_disk_compression4/trace.dat \
    > static_disk_compression_results4.dat

cat RESLT_static_disk_compression5/soln2.dat \
    RESLT_static_disk_compression5/trace.dat \
    > static_disk_compression_results5.dat



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/static_disk_compression_results0.dat.gz \
    static_disk_compression_results0.dat  >> validation.log

../../../../bin/fpdiff.py ../validata/static_disk_compression_results1.dat.gz \
    static_disk_compression_results1.dat  >> validation.log

../../../../bin/fpdiff.py ../validata/static_disk_compression_results2.dat.gz \
    static_disk_compression_results2.dat  >> validation.log

../../../../bin/fpdiff.py ../validata/static_disk_compression_results3.dat.gz \
    static_disk_compression_results3.dat  >> validation.log

../../../../bin/fpdiff.py ../validata/static_disk_compression_results4.dat.gz \
    static_disk_compression_results4.dat  >> validation.log

../../../../bin/fpdiff.py ../validata/static_disk_compression_results5.dat.gz \
    static_disk_compression_results5.dat  >> validation.log
fi


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' Validation/validation.log`
if  [ $OK_COUNT -eq $NUM_TESTS ]; then
 echo " "
 echo "======================================================================"
 echo " " 
 echo "All tests in" 
 echo " " 
 echo "    `pwd`    "
 echo " "
 echo "passed successfully."
 echo " "
 echo "======================================================================"
 echo " " 
 exit 0
else
  if [ $OK_COUNT -lt $NUM_TESTS ]; then
   echo " "
   echo "======================================================================"
   echo " " 
   echo "Only $OK_COUNT of $NUM_TESTS test(s) passed; see"
   echo " " 
   echo "    `pwd`/Validation/validation.log"
   echo " " 
   echo "for details" 
   echo " " 
   echo "======================================================================"
   echo " "
   exit 1
  else 
   echo " "
   echo "======================================================================"
   echo " " 
   echo "More OKs than tests! Need to update NUM_TESTS in"
   echo " " 
   echo "    `pwd`/validate.sh"
   echo " "
   echo "======================================================================"
   echo " "
  exit 2
  fi
fi

# Never get here
exit 10
