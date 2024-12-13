#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

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
ln -sf RESLT_static_disk_compression0 RESLT0

mkdir RESLT_static_disk_compression1
ln -sf RESLT_static_disk_compression1 RESLT1

mkdir RESLT_static_disk_compression2
ln -sf RESLT_static_disk_compression2 RESLT2

mkdir RESLT_static_disk_compression3
ln -sf RESLT_static_disk_compression3 RESLT3

mkdir RESLT_static_disk_compression4
ln -sf RESLT_static_disk_compression4 RESLT4

mkdir RESLT_static_disk_compression5
ln -sf RESLT_static_disk_compression5 RESLT5


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
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
