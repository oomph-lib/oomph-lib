#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=6


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for adaptive spin up
#--------------------------------
cd Validation

echo "Running spherical (spherical polar) spin-up validation "
mkdir RESLT_CR
mkdir RESLT_TH
../spin_up > OUTPUT_spin_up
echo "done"
echo " " >> validation.log
echo "Spherical polar spherical spin-up validation" >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_CR/soln3.dat RESLT_CR/time_trace.dat \
 > sphere_CR.dat
cat  RESLT_TH/soln3.dat RESLT_TH/time_trace.dat \
 > sphere_TH.dat

mkdir RESLT_sphere
mv RESLT_CR RESLT_TH RESLT_sphere

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Crouzeix-Raviart elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/sphere_CR.dat.gz  \
         sphere_CR.dat 0.1 5.0e-10 >> validation.log
echo "Taylor-Hood elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/sphere_TH.dat.gz  \
         sphere_TH.dat 0.1 1.0e-12 >> validation.log
fi


echo "Running spherical (cylindrical polar) spin-up validation "
mkdir RESLT_CR
mkdir RESLT_TH
../spin_up_cyl > OUTPUT_spin_up_cyl
echo "done"
echo " " >> validation.log
echo "Cylindrical polar spherical spin-up  validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_CR/soln3.dat RESLT_CR/time_trace.dat \
 > cyl_CR.dat
cat  RESLT_TH/soln3.dat RESLT_TH/time_trace.dat \
 > cyl_TH.dat

mkdir RESLT_cyl
mv RESLT_CR RESLT_TH RESLT_cyl

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Crouzeix-Raviart elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/cyl_CR.dat.gz  \
         cyl_CR.dat 0.1 1.0e-12 >> validation.log
echo "Taylor-Hood elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/cyl_TH.dat.gz  \
         cyl_TH.dat 0.1 1.0e-12 >> validation.log
fi



echo "Running refineable spherical (spherical polar) spin-up validation "
mkdir RESLT
../refineable_spin_up > OUTPUT_ref_spin_up
echo "done"
echo " " >> validation.log
echo "Refineable Spherical polar spherical spin-up validation" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln1.dat > ref_sphere_CR.dat

mkdir RESLT_ref_sphere
mv RESLT RESLT_ref_sphere

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Crouzeix-Raviart elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/ref_sphere_CR.dat.gz  \
         ref_sphere_CR.dat 0.1 1.0e-12 >> validation.log
#echo "Taylor-Hood elements" >> validation.log
#echo " " >> validation.log
#../../../../bin/fpdiff.py ../validata/sphere_TH.dat.gz  \
#         sphere_TH.dat 0.1 1.0e-12 >> validation.log
fi



echo "Running refineable spherical (cylindrical polar) spin-up validation "
mkdir RESLT
../refineable_spin_up_cyl > OUTPUT_ref_spin_up_cyl
echo "done"
echo " " >> validation.log
echo "Refineable Cylindrical polar spherical spin-up validation" >> validation.log
echo "----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln1.dat > ref_cyl_CR.dat

mkdir RESLT_ref_cyl
mv RESLT RESLT_ref_cyl

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Crouzeix-Raviart elements" >> validation.log
echo " " >> validation.log
../../../../bin/fpdiff.py ../validata/ref_cyl_CR.dat.gz  \
         ref_cyl_CR.dat 0.1 1.0e-12 >> validation.log
#echo "Taylor-Hood elements" >> validation.log
#echo " " >> validation.log
#../../../../bin/fpdiff.py ../validata/sphere_TH.dat.gz  \
#         sphere_TH.dat 0.1 1.0e-12 >> validation.log
fi



# Append log to main validation log
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
