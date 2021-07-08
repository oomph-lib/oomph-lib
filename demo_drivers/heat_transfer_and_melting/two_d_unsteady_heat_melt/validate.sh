#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=14


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation
cd Validation


# Validation for linear solid contact with gravity (unstructured)
#----------------------------------------------------------------

echo "Running solid linear contact with gravity (unstructured) validation"
mkdir RESLT
../linear_solid_contact_with_gravity_unstructured  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "Linear solid contact with gravity (unstructured) validation" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_linear_solid_contact_with_gravity_unstructured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_linear_solid_contact_with_gravity_unstructured.dat.gz \
    result_linear_solid_contact_with_gravity_unstructured.dat 2.5 3.1e-6 >> validation.log
fi

mv RESLT RESLT_linear_solid_contact_with_gravity_unstructured


# Validation for linear solid contact with gravity (structured)
#--------------------------------------------------------------

echo "Running linear solid contact with gravity (structured) validation"
mkdir RESLT
../linear_solid_contact_with_gravity_structured  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "Linear solid contact with gravity (structured) validation" >> validation.log
echo "---------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_linear_solid_contact_with_gravity_structured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_linear_solid_contact_with_gravity_structured.dat.gz \
    result_linear_solid_contact_with_gravity_structured.dat  >> validation.log
fi

mv RESLT RESLT_linear_solid_contact_with_gravity_structured


# Validation for solid contact with gravity (unstructured)
#---------------------------------------------------------

echo "Running solid contact with gravity (unstructured) validation"
mkdir RESLT
../solid_contact_with_gravity_unstructured  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "Solid contact with gravity (unstructured) validation" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_solid_contact_with_gravity_unstructured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_solid_contact_with_gravity_unstructured.dat.gz \
    result_solid_contact_with_gravity_unstructured.dat 1.1 2.0e-7 >> validation.log
fi

mv RESLT RESLT_solid_contact_with_gravity_unstructured


# Validation for solid contact with gravity (structured)
#-------------------------------------------------------

echo "Running solid contact with gravity (structured) validation"
mkdir RESLT
../solid_contact_with_gravity_structured  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "Solid contact with gravity (structured) validation" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_solid_contact_with_gravity_structured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_solid_contact_with_gravity_structured.dat.gz \
    result_solid_contact_with_gravity_structured.dat  >> validation.log
fi

mv RESLT RESLT_solid_contact_with_gravity_structured


# Validation for solid contact
#-----------------------------

echo "Running solid contact validation"
mkdir RESLT
../solid_contact --validate  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "Solid contact validation" >> validation.log
echo "------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_solid_contact.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_solid_contact.dat.gz \
    result_solid_contact.dat  >> validation.log
fi

mv RESLT RESLT_solid_contact


# Validation for Stefan Boltzmann melting
#----------------------------------------

echo "Running Stefan Boltzmann melt validation"
mkdir RESLT
../stefan_boltzmann_melt --dt 0.05 --t_max 0.2  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D Stefan Boltzmann melt validation" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_sb_melt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_sb_melt.dat.gz \
    result_sb_melt.dat 1.0e-14 0.5  >> validation.log
fi

mv RESLT RESLT_sb_melt





# Validation 2D unsteady heat
#----------------------------
echo "Running 2D unsteady heat validation "
mkdir RESLT
../two_d_unsteady_heat  > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady heat validation " >> validation.log
echo "----------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat \
    > result_unsteady_heat.dat
mv RESLT RESLT_unsteady_heat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_unsteady_heat.dat.gz \
    result_unsteady_heat.dat  >> validation.log
fi


# Validation for demo Stefan Boltzmann problem
#---------------------------------------------

echo "Running 2D Stefan Boltzmann validation"
mkdir RESLT
../stefan_boltzmann > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D Stefan Boltzmann validation" >> validation.log
echo "------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_sb.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_sb.dat.gz \
    result_sb.dat  >> validation.log
fi

mv RESLT RESLT_sb


# Validation for demo unsteady heat with pretend melting
#-------------------------------------------------------

echo "Running 2D pretend melt validation "
mkdir RESLT
../pretend_melt --validate > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady pretend melt validation " >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_pretend_melt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_pretend_melt.dat.gz \
    result_pretend_melt.dat  >> validation.log
fi

mv RESLT RESLT_pretend_melt


# Validation for demo unsteady heat with pretend melting (2)
#-----------------------------------------------------------

echo "Running 2D pretend melt validation (2)"
mkdir RESLT
../melt --disable_melting --validate > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady pretend melt validation (2)" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_pretend_melt2.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_pretend_melt2.dat.gz \
    result_pretend_melt2.dat  >> validation.log
fi

mv RESLT RESLT_pretend_melt2




# Validation for demo unsteady heat with actual melting
#------------------------------------------------------

echo "Running 2D actual melt validation"
mkdir RESLT
../melt --validate > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D unsteady actual melt validation" >> validation.log
echo "----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_melt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_melt.dat.gz \
    result_melt.dat  >> validation.log
fi

mv RESLT RESLT_melt



# Validation for spring contact
#------------------------------

echo "Running spring contact (single kink) validation "
mkdir RESLT
cd RESLT
../../spring_contact --single_kink > OUTPUT
cd ..
echo "done"
echo " " >> validation.log
echo "Spring contact (single kink) validation " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_spring_contact_single_kink.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_spring_contact_single_kink.dat.gz \
    result_spring_contact_single_kink.dat  >> validation.log
fi

mv RESLT RESLT_spring_contact_single_kink


# Validation for spring contact
#------------------------------

echo "Running spring contact (kuhn tucker) validation "
mkdir RESLT
cd RESLT
../../spring_contact --kuhn_tucker > OUTPUT
cd ..
echo "done"
echo " " >> validation.log
echo "Spring contact (kuhn tucker) validation " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_spring_contact_kuhn_tucker.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_spring_contact_kuhn_tucker.dat.gz \
    result_spring_contact_kuhn_tucker.dat  >> validation.log
fi

mv RESLT RESLT_spring_contact_kuhn_tucker



# Validation for spring contact
#------------------------------

echo "Running spring contact (old version) validation "
mkdir RESLT
cd RESLT
../../spring_contact --old_version > OUTPUT
cd ..
echo "done"
echo " " >> validation.log
echo "Spring contact (old version) validation " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat  > result_spring_contact_old_version.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_spring_contact_old_version.dat.gz \
    result_spring_contact_old_version.dat  >> validation.log
fi

mv RESLT RESLT_spring_contact_old_version




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
