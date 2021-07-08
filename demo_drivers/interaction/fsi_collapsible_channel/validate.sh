#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=12

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Validation for FSI collapsible channel (non-adapt, pseudo-elastic node update)
#--------------------------------------------------------------------
cd Validation


echo "Running validation for FSI collapsible channel (pseudo elastic, non-adapt TH)"
mkdir RESLT

# Do validation run
../fsi_pseudo_solid_collapsible_channel_TH blabla > OUTPUT_pseudo_elastic_TH
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (pseudo elastic, non-adapt TH)" \
>> validation.log
echo "-----------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_pseudo_elastic_non_adapt_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pseudo_elastic_non_adapt_TH.dat.gz \
    result_pseudo_elastic_non_adapt_TH.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_pseudo_elastic_non_adapt_TH


echo "Running validation for FSI collapsible channel (pseudo elastic, non-adapt CR)"
mkdir RESLT

# Do validation run
../fsi_pseudo_solid_collapsible_channel_CR blabla > OUTPUT_pseudo_elastic_CR

echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (pseudo elastic, non-adapt CR)" \
>> validation.log
echo "-----------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_pseudo_elastic_non_adapt_CR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pseudo_elastic_non_adapt_CR.dat.gz \
    result_pseudo_elastic_non_adapt_CR.dat 0.1 1.0e-8 >> validation.log
fi

mv RESLT RESLT_pseudo_elastic_non_adapt_CR



# Validation for FSI collapsible channel (pseudo-elastic node update)
#--------------------------------------------------------------------

echo "Running validation for FSI collapsible channel (pseudo elastic, TH)"
mkdir RESLT

# Do validation run
../fsi_pseudo_solid_collapsible_channel_adapt_TH blabla > OUTPUT_pseudo_elastic_adapt_TH
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (pseudo elastic, TH)" \
>> validation.log
echo "-------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_pseudo_elastic_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pseudo_elastic_TH.dat.gz \
    result_pseudo_elastic_TH.dat 0.5 1.0e-8 >> validation.log
fi


mv RESLT RESLT_pseudo_elastic_TH


echo "Running validation for FSI collapsible channel (pseudo elastic, CR)"
mkdir RESLT

# Do validation run
../fsi_pseudo_solid_collapsible_channel_adapt_CR blabla > OUTPUT_pseudo_elastic_adapt_CR

echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (pseudo elastic, CR)" \
>> validation.log
echo "-------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_pseudo_elastic_CR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pseudo_elastic_CR.dat.gz \
    result_pseudo_elastic_CR.dat 7.0 1.50e-6 >> validation.log
fi

mv RESLT RESLT_pseudo_elastic_CR


# Validation for FSI collapsible channel (macro-element update)
#--------------------------------------------------------------------


echo "Running validation for FSI collapsible channel (macro elements, TH)"
mkdir RESLT

# Do validation run
../fsi_collapsible_channel_macro_TH blabla > OUTPUT_macro_no_adapt_TH
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (macro elements, TH)" \
>> validation.log
echo "---------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_macro_no_adapt_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/macro_no_adapt_TH.dat.gz \
    result_macro_no_adapt_TH.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_macro_no_adapt_TH

echo "Running validation for FSI collapsible channel (macro elements, CR)"
mkdir RESLT

# Do validation run
../fsi_collapsible_channel_macro_CR blabla > OUTPUT_macro_no_adapt_CR
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (macro elements, CR)" \
>> validation.log
echo "---------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_macro_no_adapt_CR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/macro_no_adapt_CR.dat.gz \
    result_macro_no_adapt_CR.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_macro_no_adapt_CR



# Validation for FSI collapsible channel (algebraic elements)
#------------------------------------------------------------


echo "Running validation for FSI collapsible channel (alg elements, TH)"
mkdir RESLT

# Do validation run
../fsi_collapsible_channel_algebraic_TH blabla > OUTPUT_algebraic_no_adapt_TH
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (alg elements, TH)" \
>> validation.log
echo "---------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_algebraic_no_adapt_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg_no_adapt_TH.dat.gz \
    result_algebraic_no_adapt_TH.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_algebraic_no_adapt_TH


echo "Running validation for FSI collapsible channel (alg elements, CR)"
mkdir RESLT

# Do validation run
../fsi_collapsible_channel_algebraic_CR blabla > OUTPUT_algebraic_no_adapt_CR
echo "done"
echo " " >> validation.log
echo "FSI collapsible channel validation (alg elements, CR)" \
>> validation.log
echo "---------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_algebraic_no_adapt_CR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg_no_adapt_CR.dat.gz \
    result_algebraic_no_adapt_CR.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_algebraic_no_adapt_CR




# Validation for adaptive FSI collapsible channel (macro elements)
#-----------------------------------------------------------------

echo "Running validation for adaptive FSI collapsible channel (macro elements, TH)"
mkdir RESLT 

# Do validation run
../fsi_collapsible_channel_macro_adapt_TH blabla > OUTPUT_macro_adapt_TH
echo "done"
echo " " >> validation.log
echo "Adaptive FSI collapsible channel validation (macro elements, TH)" \
>> validation.log
echo "------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_macro_adapt_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/macro_adapt_TH.dat.gz \
    result_macro_adapt_TH.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_macro_adapt_TH



echo "Running validation for adaptive FSI collapsible channel (macro elements, CR)"
mkdir RESLT 

# Do validation run
../fsi_collapsible_channel_macro_adapt_CR blabla > OUTPUT_macro_adapt_CR
echo "done"
echo " " >> validation.log
echo "Adaptive FSI collapsible channel validation (macro elements, CR)" \
>> validation.log
echo "------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_macro_adapt_CR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/macro_adapt_CR.dat.gz \
    result_macro_adapt_CR.dat 7.0 1.5e-6 >> validation.log
fi


mv RESLT RESLT_macro_adapt_CR


# Validation for adaptive FSI collapsible channel (algebraic elements)
#-----------------------------------------------------------------

echo "Running validation for adaptive FSI collapsible channel (algebraic elements, TH)"
mkdir RESLT 

# Do validation run
../fsi_collapsible_channel_algebraic_adapt_TH blabla > OUTPUT_algebraic_adapt_TH
echo "done"
echo " " >> validation.log
echo "Adaptive FSI collapsible channel validation (algebraic elements, TH)" \
>> validation.log
echo "------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_algebraic_adapt_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg_adapt_TH.dat.gz \
    result_algebraic_adapt_TH.dat 0.1 1.0e-8 >> validation.log
fi


mv RESLT RESLT_algebraic_adapt_TH



echo "Running validation for adaptive FSI collapsible channel (algebraic elements, CR)"
mkdir RESLT 

# Do validation run
../fsi_collapsible_channel_algebraic_adapt_CR blabla > OUTPUT_algebraic_adapt_CR
echo "done"
echo " " >> validation.log
echo "Adaptive FSI collapsible channel validation (algebraic elements, CR)" \
>> validation.log
echo "------------------------------------------------------------" \
>> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln1.dat  \
    RESLT/soln2.dat RESLT/soln3.dat  \
    > result_algebraic_adapt_CR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg_adapt_CR.dat.gz \
    result_algebraic_adapt_CR.dat 7.0 1.5e-6 >> validation.log
fi


mv RESLT RESLT_algebraic_adapt_CR



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
