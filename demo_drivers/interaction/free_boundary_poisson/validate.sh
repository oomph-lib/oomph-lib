#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=7


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Validation for Circle
#----------------------
cd Validation

echo "Running circle validation "
../circle > OUTPUT_circle
echo "done"
echo " " >> validation.log
echo "Circle validation" >> validation.log
echo "-----------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat circle0.dat circle1.dat > circle_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/circle_results.dat.gz \
    circle_results.dat  >> validation.log
fi


# Validation for Circle as GeneralisedElement
#--------------------------------------------

echo "Running circle as GeneralisedElement validation "
mkdir RESLT
../geom_object_element > OUTPUT_circle_element
echo "done"
echo " " >> validation.log
echo "Circle as GeneralisedElement validation" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > circle_element_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/circle_element_results.dat.gz \
    circle_element_results.dat  >> validation.log
fi
mv RESLT RESLT_circle


# Validation for doc of macro-element-based node update
#------------------------------------------------------

echo "Running doc of macro-element-based node update validation "
mkdir RESLT
../doc_sparse_macro_node_update > OUTPUT_doc_sparse_macro_node_update
echo "done"
echo " " >> validation.log
echo "Doc of macro-element-based node update validation" >> validation.log
echo "-------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln10.dat RESLT/soln20.dat \
  > doc_sparse.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/doc_sparse.dat.gz \
    doc_sparse.dat  >> validation.log
fi
mv RESLT RESLT_doc_sparse_macro_node_update



# Validation for free_boundary poisson with algebraic mesh update
#----------------------------------------------------------

echo "Running free_boundary poisson with algebraic mesh update validation "
mkdir RESLT
mkdir RESLT_coupled
../algebraic_free_boundary_poisson lalal > OUTPUT_algebraic_free_boundary_poisson
echo "done"
echo " " >> validation.log
echo "Free-boundary poisson with algebraic mesh update validation" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT_coupled/soln0.dat RESLT_coupled/soln_1_0.dat  RESLT_coupled/soln_2_0.dat   RESLT_coupled/soln_3_0.dat   RESLT_coupled/soln_4_0.dat  > alg.dat
mv RESLT RESLT_algebraic
mv RESLT_coupled RESLT_coupled_algebraic

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/alg.dat.gz\
     alg.dat >> validation.log
fi




# Validation for free_boundary poisson with macro element mesh update
#--------------------------------------------------------------------

echo "Running free_boundary poisson with macro element mesh update validation "
mkdir RESLT
../macro_element_free_boundary_poisson lalal > OUTPUT_macro_element_free_boundary_poisson
echo "done"
echo " " >> validation.log
echo "Free-boundary poisson with macro_element mesh update validation" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  > macro.dat
mv RESLT RESLT_coupled_macro

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/macro.dat.gz\
     macro.dat >> validation.log
fi




# Validation for old free_boundary poisson with macro element mesh update
#--------------------------------------------------------------------

echo "Running old free_boundary poisson with macro element mesh update validation "
mkdir RESLT
mkdir RESLT_coupled
../old_macro_element_free_boundary_poisson_for_doc lalal > OUTPUT_old_macro_element_free_boundary_poisson_for_doc
echo "done"
echo " " >> validation.log
echo "Old free-boundary poisson with macro_element mesh update validation" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT_coupled/soln0.dat > old_macro.dat
mv RESLT RESLT_old_macro
mv RESLT_coupled RESLT_coupled_old_macro

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/old_macro.dat.gz\
     old_macro.dat >> validation.log
fi


# Validation for free_bound poisson w/ macro el mesh update (non-refineable)
#---------------------------------------------------------------------------

echo "Running free_boundary poisson with macro element mesh update "
echo "(non-refineable) validation "
mkdir RESLT
../macro_element_free_boundary_poisson_non_ref lalal > OUTPUT_macro_element_free_boundary_poisson_non_ref
echo "done"
echo " " >> validation.log
echo "Free-boundary poisson w/ macro_el mesh update (non-ref) validation" >> validation.log
echo "------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  > macro_non_ref.dat
mv RESLT RESLT_coupled_macro_non_ref

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/macro_non_ref.dat.gz\
     macro_non_ref.dat >> validation.log
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
