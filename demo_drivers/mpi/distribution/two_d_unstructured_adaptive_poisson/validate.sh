#! /bin/sh

sample_point_container_flag=""
if [ $# -eq 1 ]; then
    MPI_VARIABLENP_RUN_COMMAND="mpirun -np OOMPHNP "
    case $1 in
        --ref_bin ) sample_point_container_flag=" --ref_bin "
            ;;
        --non_ref_bin ) sample_point_container_flag=" --non_ref_bin "
            ;;
        --cgal ) sample_point_container_flag=" --cgal "
            ;;
        * ) echo "Wrong command line argument: $1"; exit 1
    esac
elif [ $# -gt 1 ]; then
    echo "Can only run this with one command line argument which "
    echo " is assumed to be sample point container flag: "
    echo " --ref_bin, --non_ref_bin or --cgal. "
    exit
fi
echo "sample point container flag: "$sample_point_container_flag


#Set the number of tests to be checked
NUM_TESTS=14


# Doc what we're using to run tests on variable processors
echo " " 
echo "Running mpi tests with mpi run command: " $MPI_VARIABLENP_RUN_COMMAND
echo "OOMPHNP = 1, 2, 3, 4"
echo " " 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation
#cp ../*partition.dat .

validation_flag=--validation

#----------------------------------------------------------------------

# Validation for unstructured adaptive two outer boundaries
#------------------------------------------------

echo "Running unstructured adaptive two outer boundaries validation "
mkdir RESLT


# Wait for a bit to allow creation of directory
sleep 5

../unstructured_adaptive_mesh_two_outer_boundaries $sample_point_container_flag > OUTPUT_unstructured_adaptive_mesh_two_outer_boundaries
echo "done"
echo " " >> validation.log
echo "unstructured_adaptive_mesh_two_outer_boundaries validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > TRACE_unstructured_adaptive_mesh_two_outer_boundaries_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_unstructured_adaptive_mesh_two_outer_boundaries_results.dat.gz  \
         TRACE_unstructured_adaptive_mesh_two_outer_boundaries_results.dat >> validation.log
fi

mv RESLT RESLT_UNSTRUCTURED_ADAPTIVE_MESH_TWO_OUTER_BOUNDARIES

#----------------------------------------------------------------------

# Validation for unstructured adaptive two outer boundaries crossed
#------------------------------------------------

echo "Running unstructured adaptive two outer boundaries crossed validation "
mkdir RESLT


# Wait for a bit to allow creation of directory
sleep 5

../unstructured_adaptive_mesh_two_outer_boundaries_crossed  $sample_point_container_flag > OUTPUT_unstructured_adaptive_mesh_two_outer_boundaries_crossed
echo "done"
echo " " >> validation.log
echo "unstructured_adaptive_mesh_two_outer_boundaries_crossed validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > TRACE_unstructured_adaptive_mesh_two_outer_boundaries_crossed_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_unstructured_adaptive_mesh_two_outer_boundaries_crossed_results.dat.gz  \
         TRACE_unstructured_adaptive_mesh_two_outer_boundaries_crossed_results.dat >> validation.log
fi

mv RESLT RESLT_UNSTRUCTURED_ADAPTIVE_MESH_TWO_OUTER_BOUNDARIES_CROSSED

# --------------------------------------------------------------------------
# Validation flags for SQUARE domain
# --------------------------------------------------------------------------
domain_configuration_flag=--domain_configuration
domain_configuration_value=1
element_size_flag=--element_size
element_size_value=0.01
max_adapt_flag=--max_adapt
max_adapt_value=3
max_permitted_error_flag=--max_permitted_error
max_permitted_error_value=1.0e-3
min_permitted_error_flag=--min_permitted_error
min_permitted_error_value=1.0e-5
max_element_size_flag=--max_element_size
max_element_size_value=$element_size_value
min_element_size_flag=--min_element_size
min_element_size_value=1.0e-14
load_balance_flag=--load_balance
load_balance_value=0

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (2 processors) (SQUARE domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_SQUARE

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 2 processors for SQUARE domain"
mkdir RESLT

# two processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/2/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_square ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_square_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (2 processors) SQUARE domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_square_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_square_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_square_results.dat 8.0 1.0e-10 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_SQUARE

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (3 processors) (SQUARE domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_SQUARE

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 3 processors for SQUARE domain"
mkdir RESLT

# three processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/3/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_square ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_square_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (3 processors) SQUARE domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_square_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_square_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_square_results.dat 5.0 1.0e-10 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_SQUARE

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (4 processors) (SQUARE domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_SQUARE

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 4 processors for SQUARE domain"
mkdir RESLT

# four processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/4/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_square ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_square_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (4 processors) SQUARE domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_square_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_square_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_square_results.dat 6.5 1.0e-10 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_SQUARE

# --------------------------------------------------------------------------
# Validation flags for HALF CIRCLE domain
# --------------------------------------------------------------------------
domain_configuration_flag=--domain_configuration
domain_configuration_value=2
element_size_flag=--element_size
element_size_value=0.05
max_adapt_flag=--max_adapt
max_adapt_value=3
max_permitted_error_flag=--max_permitted_error
max_permitted_error_value=1.0e-3
min_permitted_error_flag=--min_permitted_error
min_permitted_error_value=1.0e-5
max_element_size_flag=--max_element_size
max_element_size_value=$element_size_value
min_element_size_flag=--min_element_size
min_element_size_value=1.0e-14
load_balance_flag=--load_balance
load_balance_value=0

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (2 processors) (HALF CIRCLE domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 2 processors for HALF CIRCLE domain"
mkdir RESLT

# two processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/2/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (2 processors) HALF CIRCLE domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_results.dat 5.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (3 processors) (HALF CIRCLE domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 3 processors for HALF CIRCLE domain"
mkdir RESLT

# three processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/3/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (3 processors) HALF CIRCLE domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (4 processors) (HALF CIRCLE domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 4 processors for HALF CIRCLE domain"
mkdir RESLT

# four processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/4/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (4 processors) HALF CIRCLE domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_results.dat 5.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE

# --------------------------------------------------------------------------
# Validation flags for HALF CIRCLE domain with INTERNAL BOUNDARIES
# --------------------------------------------------------------------------
domain_configuration_flag=--domain_configuration
domain_configuration_value=3
element_size_flag=--element_size
element_size_value=0.1
max_adapt_flag=--max_adapt
max_adapt_value=3
max_permitted_error_flag=--max_permitted_error
max_permitted_error_value=1.0e-3
min_permitted_error_flag=--min_permitted_error
min_permitted_error_value=1.0e-5
max_element_size_flag=--max_element_size
max_element_size_value=$element_size_value
min_element_size_flag=--min_element_size
min_element_size_value=1.0e-14
load_balance_flag=--load_balance
load_balance_value=0

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (2 processors) (HALF_CIRCLE_INTERNAL_BOUNDARIES domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE_INTERNAL_BOUNDARIES

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 2 processors for HALF_CIRCLE_INTERNAL_BOUNDARIES domain"
mkdir RESLT

# two processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/2/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_internal_boundaries ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_internal_boundaries_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (2 processors) HALF_CIRCLE_INTERNAL_BOUNDARIES domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_internal_boundaries_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_internal_boundaries_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_half_circle_internal_boundaries_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE_INTERNAL_BOUNDARIES

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (3 processors) (HALF_CIRCLE_INTERNAL_BOUNDARIES domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE_INTERNAL_BOUNDARIES

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 3 processors for HALF_CIRCLE_INTERNAL_BOUNDARIES domain"
mkdir RESLT

# three processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/3/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_internal_boundaries ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_internal_boundaries_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (3 processors) HALF_CIRCLE_INTERNAL_BOUNDARIES domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_internal_boundaries_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_internal_boundaries_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_half_circle_internal_boundaries_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE_INTERNAL_BOUNDARIES

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (4 processors) (HALF_CIRCLE_INTERNAL_BOUNDARIES domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE_INTERNAL_BOUNDARIES

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 4 processors for HALF_CIRCLE_INTERNAL_BOUNDARIES domain"
mkdir RESLT

# four processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/4/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_internal_boundaries ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_internal_boundaries_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (4 processors) HALF_CIRCLE_INTERNAL_BOUNDARIES domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_internal_boundaries_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_internal_boundaries_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_half_circle_internal_boundaries_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE_INTERNAL_BOUNDARIES

# --------------------------------------------------------------------------
# Validation flags for COMPLEX HOLES domain
# --------------------------------------------------------------------------
domain_configuration_flag=--domain_configuration
domain_configuration_value=4
element_size_flag=--element_size
element_size_value=0.1
max_adapt_flag=--max_adapt
max_adapt_value=3
max_permitted_error_flag=--max_permitted_error
max_permitted_error_value=1.0e-3
min_permitted_error_flag=--min_permitted_error
min_permitted_error_value=1.0e-5
max_element_size_flag=--max_element_size
max_element_size_value=$element_size_value
min_element_size_flag=--min_element_size
min_element_size_value=1.0e-14
load_balance_flag=--load_balance
load_balance_value=0

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (2 processors) (COMPLEX_HOLES domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_COMPLEX_HOLES

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 2 processors for COMPLEX_HOLES domain"
mkdir RESLT

# two processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/2/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_complex_holes ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np2_complex_holes_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (2 processors) COMPLEX_HOLES domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_complex_holes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_complex_holes_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np2_complex_holes_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_COMPLEX_HOLES

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (3 processors) (COMPLEX_HOLES domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_COMPLEX_HOLES

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 3 processors for COMPLEX_HOLES domain"
mkdir RESLT

# three processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/3/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_complex_holes ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np3_complex_holes_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (3 processors) COMPLEX_HOLES domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_complex_holes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_complex_holes_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np3_complex_holes_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_COMPLEX_HOLES

#------------------------------------------------
# Validation for two_d_parallel_unstructured_adaptive_poisson (4 processors) (COMPLEX_HOLES domain)
#------------------------------------------------

folder_distribution_file_flag=--folder_distribution_file
folder_distribution_file_value=../validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_COMPLEX_HOLES

echo "Running two_d_parallel_unstructured_adaptive_poisson validation with 4 processors for COMPLEX_HOLES domain"
mkdir RESLT

# four processor
MPI_RUN_ON_NP_COMMAND=`echo $MPI_VARIABLENP_RUN_COMMAND | sed -e "s/OOMPHNP/4/g"`
echo "Flags current state"
echo $domain_configuration_flag $domain_configuration_value
echo $element_size_flag $element_size_value
echo $max_adapt_flag $max_adapt_value
echo $max_permitted_error_flag $max_permitted_error_value
echo $min_permitted_error_flag $min_permitted_error_value
echo $max_element_size_flag $max_element_size_value
echo $min_element_size_flag $min_element_size_value
echo $load_balance_flag $load_balance_value
echo $folder_distribution_file_flag $folder_distribution_file_value

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

$MPI_RUN_ON_NP_COMMAND --output-filename OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_complex_holes ../two_d_parallel_unstructured_adaptive_poisson $sample_point_container_flag $domain_configuration_flag $domain_configuration_value $element_size_flag $element_size_value $max_adapt_flag $max_adapt_value $max_permitted_error_flag $max_permitted_error_value $min_permitted_error_flag $min_permitted_error_value $max_element_size_flag $max_element_size_value $min_element_size_flag $min_element_size_value $load_balance_flag $load_balance_value $folder_distribution_file_flag $folder_distribution_file_value > OUTPUT_two_d_parallel_unstructured_adaptive_poisson_np4_complex_holes_from_screen
echo "done"
echo " " >> validation.log
echo "two_d_parallel_unstructured_adaptive_poisson (4 processors) COMPLEX_HOLES domain validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_proc0.dat > TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_complex_holes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_complex_holes_results.dat.gz  \
         TRACE_two_d_parallel_unstructured_adaptive_poisson_np4_complex_holes_results.dat 3.0 1.0e-12 >> validation.log
fi

mv RESLT RESLT_TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_COMPLEX_HOLES

#----------------------------------------------------------------------

# Append log to main validation log
cat validation.log >> ../../../../../validation.log

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
  fi
fi
