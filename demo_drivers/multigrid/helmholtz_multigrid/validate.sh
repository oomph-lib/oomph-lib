#!/bin/bash

# Get the OOMPH-LIB root directory from a makefile
#-------------------------------------------------
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# The base directory
#-------------------
base_dir=$(pwd)

# Set the number of tests to be checked (checking convergence AND trace file)
#----------------------------------------------------------------------------
NUM_TESTS=4

# Name the directory to work in
#------------------------------
temp_dir="Validation"

# Name the storage directory
#---------------------------
store="validata"

# Results directory
#------------------
result_dir="RESLT"

# Make the working directory
#---------------------------
# If the directory doesn't exist
if [ -d $temp_dir ]; then
    rm -rf $temp_dir
fi

# Now make it
mkdir $temp_dir

# Make the validata directory
#----------------------------
# If the directory doesn't exist
if [ ! -d $store ]; then
    # If the directory doesn't exist, make it
    mkdir $store
fi

# Which code do you want to run?
#-------------------------------
code_stem=(unstructured_two_d_helmholtz structured_cubic_point_source)

# Jump in to the validation directory
#------------------------------------
cd $temp_dir

# Copy the executable into this directory
#----------------------------------------
for code in ${code_stem[@]}; do
    cp ../$code ./
done

# Set array of polynomial degrees to be used
#-------------------------------------------
nnode_1d=2

# Refinement levels
#------------------
add_refine=0
min_refine=1

# Linear solver
#--------------
# For the purposes of validation we only use the MG solver
linear_solver=1

# MG smoother
#------------
pre_smoother=1
post_smoother=1

# Value of k^2
#-------------
k_sq=0.0

# Number of PML elements
#-----------------------
npml_element=1

# Use the test PML mapping?
#--------------------------
test_pml_mapping=1

# Alpha value
#------------
alpha=0.5

# Ask for the convergence information
#------------------------------------
conv_flag=1

# Set the output file name
#-------------------------
output_file=output_file.txt

# Delete it if it already exists
#-------------------------------
if [ -e $output_file ]; then
    rm -f $output_file
fi

# Create the output file
#-----------------------
touch $output_file

# Create the gzipped filename suffix
#-----------------------------------
gzip_suffix=".gz"

# Create the file names:
#-----------------------
# The convergence file
conv_file="/conv.dat"

# The trace file
trace_file="/trace.dat"

# Threshold for number of iterations in comparison of convergence histories:
#---------------------------------------------------------------------------
threshold_n_its=3

# The address of the bash script to compare file lengths
#-------------------------------------------------------
compare_length="../../../../bin/compare_file_length_with_tolerance.bash"

# The address of the bash script to compare file data
#----------------------------------------------------
compare_data="../../../../bin/fpdiff.py"

# Run the program with the specified properties
#----------------------------------------------
for code in ${code_stem[@]}; do
    if [ $code == ${code_stem[0]} ]; then
        echo "Running multigrid validation with 2D Helmholtz problem"
    else
        echo "Running multigrid validation with 3D Helmholtz problem"
    fi

    # Make it
    #--------
    if [ ! -d $result_dir ]; then
        mkdir $result_dir
    fi

    # If we're working in 2D we can deal with a higher value of k^2
    # but it takes much longer in 3D so use a lower wavenumber in 3D
    #---------------------------------------------------------------
    if [ $code == ${code_stem[0]} ]; then
        k_sq=30.0
    else
        k_sq=15.0
    fi

    # If we're working in 2D use quadratic elements otherwise use linear
    # elements since the cost of quadratic elements in 3D is extremely large
    #-----------------------------------------------------------------------
    if [ $code == ${code_stem[0]} ]; then
        nnode_1d=3
    else
        nnode_1d=2
    fi

    # Set the inputs
    #---------------
    command_line_flag="--add_ref $add_refine "
    command_line_flag+="--min_ref $min_refine "
    command_line_flag+="--linear_solver $linear_solver "
    command_line_flag+="--alpha $alpha "
    command_line_flag+="--k_sq $k_sq "
    command_line_flag+="--nnode_1d $nnode_1d "
    command_line_flag+="--conv_flag $conv_flag "
    command_line_flag+="--presmoother $pre_smoother "
    command_line_flag+="--postsmoother $post_smoother "
    command_line_flag+="--test_pml_mapping $test_pml_mapping "

    # Extra inputs which are only relevant to one specific code
    #----------------------------------------------------------
    if [ $code == ${code_stem[0]} ]; then
        # Can only set the number of PML layers in the 3D case
        command_line_flag+="--npml_element $npml_element "
    fi

    # Run the code and direct the output into the output file
    #--------------------------------------------------------
    ./$code $command_line_flag >>$output_file

    # Create a little pause on the 3D problem
    #----------------------------------------
    if ((${#code_stem[@]} > 1)); then
        if [ $code == ${code_stem[1]} ]; then
            sleep 3s
        fi
    fi

    # Create some output in the validation.log file:
    #-----------------------------------------------
    echo " " >>validation.log
    echo "-------------------------------------------------------------" >>validation.log
    if [ $code == ${code_stem[0]} ]; then
        echo "Complex-valued multigrid validation with 2D Helmholtz problem" >>validation.log
    else
        echo "Complex-valued multigrid validation with 3D Helmholtz problem" >>validation.log
    fi
    echo "-------------------------------------------------------------" >>validation.log
    echo "Problem properties:" >>validation.log

    # Interpolation scheme:
    #----------------------
    string="Interpolation scheme:"
    if [ $nnode_1d == 2 ]; then
        echo "$string Linear" >>validation.log
    elif [ $nnode_1d == 3 ]; then
        echo "$string Quadratic" >>validation.log
    else
        echo "$string Cubic" >>validation.log
    fi

    # Value of the square of the wavenumber:
    #---------------------------------------
    echo "Value of k^2: $k_sq" >>validation.log

    # Value of the shift:
    #--------------------
    echo "Value of alpha: $alpha" >>validation.log

    # Refinement scheme:
    #-------------------
    echo "Refinement scheme: Uniform refinement" >>validation.log
    echo "--------------------------------------" >>validation.log
    echo " " >>validation.log
    echo "Validation directory: " >>validation.log
    echo " " >>validation.log
    echo "  " $(pwd) >>validation.log
    echo " " >>validation.log

    # Create the unique output file using the identifying parameters
    #---------------------------------------------------------------
    filename="/$code"
    filename+="_ksq_"$k_sq
    filename+="_alpha_"$alpha

    # We need a unique filename for the convergence history
    #------------------------------------------------------
    conv_filename=$filename"_conv.dat"
    trace_filename=$filename"_trace.dat.gz"

    # Check trace.dat and conv.dat
    #-----------------------------
    if test "$1" = "no_fpdiff"; then
        echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
    else
        ./$compare_data ../$store$trace_filename $result_dir$trace_file >>validation.log

        # Compare number of iterations against reference data and append
        ./$compare_length $result_dir$conv_file ../$store$conv_filename $threshold_n_its >>validation.log
    fi

    # Move the result directory into storage
    #---------------------------------------
    new_result_dir=$result_dir"_"$code
    mv $result_dir $new_result_dir

    # Notify the user that we've finished the test
    echo "done"
done

# Append output to global validation log file
#--------------------------------------------
cat validation.log >>../../../../validation.log

# Now jump back to the base directory
#------------------------------------
cd $base_dir

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
