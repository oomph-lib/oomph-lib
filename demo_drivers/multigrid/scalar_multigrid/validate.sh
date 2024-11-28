#!/bin/bash

# Get the OOPMH-LIB root directory from a makefile
#-------------------------------------------------
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# The base directory
#-------------------
base_dir=`pwd`;

# Set the number of tests to be checked
#--------------------------------------
NUM_TESTS=148

# Which code do you want to run?
#-------------------------------
code_stem=(two_d_poisson_tanh_flux_bc_validate unit_cube_poisson_validate);
demo_code_stem=(two_d_poisson_tanh_flux_bc unit_cube_poisson);

# Name the directory to work in
#------------------------------
val_dir="Validation"

# Name the storage directory
#---------------------------
store="validata"

# Delete it if it already exists
#-------------------------------
# If the directory already exists
if [ -d $val_dir ];
then
    # If the directory exists, kill it
    rm -rf $val_dir
fi

# Make it
#--------
mkdir $val_dir;

# Make the validata directory
#----------------------------
# If the directory doesn't exist
if [ ! -d $store ];
then
    # If the directory doesn't exist, make it
    mkdir $store
fi

# Jump in to the directory
#-------------------------
cd $val_dir

# In the Validation directory we need a validation.log file; create it!
#----------------------------------------------------------------------
touch validation.log

# Copy the executable into this directory
#----------------------------------------
for code in ${code_stem[@]}
do
    cp ../$code ./
done
for code in ${demo_code_stem[@]}
do
    cp ../$code ./
done

# Results directory
#------------------
result_dir="RESLT";

# Set array of polynomial degrees to be used
#-------------------------------------------
nnode_1d_list=(2 3 4);

# Refinement levels
#------------------
min_refine=1;
add_refine=(1);

# Number of adaptations to be used
#---------------------------------
n_adapt=(1);

# Refinement array
#-----------------
# Needs to be the same length as add_refine and n_adapt
# (used to avoid vagueness in the innermost for loop)
refine_array=(0);

# Whether or not we use adaptative refinement
#--------------------------------------------
use_adapt=(0 1);

# Linear solver
#--------------
# We only use the MG solver (not SuperLU)
solver=1;

# Ask for the convergence information
#------------------------------------
conv_flag=1;

# Pre- and post-smoother list
#----------------------------
presmoother_list=(0 1);
postsmoother_list=(0 1);

# Set the output file name
#-------------------------
output_file=output_file.txt;

# Delete it if it already exists
#-------------------------------
if [ -e $output_file ]
then
    rm -f $output_file
fi

# Create the output file
#-----------------------
touch $output_file

# We will have three files to move
#---------------------------------
file_0="/soln0.dat";
file_1="/soln1.dat";
file_2="/conv.dat";	

# Create the gzipped filename suffix
#-----------------------------------
gzip_suffix=".gz";

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
for code in ${code_stem[@]}
do
    if [ $code == ${code_stem[0]} ]
    then
	echo "Running multigrid validation with first 2D Poisson problem"
    else
	echo "Running multigrid validation with first 3D Poisson problem"
    fi
    for nnode in ${nnode_1d_list[@]}
    do
	for pre in ${presmoother_list[@]}
	do
	    for post in ${postsmoother_list[@]}
	    do
		for adapt in ${use_adapt[@]}
		do
		    # Make it
		    #--------
		    if [ ! -d $result_dir ];
		    then
			mkdir $result_dir;
		    fi
		    
		    for ref in ${refine_array[@]}
		    do
			# Set the inputs
			#---------------
			command_line_flag="--min_ref $min_refine"
			command_line_flag+=" --linear_solver $solver"
			command_line_flag+=" --use_adapt $adapt"
			command_line_flag+=" --n_adapt ${n_adapt[$ref]}"
			command_line_flag+=" --add_ref ${add_refine[$ref]}"
			command_line_flag+=" --nnode_1d $nnode "
			command_line_flag+=" --presmoother $pre"
			command_line_flag+=" --postsmoother $post"
			command_line_flag+=" --conv_flag $conv_flag"
			
			# Run the code and direct the output into the output file
			#--------------------------------------------------------
			command_line_prompt="./$code $command_line_flag"
			`echo $command_line_prompt` >> $output_file
			
			# Create a little pause on the 3D problem
			#----------------------------------------
			if [ $code != ${code_stem[0]} ]
			then
			    sleep 3s
			fi
		    done

		    # Create some output in the validation.log file
		    echo " " >> validation.log		    
		    if [ $code == ${code_stem[0]} ]
		    then
			echo "Multigrid validation with 2D Poisson problem" >> validation.log
		    else
			echo "Multigrid validation with 3D Poisson problem" >> validation.log
		    fi		    
		    echo "--------------------------------------------" >> validation.log
		    echo "Problem properties:" >> validation.log

		    # Interpolation scheme
		    string="Interpolation scheme:"
		    if [ $nnode == 2 ]
		    then
			echo "$string Linear" >> validation.log
		    elif [ $nnode == 3 ]
		    then
			echo "$string Quadratic" >> validation.log
		    else
			echo "$string Cubic" >> validation.log		
		    fi

		    # Pre-smoother
		    string="Pre-smoother:"
		    if [ $pre == 0 ]
		    then
			echo "$string Damped Jacobi" >> validation.log
		    else
			echo "$string Gauss-Seidel" >> validation.log
		    fi
		    
		    # Post-smoother
		    string="Post-smoother:"
		    if [ $pre == 0 ]
		    then
			echo "$string Damped Jacobi" >> validation.log
		    else
			echo "$string Gauss-Seidel" >> validation.log
		    fi
		    
		    # Refinement scheme
		    string="Refinement scheme:"
		    if [ $adapt == 0 ]
		    then
			echo "$string Uniform refinement" >> validation.log
		    else
			echo "$string Adaptive refinement" >> validation.log
		    fi
		    
		    echo "--------------------------------------" >> validation.log
		    echo " " >> validation.log
		    echo "Validation directory: " >> validation.log
		    echo " " >> validation.log
		    echo "  " `pwd` >> validation.log
		    echo " " >> validation.log
		    
		    # Create the unique output file using the for loop names
		    #-------------------------------------------------------
		    filename="/$code";
		    filename+="_nnode_"$nnode;
		    filename+="_pre_"$pre;
		    filename+="_post_"$post;
		    filename+="_adapt_"$adapt;
		    
		    # We need two unique filenames
		    #-----------------------------
		    validata_file_0=$filename"_0.dat.gz"
		    validata_file_1=$filename"_1.dat.gz"
		    validata_file_2=$filename"_conv.dat"
		    
		    # Check conv.dat
		    #----------------
		    # Compare number of iterations against reference data and append
		    ./$compare_length $result_dir$file_2 ../$store$validata_file_2 $threshold_n_its >> validation.log
		    
		    # Check soln0.dat
		    #----------------
		    if test "$1" = "no_fpdiff";
		    then
			echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
		    else
			./$compare_data ../$store$validata_file_0 $result_dir$file_0 0.1 1e-13 >> validation.log
		    fi
		    
		    # Check soln1.dat
		    #----------------
		    if test "$1" = "no_fpdiff";
		    then
			echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
		    else
			./$compare_data ../$store$validata_file_1 $result_dir$file_1 0.1 1e-13 >> validation.log
		    fi
		    
		    # Move the result directory into storage
		    #---------------------------------------
		    mv $result_dir ./$result_dir"_"$code$nnode$pre$post$adapt
		done
	    done
	done    
    done
    # Notify the user that we've finished the test
    echo "done"
done

# Run the program with the specified properties
#----------------------------------------------
for code in ${demo_code_stem[@]}
do
    if [ $code == ${demo_code_stem[0]} ]
    then
	echo "Running multigrid validation with second 2D Poisson problem"
    else
	echo "Running multigrid validation with second 3D Poisson problem"
    fi
    
    # Make it
    #--------
    if [ ! -d $result_dir ];
    then
	mkdir $result_dir;
    fi
    
    # Run the code and direct the output into the output file
    #-------------------------------------------------------
    command_line_prompt="./$code $command_line_flag"
    `echo $command_line_prompt` >> $output_file
    
    # Create the unique output file using the for loop names
    #-------------------------------------------------------
    filename="/$code";
    
    # We need two unique filenames
    #-----------------------------
    validata_file_0=$filename"_0.dat.gz"
    validata_file_1=$filename"_conv.dat"

    # Create some output in the validation.log file
    echo " " >> validation.log		    
    if [ $code == ${demo_code_stem[0]} ]
    then
	echo "Multigrid validation with 2D Poisson demo problem" >> validation.log
    else
	echo "Multigrid validation with 3D Poisson demo problem" >> validation.log
    fi		    
    echo "--------------------------------------------" >> validation.log
    echo "Problem properties:" >> validation.log

    # Interpolation scheme
    string="Interpolation scheme:"
    echo "$string Quadratic" >> validation.log

    # Pre-smoother
    string="Pre-smoother:"
    echo "$string Gauss-Seidel" >> validation.log
    
    # Post-smoother
    string="Post-smoother:"
    echo "$string Gauss-Seidel" >> validation.log
    
    # Refinement scheme
    string="Refinement scheme:"
    echo "$string Uniform refinement" >> validation.log
    
    echo "--------------------------------------" >> validation.log
    echo " " >> validation.log
    echo "Validation directory: " >> validation.log
    echo " " >> validation.log
    echo "  " `pwd` >> validation.log
    echo " " >> validation.log
    
    # Check conv.dat
    #----------------
    # Compare number of iterations against reference data and append
    ./$compare_length $result_dir$file_2 ../$store$validata_file_1 $threshold_n_its >> validation.log
    
    # Check soln0.dat
    #----------------
    if test "$1" = "no_fpdiff";
    then
	echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    else
	./$compare_data ../$store$validata_file_0 $result_dir$file_0 0.1 1e-13 >> validation.log
    fi
    
    # Move the result directory into storage
    #---------------------------------------
    mv $result_dir ./$result_dir"_"$code$nnode$pre$post$adapt
    
    # Notify the user that we've finished the test
    echo "done"
done

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log

# Jump to the base directory
#---------------------------
cd ..

# Delete the useless results directory
#-------------------------------------
rm -rf $result_dir

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
