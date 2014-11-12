#! /bin/bash

. bin/autogen_helpers.sh

set -o errexit


#====================================================================
# Start Q/A session
#====================================================================

echo " "
echo "============================================================= "
echo "              oomph-lib installation script" 
echo "============================================================= "
echo " "


# Set the script to crash if any un set variables are used (we put this after
# the options processsing since some command line arguments may legitimately not
# exist).
set -o nounset

# Read out root install directory
#--------------------------------
MY_HOME_WD=`pwd`


# Choose build directory (for lib,include), relative to root
#------------------------------------------------------------
build_dir=$MY_HOME_WD/build

echo " "
echo " "
echo "I'm going to install the distribution (the lib and include directories)"
echo "in:"
echo " "
echo "    " $build_dir
echo " "
echo " "
if ! YesNoRead "Is this OK?" "y"; then
    OptionPrompt "Specify build directory [e.g. /home/joe_user/build] :"
    build_dir=`OptionRead`
fi

echo " "
echo "============================================================= "
echo " "
echo "Build directory is: " 
echo " " 
echo "     " $build_dir
echo " " 
echo "--> The include directory will be in: "
echo " " 
echo "    " $build_dir"/include" 
echo " " 
echo "--> The lib directory will be in: "
echo " " 
echo "    " $build_dir"/lib" 
echo " "
echo "etc.       " 
echo " "
echo "============================================================= "
echo " "


# Choose configure options file
#------------------------------

# Ask if the initial options are OK
echo " "
echo "Configure options are: "
cat "config/configure_options/current"
echo 
if YesNoRead "Is this OK?" "y"; then
    accept_configure_options="true"
else
    accept_configure_options="false"
fi

# Continue asking if the options are OK until approved
while [[ $accept_configure_options != "true" ]]; do

    # Get list of options files
    configure_option_files="$(find config/configure_options -type f)"

    echo " "
    echo "======================================================================"
    echo 
    echo "Choose an alternative configuration file "
    # Loop over files and display a menu
    count=0
    for file in $configure_option_files
    do
        #Increase the counter
        count=`expr $count + 1`
        echo $count ": " $(basename $file)
    done

    echo
    echo "Enter the Desired configuration file [1-"$count"]"
    echo "Enter 0 to specify the options on the command line"

    # Read in the Desired File and validate it
    file_number=$(OptionRead)
    if (( $file_number >= $count )) || (( $file_number < 0 )); then
        # Error and go to start of loop
        echo "File number out of range, trying again." 1>&2
        continue
    fi


    # If options are to be read from the command line then store the
    # options in the file config/configure_options/current
    if [[ "$file_number" == "0" ]]; then
        echo 
        echo "Enter options"
        configure_options=`OptionRead`  
        echo $configure_options > "config/configure_options/current"

    # Otherwise copy the desired options file to config/configure_options/current
    else 
        # Use cut to extract the nth entry in the list
        configure_options_file="$(echo $configure_option_files | cut -d \  -f $file_number)"

        # Copy to current
        cp -f "$configure_options_file" "config/configure_options/current"
    fi

    # Check that the options are in the correct order
    configure_options_are_ok="$(CheckOptions config/configure_options/current)"
    if test "$configure_options_are_ok" != ""; then

        echo " " 1>&2
        echo "===============================================================" 1>&2
        echo "Error message from autogen.sh:" 1>&2
        echo " "  1>&2
        echo $configure_options_are_ok 1>&2
        echo " "  1>&2
        echo "===============================================================" 1>&2
        
        # Fail, go back to start of while loop
        continue
    fi

    # Ask if these options are OK
    echo " "
    echo "Configure options are: "
    cat "config/configure_options/current"
    echo 
    if YesNoRead "Is this OK?" "y"; then
        accept_configure_options="true"
    else
        accept_configure_options="false"
    fi

done


#====================================================================
# Start actual build process
#====================================================================

./non_interactive_autogen.sh "$@"

# echo "calling autogen as"
# echo "./non_interactive_autogen.sh" "$@"

echo " "
echo "autogen.sh has finished! If you can't spot any error messages" 
echo "above this, oomph-lib should now be ready to use... " 
echo " " 
echo "If you encounter any problems, please study the installation" 
echo "instructions and the FAQ before contacting the developers. " 
echo " "
echo "To run self tests use \"make check -k\" or ./bin/parallel_self_test.py"
