#! /bin/bash

# Crash if any sub command crashes
set -o errexit

# Crash if any unset variables are used
set -o nounset

# Get the directory that autogen.sh is in (stolen from stackoverflow:
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
# ), this is the oomph-lib root directory. Doesn't follow symlinks to the
# script itself, should be robust for anything else. If you move autogen.sh
# this will need to change a little.
oomph_root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$oomph_root"

# Load helper functions
source bin/autogen_helpers.sh

# Initialise command line flags (single dash, single letter!)
# for David's non-interactive autogen.sh
non_interactive_flags=""

# Do you want to rebuild from scratch?
#-------------------------------------
# If so specify --rebuild as command line argument. Default is 
# to just do the normal ".configure, make, make install, make check" sequence.

# Bail out if more than two command line arguments
if (test $# -gt 2); then 
  echo " "
  echo "ERROR: Too many command line options!"
  echo "------"
  echo " "
  EchoUsageForInteractiveScript
fi   

# Process the command line options
raw_build=false
while (test $# -gt 0)
do
   case "$1" in
     # Set the rebuild flag
     --rebuild) 
      echo "             [Doing complete rebuild from scratch.]"
      raw_build=true
      non_interactive_flags=`echo $non_interactive_flags`" -r ";;
     # Set the jobs flag
     --jobs*)
      non_interactive_flags=`echo $non_interactive_flags`" -j "`echo $1| awk '{print substr($0,8)}'`;;
     # Anything else bail out     
      *)  
       echo " "
       echo "ERROR: Unrecognised command line option!"
       echo "------"
       echo " "
       EchoUsageForInteractiveScript;;
   esac
   shift
done

if (test "$raw_build" = "false"); then
   echo "                     [Doing normal build.]"
fi   

#====================================================================
# Start Q/A session
#====================================================================

echo " "
echo "============================================================= "
echo "        oomph-lib interactive installation script" 
echo "============================================================= "
echo " "

# Choose build directory (for lib,include)
build_dir="$oomph_root/build"

echo " "
echo " "
echo "I'm going to install the distribution (the lib and include directories)"
echo "in:"
echo " "
echo "    " $build_dir
echo " "
echo " "
if ! YesNoRead "Is this OK?" "y"; then
    printf "Specify build directory [e.g. /home/joe_user/build] :"
    build_dir=$(OptionRead)
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

# Wipe previous builds
#---------------------
if (test -d  $build_dir); then 
    echo " "
    echo "Note: Build directory " $build_dir " exists."
    echo " "
    #OptionPrompt "Do you want to wipe it [y/n -- default: n]"
    #reply=`OptionRead`
    #if test "$reply" = "y" -o "$reply" = "Y" ; then 
    if YesNoRead "Do you want to wipe it" "n"; then
       echo " "
       echo "Sorry to be over-cautious here, but we'd better double-check"
       echo "before we delete some of your precious files..."
       echo " "
       echo "The contents of " $build_dir " are:"
       echo " "
       ls -l  $build_dir
       echo " "
       if YesNoRead "Are you still sure you want to wipe it" "n"; then
          echo " "
          echo "Wiping it..."
          rm -f -r $build_dir
          echo "Done"
       fi
    fi
fi

echo "Self tests"
echo "=========="
echo "Following the installation of oomph-lib you can run a comprehensive set of"
echo "self tests with 'make check -k' or with './bin/parallel_self_test.py'. The"
echo "latter version tends to be much faster because it performs multiple"
echo "self-tests at the same time."
echo " " 
echo "Would you like to automatically run self tests "
if YesNoRead "(serially) if the build is successful?" "n"; then
    run_self_tests="true"
else
    run_self_tests="false"
fi

# Choose configure options file
#------------------------------

configure_options_file="config/configure_options/current"

# If the current file exists then ask if it is ok
if [[ -f $configure_options_file ]]; then

    echo " "
    echo " "
    echo "Configure options"
    echo "================="
    echo "Configure options are: "
    cat "$configure_options_file" | ProcessOptionsFile
    echo " "
    echo " " 
    if YesNoRead "Is this OK?" "y"; then
        accept_configure_options="true"
    else
        accept_configure_options="false"
    fi

    # Otherwise have to choose one interactively
else
    accept_configure_options="false"
fi

# Continue asking if the options are OK until approved
while [[ $accept_configure_options != "true" ]]; do

    # Get list of options files
    configure_option_files="$(find config/configure_options -type f | sort)"

    echo " "
    echo "======================================================================"
    echo 
    echo "Choose an alternative configuration file "
    # Loop over files and display a menu
    count=0
    for file in $configure_option_files
    do
        #Increase the counter
        count=$(expr $count + 1)
        echo $count ": " $(basename $file)
    done

    echo
    echo "Enter the Desired configuration file [1-"$count"]"
    echo "Enter 0 to specify the options on the command line"

    # Read in the Desired File and validate it
    file_number=$(OptionRead)
	if ! [[ $file_number =~ ^[0-9]+$ ]]; then
		echo "File number not readable as an integer" 1>&2
		continue
	fi

    if (( $file_number > $count )) || (( $file_number < 0 )); then
        # Error and go to start of loop
        echo "File number out of range, trying again." 1>&2
        continue
    fi

    # If options are to be read from the command line then store the
    # options in the file config/configure_options/current
    if [[ "$file_number" == "0" ]]; then
        echo 
        echo "Enter options"
        configure_options=$(OptionRead)
        echo $configure_options > "config/configure_options/new_options_file"
        configure_options_file="config/configure_options/new_options_file"

    # Otherwise copy the desired options file to config/configure_options/current
    else 
        # Use cut to extract the nth entry in the list
        configure_options_file="$(echo $configure_option_files | cut -d \  -f $file_number)"
    fi

    # Check that the options are in the correct order
    configure_options_are_ok="$(CheckOptions $configure_options_file)"
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
    cat "$configure_options_file" | ProcessOptionsFile
    echo 
    if YesNoRead "Is this OK?" "y"; then
        accept_configure_options="true"
    else
        accept_configure_options="false"
    fi
done

echo
echo
echo "Build behaviour after failure"
echo "============================="
echo "You can build/install/test the library using \"make -k\". This "
echo "keeps going after a failure and therefore builds/tests whatever "
echo "it can, but makes it harder to spot errors."
echo " "
if YesNoRead "Build with \"make -k\"?" "y"; then
      non_interactive_flags=`echo $non_interactive_flags`" -k "
fi
echo " "

echo
echo
echo "Minimise output from make process"
echo "================================="
echo "You can run ask for make to be run (more) silently, (via the "
echo "\"make -s\" option) which doesn't echo commands. "
echo " "
if YesNoRead "Build with \"make -s\"?" "n"; then
      non_interactive_flags=`echo $non_interactive_flags`" -s "
fi
echo " "

#====================================================================
# Start actual build process
#====================================================================

# Start the timer
start=`date +%s`

# Call non-interactive autogen
build_command="./non_interactive_autogen.sh -b $build_dir -c ${oomph_root}/$configure_options_file $non_interactive_flags"
echo "The interactive part of the build process is over."
echo "Running $build_command"

$build_command

# Run tests if requested
if test "$run_self_tests" == "true"; then
    self_test_command="make check -k"
    echo "Running self test command: $self_test_command"
    $self_test_command
fi

echo " "
echo "=============================================================== "
echo " "
echo "autogen.sh has finished! If you can't spot any error messages" 
echo "above this, oomph-lib should now be ready to use... " 
echo " " 
echo "If you encounter any problems, please study the installation" 
echo "instructions and the FAQ before contacting the developers. " 
echo
echo "To run self tests you can use 'make check -k' or './bin/parallel_self_test.py'"
echo " "
echo "=============================================================== "
echo " "

# End the timer
end=`date +%s`

# Calculate the time taken to run the non-interactive part of the script
time_taken=$((end-start))

# Output
echo "=============================================================== "
if test "$run_self_tests" == "true";
then
    echo "Time taken for setup and self-tests: $time_taken"
else    
    echo "Time taken for complete setup: $time_taken"
fi
echo "=============================================================== "
echo " "

