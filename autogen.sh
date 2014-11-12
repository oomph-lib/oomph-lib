#! /bin/bash

set -o errexit


#====================================================================
# A few helper functions
#====================================================================

# A little function 'borrowed' from the tecplot installation script...
OptionPrompt() 
{ 
    printf "%s " "$1" 
}


# Another little function 'borrowed' from the tecplot installation script...
OptionRead()
{
    read Opt
    echo $Opt
}


# Convert a string passed as argument to return true/false (ie a bool),
# return value 255 indicates error.
YNToBool()
{
    if [[ $1 == "y" || $1 == "Y" ]]; then
        return $(true) # 0 == true in bash
    elif [[ $1 == "n" || $1 == "N" ]]; then
        return $(false) # 1 is false in bash
    else
        echo "I don't understand \"$1\", y/n only please." 1>&2
        return 255
    fi
}


# Read a y/n answer from input, if the answer is not y/n then repeat until
# it is. Turns out this can be nasty when combined with set -o errexit: you
# can only call this inside an if statement otherwise bash thinks any false
# return value is an error!
YesNoRead()
{
    prompt="$1"
    default="$2"

    # Read an answer
    printf "%s y/n [default %s]: " "$1" "$2"
    read Opt

    # Convert to a bool
    bool=""
    if [[ $Opt == "" ]]; then
        YNToBool $default
        bool=$?
    else
        YNToBool $Opt
        bool=$?
    fi

    # If we didn't recognise it then try again
    if [[ $bool != 1 && $bool != 0 ]]; then
        YesNoRead $prompt $default
        return $?
    fi

    return $bool
}


# This little function takes the input, removes anything following a #
# deletes blanks lines and then replaces all newlines by spaces
ProcessOptionsFile()
{
    sed < $1 -e 's/#.*$//' -e '/^[ \t]*$/d' -e 's/\n/ /g'
}


#Check that the "--" options preceed other configure options.
CheckOptions()
{
    awk '
   BEGIN{encountered_first_minus_minus=0
         encountered_first_non_minus_minus_after_first_minus_minus=0}
   NF {# pattern NF ignores blank lines since it expands to 0 for empty lines!
   # Ignore any comments (first entry in row starts with "#")
   if (substr($1,1,1)!="#")
    { 
     # Does the first entry in the line start with "--"?
     if (substr($1,1,2)=="--"){encountered_first_minus_minus=1}

     # Have we encountered the first "--" entry already?
     if (encountered_first_minus_minus==1)
       {
        # Does the first entry in the line not start with "--"?
        if (substr($1,1,2)!="--")
         {
          encountered_first_non_minus_minus_after_first_minus_minus=1
         }
       }
     # Should if this is followed by another "--" entry!
     if ((encountered_first_minus_minus==1)&&
         (encountered_first_non_minus_minus_after_first_minus_minus==1))
      {
       if (substr($1,1,2)=="--")
        {
         ok=0
         print "ERROR: The entry\n\n" $0 "\n"
         print "is in the wrong place. All the \"--\" prefixed options should go first!\n\n"
        }
      }
    }
    }' `echo $1`
}


#This little function echo's the usage information
EchoUsage()
{
    echo "Usage: "
    echo "------ "
    echo " "
    echo "[without flags]: Normal \"./configure; make; make install; make check -k\" sequence."
    echo " "
    echo " --rebuild     : Complete re-configure, followed by normal build sequence."
    echo " "
    echo "--jobs[=N]     :  Run N make jobs simultaneously."
    echo "                  Useful for speeding up the build on multi-core processors." 
    exit
}




#====================================================================
# Start Q/A session
#====================================================================

echo " "
echo "============================================================= "
echo "              oomph-lib installation script" 
echo "============================================================= "
echo " "

# Do we want to rebuild from scratch?
#-------------------------------------

#Bail out if more than two command line arguments
if (test $# -gt 2); then 
    EchoUsage 
fi   

#Process the command line options
raw_build=false;
make_options=" ";
while (test $# -gt 0)
do
    case "$1" in
        #Set the rebuild flag
        --rebuild) 
            echo "             [Doing complete rebuild from scratch.]"
            raw_build=true;;
        #Set the jobs flag
        --jobs*)
            make_options="$1";;
        #Anything else bail out     
        *)  
            EchoUsage;;
    esac
    shift
done

if (test "$raw_build" = "false"); then
    echo "                     [Doing normal build.]"
fi   


# Set the script to crash if any un set variables are used (we put this after
# the options processsing since some command line arguments may legitimately not
# exist).
set -o nounset

# Read out root install directory
#--------------------------------
MY_HOME_WD=`pwd`



# If this is a rebuild: Check for helper scripts 
#-----------------------------------------------
if $raw_build; then

    SCRIPT_LIST=`echo config.guess config.sub depcomp install-sh ltmain.sh missing aclocal.m4 mkinstalldirs `
    SCRIPTS_EXIST="no"
    for script in $SCRIPT_LIST
    do
        if (test -e $script); then
            SCRIPTS_EXIST="yes"
        fi 
    done
    if test "$SCRIPTS_EXIST" = "yes" ; then 
        echo " "
        echo "You may wipe the symbolic links to the autoconf/automake helper scripts"
        echo " "
        for script in $SCRIPT_LIST
        do
            if (test -e $script); then
                echo "   " $script
            fi 
        done
        echo " "
        echo "[This is recommended if you have moved the sources to a different"
        echo " machine without packaging them up with make dist. The symbolic "
        echo " links tend to be machine-specific so it's best to force "
        echo " autoconf/automake to rebuild them on the new machine]."
        echo " "
        
        if YesNoRead "Do you want to wipe the helper scripts?" "n"; then
            echo " "
            echo "As a backup: Here are the old symbolic links:"
            echo " "
            for script in $SCRIPT_LIST
            do
                if (test -L $script); then
                    ls -L $script
                    ls -l $script > old_symbolic_links.txt
                fi
            done
            echo " "
            echo "We have stored this information in old_symbolic_links.txt"
            echo " "
            echo "Wiping them..."
            rm -f  $SCRIPT_LIST
            echo "Done"
        fi   
    else
        echo " "
        echo "[No autoconf/automake helper scripts to be wiped...]"
        echo " "
    fi
fi


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

# If "current" configure options file does not exist then copy in the
# default one:
if (test ! -f config/configure_options/current); then
    cp config/configure_options/default config/configure_options/current
fi

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


echo " "
echo "==================================================================="
echo " "
echo " "
echo "End of customisation -- the actual build process is about to start."
echo "This may take a while... No user intervention is required during"
echo "the build procedure, so go and take the dog for a walk..."
echo " " 
OptionPrompt "Hit enter to continue."
tmp=`OptionRead`





#====================================================================
# Start actual build process
#====================================================================


# Autodetect folders in user_drivers
#-----------------------------------

# Backup old file (use -f so it doesn't give an error if the file doesn't exist)
touch ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list
mv -f ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list.backup

# Get a list of locations of files named Makefile.am, modify a little and write to user_drivers.dir_list.
find ${MY_HOME_WD}/user_drivers -type f -name "Makefile.am" \
    | grep -v "^${MY_HOME_WD}/user_drivers/Makefile.am" \
    | sed 's:/Makefile.am$::' \
    | sed "s:^${MY_HOME_WD}/::" \
    > ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list

# The grep and sed commands above do the following: 1) Remove the line that
# corresponds to the Makefile.am in user_drivers itself. 2) Remove
# "/Makefile.am" from each line leaving only the directory (dirname doesn't work
# with pipes). 3) Remove the start of the path from each line leaving only the
# location relative to the oomph-lib root directory.

echo
echo "User driver folders included are:"
cat ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list
echo


# If we are doing a raw build or if ./configure does not yet exist then generate
# all config files needed.
#--------------------------------------------------------
if [ $raw_build -o ! -e ./configure ]; then
    $MY_HOME_WD/bin/regenerate_config_files.sh $MY_HOME_WD
fi

# Read the options from the file and convert them into a single one-line string
configure_options=`ProcessOptionsFile config/configure_options/current`

# Run configure command
echo " "
echo "Running ./configure --prefix $build_dir $configure_options"
echo " " 
/bin/sh -c "./configure --prefix $build_dir $configure_options"

echo " " 
echo " " 
echo "done running ./configure"
echo " " 
echo " " 

# Test that the mpi commands work (automatically passes if no variable
# MPI_RUN_COMMAND in makefile). This needs to go after configure so that we
# can use the generated Makefile to (robustly) get the run and compile
# commands.
$MY_HOME_WD/bin/check_mpi_command.sh $MY_HOME_WD/Makefile


# Make all libraries
echo " "
echo "Running make $make_options" 
make $make_options
echo "done"


# Install the libraries (in build directory specified above)
echo " "
echo "running make $make_options install"
make $make_options install
echo "done" 


echo " "
echo "autogen.sh has finished! If you can't spot any error messages" 
echo "above this, oomph-lib should now be ready to use... " 
echo " " 
echo "If you encounter any problems, please study the installation" 
echo "instructions and the FAQ before contacting the developers. " 
echo " "
echo "To run self tests use \"make check -k\" or ./bin/parallel_self_test.py"
