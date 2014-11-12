#! /bin/bash

. bin/autogen_helpers.sh

set -o errexit
set -o nounset

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


# Temporary defaults
raw_build="false"
oomph_root=$(pwd)
build_dir="${oomph_root}/build"
make_options=""
configure_options_file="config/configure_options/current"



# If this is a rebuild then clean up the helper scripts
#-----------------------------------------------
if $raw_build != "true"; then

    SCRIPT_LIST="config.guess config.sub depcomp install-sh ltmain.sh missing aclocal.m4 mkinstalldirs"

    echo " "
    echo "Wiping the symbolic links to the autoconf/automake helper scripts"
    echo "[This is recommended if you have moved the sources to a different"
    echo " machine without packaging them up with make dist. The symbolic "
    echo " links tend to be machine-specific so it's best to force "
    echo " autoconf/automake to rebuild them on the new machine]."
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



# Autodetect folders in user_drivers
#-----------------------------------

# Backup old file (use -f so it doesn't give an error if the file doesn't exist)
touch ${oomph_root}/config/configure.ac_scripts/user_drivers.dir_list
mv -f ${oomph_root}/config/configure.ac_scripts/user_drivers.dir_list ${oomph_root}/config/configure.ac_scripts/user_drivers.dir_list.backup

# Get a list of locations of files named Makefile.am, modify a little and write to user_drivers.dir_list.
find ${oomph_root}/user_drivers -type f -name "Makefile.am" \
    | grep -v "^${oomph_root}/user_drivers/Makefile.am" \
    | sed 's:/Makefile.am$::' \
    | sed "s:^${oomph_root}/::" \
    > ${oomph_root}/config/configure.ac_scripts/user_drivers.dir_list

# The grep and sed commands above do the following: 1) Remove the line that
# corresponds to the Makefile.am in user_drivers itself. 2) Remove
# "/Makefile.am" from each line leaving only the directory (dirname doesn't work
# with pipes). 3) Remove the start of the path from each line leaving only the
# location relative to the oomph-lib root directory.

echo
echo "User driver folders included are:"
cat ${oomph_root}/config/configure.ac_scripts/user_drivers.dir_list
echo


# Set up configure options
#--------------------------------------------------------

# If we are doing a raw build or if ./configure does not yet exist then generate
# the config files needed.
if [ $raw_build -o ! -e ./configure ]; then
    $oomph_root/bin/regenerate_config_files.sh $oomph_root
fi

# If "current" configure options file does not exist then copy in the
# default one:
if [[ "$configure_options_file" != "config/configure_options/current" ]]; then
    cp "$configure_options_file" "config/configure_options/current"
elif test ! -f config/configure_options/current; then
    cp "config/configure_options/default" "config/configure_options/current"
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
    
    # Failed
    exit 4
fi

# Read the options from the file and convert them into a single one-line string
configure_options=$(ProcessOptionsFile config/configure_options/current)



# Go!
#--------------------------------------------------------

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
$oomph_root/bin/check_mpi_command.sh $oomph_root/Makefile


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
