#! /bin/bash

# TODO: automatically remove config.cache if options change?
#       automatically detect if user_drivers has changed and regen config files?
#       automatic hypre/trilinos pull source + install?
#       move user driver detection to bin/regenerate_config_files.sh?
#       merge scripts again?

set -o errexit
set -o nounset

# Default values for arguments
generate_config_files="false"
oomph_root=""
build_dir=""
make_options=""
extra_configure_options=""
configure_options_file="config/configure_options/current"

while getopts ":hrd:c:b:j:sk" opt; do
  case $opt in
      h)
          echo "Options for autogen.sh:"
          echo
          EchoUsage
          exit 0
          ;;

      r)
          generate_config_files="true"
          echo "Doing a complete rebuild from scratch."
          ;;
      d)
          oomph_root="$OPTARG"
          ;;
      c)
          configure_options_file="$OPTARG"
          ;;
      b)
          build_dir="$OPTARG"
          ;;

      # flags for make
      j)
          job_option="--jobs $OPTARG"
          make_options="$make_options $job_option"
          echo "Added make option $job_option"
          ;;
      k)
          k_option="--keep-going"
          make_options="$make_options $k_option"
          echo "Added make option $k_option"
          ;;
      s)
          silent_option="--silent LIBTOOLFLAGS=--silent"
          make_options="$make_options $silent_option"
          extra_configure_options="$extra_configure_options -q"
          echo "Added make option $silent_option, configure option -q"
          ;;


      \?)
          echo "Invalid option: -$OPTARG" >&2
          exit 3
          ;;
  esac
done

# Default value for oomph_root if not set by args
if [[ $oomph_root == "" ]]; then
    oomph_root=$(pwd)
fi

# and for build dir 
if [[ $build_dir == "" ]]; then
    build_dir=${oomph_root}/build
fi

# Load helper functions
source "${oomph_root}/bin/autogen_helpers.sh"

# Convert to absolute paths
build_dir="$(AbsPath $build_dir)"
configure_options_file="$(AbsPath $configure_options_file)"

# Now move to the oomph lib directory
cd "$oomph_root"

# Print information about options
echo
echo "Using oomph lib in directory \"$PWD\""
echo "Using configure options file \"$configure_options_file\""
echo "Placing built files in \"$build_dir\""
if [[ $make_options != "" ]]; then
    echo "Using make options: \"$make_options\""
fi
if [[ $extra_configure_options != "" ]]; then
    echo "Using extra configure options \"$extra_configure_options\""
fi



# If we are regenerating config files then clean up the helper scripts
#-----------------------------------------------
if $generate_config_files == "true"; then

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
touch "config/configure.ac_scripts/user_drivers.dir_list"
mv -f "config/configure.ac_scripts/user_drivers.dir_list" "config/configure.ac_scripts/user_drivers.dir_list.backup"

# Get a list of locations of files named Makefile.am, modify a little and write to user_drivers.dir_list.
find "user_drivers" -type f -name "Makefile.am" \
    | grep -v "^user_drivers/Makefile.am" \
    | sed 's:/Makefile.am$::' \
    | sed "s:^./::" \
    > "config/configure.ac_scripts/user_drivers.dir_list"

# The grep and sed commands above do the following: 1) Remove the line that
# corresponds to the Makefile.am in user_drivers itself. 2) Remove
# "/Makefile.am" from each line leaving only the directory (dirname doesn't work
# with pipes). 3) Remove the start of the path from each line leaving only the
# location relative to the oomph-lib root directory.

echo
echo "User driver folders included are:"
cat "config/configure.ac_scripts/user_drivers.dir_list"
echo


# Set up configure options
#--------------------------------------------------------

# generate the config files if needed.
if [[ $generate_config_files == "true" ]] || [ ! -e configure ]; then
    ./bin/regenerate_config_files.sh "$PWD"
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

echo "Using configure options:"
cat "config/configure_options/current"
echo

# Read the options from the file and convert them into a single one-line string
configure_options=$(ProcessOptionsFile config/configure_options/current)



# Go!
#--------------------------------------------------------

# Run configure command
echo " "
echo "Running ./configure --prefix $build_dir $configure_options $extra_configure_options"
echo " "
/bin/sh -c "./configure --prefix $build_dir $configure_options $extra_configure_options"

echo " "
echo " "
echo "done running ./configure"
echo " "
echo " "

# Test that the mpi commands work (automatically passes if no variable
# MPI_RUN_COMMAND in makefile). This needs to go after configure so that we
# can use the generated Makefile to (robustly) get the run and compile
# commands.
./bin/check_mpi_command.sh Makefile


# Make all libraries
echo " "
echo "Running `make $make_options` in $PWD"
make $make_options
echo "done"


# Install the libraries (in build directory specified above)
echo " "
echo "running `make $make_options install` in $PWD"
make $make_options install
echo "done"
