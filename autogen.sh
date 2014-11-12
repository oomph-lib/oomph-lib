#! /bin/bash

# TODO:    automatic hypre/trilinos pull source + install?

set -o errexit
set -o nounset

# Default values for arguments
oomph_root=""
build_dir=""
make_options=""
extra_configure_options=""
configure_options_file="config/configure_options/current"
generate_config_files="false"


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

# Convert some things to absolute paths
build_dir="$(AbsPath $build_dir)"
configure_options_file="$(AbsPath $configure_options_file)"

# Now switch to the oomph lib directory
cd "$oomph_root"

# Force generate config files if this is a brand new build.
if [ ! -e configure ]; then
    echo "No ./configure found, assuming this is a new build and regenerating everything."
    generate_config_files="true"
fi

# Print information about command line options selected
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


# David Shepherd's automake compatability fix
# ============================================================================
if $generate_config_files == "true"; then

    # This is an awful hack but I can't find any other way to handle it :(

    # If we have automake version more recent than 1.13 then the default is to
    # use the new parallel self test harness which doesn't work with
    # parallel_self_tests.py (and doesn't actually run tests in parallel
    # without a major rewrite of all Makefile.am s). So we need to disable it.
    # However the command to disable the new test harness was only introduced
    # in version 1.12 which is still very new! So it looks like the only way
    # around this for now is to modify configure.ac here if the automake
    # version is greater than 1.12.

    # Version comparison function from
    # http://stackoverflow.com/questions/3511006/how-to-compare-versions-of-some-products-in-unix-shell
    # We can't use the much simpler sort -V, or bash code because might not
    # have them on some systems...
    version_greater_equal_1_12_0()
    {
        # Put each value into a variable
        v1=$(echo $1 | cut -d "." --output-delimiter=" " -f 1)
        v2=$(echo $1 | cut -d "." --output-delimiter=" " -f 2)
        v3=$(echo $1 | cut -d "." --output-delimiter=" " -f 3)

        # The values to compare against
        c1=1
        c2=12
        c3=0

        # Test each value and echo the result
        if test $v1 -gt $c1; then
            echo "1"
        elif test $v1 -eq $c1; then
            if test $v2 -gt $c2; then
                echo "1"
            elif test $v2 -eq $c2; then
                if test $v3 -ge $c3; then
                    echo "1"
                else
                    echo "0"
                fi
            else
                echo "0"
            fi
        else
            echo "0"
        fi
    }


    # Now we need to get the automake version (hopefully they don't change the
    # formatting of the --version output! Don't think there's any other way to
    # get the version).
    automake_version=$(automake --version | head -1 | tr ' ' '\n' | tail -1)

    echo "Detected automake version $automake_version"
    if test $(version_greater_equal_1_12_0 $automake_version) -gt 0; then
        echo "I'm modifying configure.ac to use serial self tests in automake because you have a recent enough version of automake."
        # Enforce serial tests and require automake version 1.12 or above (just
        # in case someone does something really weird..)
        sed 's/^[ \t]*AM_INIT_AUTOMAKE[ \t]*$/AM_INIT_AUTOMAKE([1.12 serial-tests])/' \
            -i configure.ac
    else
        echo "Not setting serial tests option in configure.ac because your version"
        echo "of automake is old enough to use it by default (older than 1.12.0)."
    fi

fi



# Autodetection of folders in user_drivers
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


confdir="config/configure.ac_scripts"
cat "$confdir/core.dir_list" \
    "$confdir/doc.dir_list" \
    "$confdir/user_drivers.dir_list" \
    "$confdir/user_src.dir_list" \
    | sort \
    | sed -e 's|\(^.*$\)|\1/Makefile|' \
    | cat <(echo "AC_CONFIG_FILES([Makefile") - <(echo "])") \
    > "$confdir/new_makefile_list"

if ! diff -q "$confdir/new_makefile_list" "$confdir/makefile_list" > /dev/null 2>&1; 
then
    echo "New/removed user dirs detected and $confdir/new_makefile_list updated, configure will be rerun automatically by make."

    echo
    echo "User driver folders included are:"
    cat "config/configure.ac_scripts/user_drivers.dir_list"
    echo

    mv "$confdir/new_makefile_list" "$confdir/makefile_list"
fi



# Set up configure options
#--------------------------------------------------------


# Read the options from the files and convert them into a single one-line string
new_configure_options=$(ProcessOptionsFile "$configure_options_file")
old_configure_options=$(ProcessOptionsFile config/configure_options/current)

# If configure options have changed then run configure command
if [[ "$new_configure_options" != "$old_configure_options" || "$generate_config_files" == "true" ]]; then

    echo "Using configure options:"
    cat "$configure_options_file"
    echo

    # Check that the options are in the correct order
    configure_options_are_ok="$(CheckOptions config/configure_options/current)"
    if test "$configure_options_are_ok" != ""; then

        echo 1>&2
        echo "===============================================================" 1>&2
        echo "Error message from autogen.sh:" 1>&2
        echo  1>&2
        echo $configure_options_are_ok 1>&2
        echo  1>&2
        echo "===============================================================" 1>&2
        
        # Failed
        exit 4
    fi

    # Update current options
    cp "$configure_options_file" "config/configure_options/current"

    # Before running automake, build the auxillary files required in
    # /src/meshes ??ds this should probably go in a makefile!
    echo
    echo "Building Auxillary Files in /src/meshes"
    ./bin/build_mesh_makefile.sh .

    # Run all the autotools and just do the right things to generate
    # configure, Makefile.in and all the dependency relationships.
    autoreconf --install

    # Finally run configure itself to convert "Makefile.in"s into "Makefile"s
    echo
    echo "Running ./configure --prefix $build_dir $new_configure_options $extra_configure_options"
    echo
    /bin/sh -c "./configure --prefix $build_dir $new_configure_options $extra_configure_options"

    # Test that the mpi commands work (automatically passes if no variable
    # MPI_RUN_COMMAND in makefile). This needs to go after configure so that we
    # can use the generated Makefile to (robustly) get the run and compile
    # commands.
    ./bin/check_mpi_command.sh Makefile
fi


# make is smart enough to automatically run automake, configure etc. if
# they are needed for other reasons (e.g if we have added new dirs to one
# of the dir lists or modified a Makefile.am).


# Make all libraries
echo
echo "Running `make $make_options` in $PWD"
make $make_options
echo "done"


# Install the libraries (in build directory specified above)
# echo
# echo "running `make $make_options install` in $PWD"
make $make_options install
echo "done"
