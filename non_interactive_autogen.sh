#! /bin/bash

# TODO:    automatic hypre/trilinos pull source + install?

set -o errexit
set -o nounset


# Get the directory that autogen.sh is in (stolen frome stackoverflow:
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
# ), this is the oomph-lib root directory. Doesn't follow symlinks to the
# script itself, should be robust for anything else. If you move autogen.sh
# this will need to change a little.
oomph_root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load helper functions
source "${oomph_root}/bin/autogen_helpers.sh"



# Handle command line input
# ============================================================


# Default values for arguments
build_dir=""
make_options=""
extra_configure_options=""
generate_config_files="false"
only_generate_config_files="false"
configure_options_file="config/configure_options/current"

# If "current" configure options does not exist then use "default".
if [[ ! -f "${oomph_root}/config/configure_options/current" ]]; then
    echo "No current configure options found, copying over the default options"
    cp "${oomph_root}/config/configure_options/default" "${oomph_root}/config/configure_options/current"
fi

run_self_tests=0

# Parse command line arguments
while getopts ":hrd:c:b:j:skoS" opt; do
    case $opt in

        S) 
            echo "Will run self tests (serially) at end of build"
            run_self_tests=1
            ;;

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
        c)
            configure_options_file="$OPTARG"
            ;;
        b)
            build_dir="$OPTARG"
            ;;

        o)
            only_generate_config_files="true"
            generate_config_files="true"
            echo "Regenerating the config files only"
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
            echo "Valid options are:"
            echo
            EchoUsage
            exit 3
            ;;
    esac
done

# and for build dir 
if [[ $build_dir == "" ]]; then
    build_dir=${oomph_root}/build
fi

# Convert some things to absolute paths
build_dir="$(AbsPath $build_dir)"
configure_options_file="$(AbsPath $configure_options_file)"

# Now switch to the oomph lib directory
cd "$oomph_root"

# Force generate config files if we are missing important files
if [ ! -e configure ]; then
    echo "No ./configure found, assuming this is a new build and regenerating everything."
    generate_config_files="true"
elif [ ! -e Makefile ]; then
    echo "No Makefile found in root, I assume you did a distclean so I'll regenerate everything."
    generate_config_files="true"
elif [ ! -e Makefile.in ]; then
    echo "No Makefile.in found in root, no idea how you managed that but I'll regenerate everything anyway"
    generate_config_files="true"
fi

if [ ! -e $configure_options_file ]; then
    echo "Configure options file $configure_options_file does not exist"
    exit 8
fi

# Print information about command line options selected
echo
echo "Building using the oomph-lib source in \"$PWD\""
echo "Using configure options file \"$configure_options_file\""
echo "Placing built files in \"$build_dir\""
if [[ $make_options != "" ]]; then
    echo "Using make options: \"$make_options\""
fi
if [[ $extra_configure_options != "" ]]; then
    echo "Using extra configure options \"$extra_configure_options\""
fi


# If this is a new build or a forced rebuild then we need to do some extra
# stuff.
if [[ $generate_config_files == "true" ]]; then

    # David Shepherd's automake compatability fix
    #=========================================================================

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
        #v1=$(echo $1 | cut -d "." --output-delimiter=" " -f 1)
        #v2=$(echo $1 | cut -d "." --output-delimiter=" " -f 2)
        #v3=$(echo $1 | cut -d "." --output-delimiter=" " -f 3)

        v1=$(echo $1 | cut -d "."  -f 1)
        v2=$(echo $1 | cut -d "."  -f 2)
        v3=$(echo $1 | cut -d "."  -f 3)

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

    # Create the file to contain the automake init command
    automake_init_command_file="config/configure.ac_scripts/automake_init_command_file"
    echo '# File generated by autogen.sh, DO NOT MODIFY' > $automake_init_command_file

    # Check automake version and pick the command to use
    echo "Detected automake version $automake_version"
    if test $(version_greater_equal_1_12_0 $automake_version) -gt 0; then
        echo "I'm modifying configure.ac to use serial self tests in automake because you have a recent enough version of automake."

        # Enforce serial tests and require automake version 1.12 or above (just
        # in case someone does something really weird..)
        echo 'AM_INIT_AUTOMAKE([1.12 serial-tests foreign])' >> $automake_init_command_file
    else
        echo "Not setting serial tests option in configure.ac because your version"
        echo "of automake is old enough to use it by default (older than 1.12.0)."

        echo 'AM_INIT_AUTOMAKE([foreign])' >> $automake_init_command_file
    fi


    # ??ds not really sure why this is here or what it does
    echo
    echo "Building Auxillary Files in /src/meshes"
    ./bin/build_mesh_makefile.sh .

fi



# Autodetection of Makefiles to generate
#============================================================================

confdir="config/configure.ac_scripts"

# Generate a sorted list of all the makefiles in the project, wrap it into
# an autoconfigure command and put it into a file.
makefile_list="$(find . -path './external_distributions' -prune \
                    -o -type f -name 'Makefile.am' -print \
                | sed -e 's:Makefile\.am:Makefile:' -e 's:^./::' \
                | sort)"

# A bit more explanation of the above command:
# First we find all Makefile.ams in the project except those in
# ./external_distributions (because that is where hypre/trilinos/etc live
# and we don't want to build those directly).

# Then we remove the .am using sed to get a list of Makefiles to create. We
# also remove the "./" here because the autotools don't like it.

# Finally we sort the output so that the order of the resulting list is
# deterministic.


# Create the file containing the list of Makefiles
cat > "$confdir/new_makefile_list" <<EOF
# GENERATED FILE, DO NOT MODIFY.
AC_CONFIG_FILES([
$makefile_list
external_distributions/Makefile
external_distributions/hypre/Makefile
external_distributions/trilinos/Makefile
external_distributions/gmp/Makefile
external_distributions/mpfr/Makefile
external_distributions/boost/Makefile
external_distributions/cgal/Makefile
external_distributions/mumps_and_scalapack/Makefile
])
EOF
# In case you haven't seen it before: this writes the lines between <<EOF
# and EOF into the file $confdir/new_makefile_list. Variables are
# substituted as normal.

# If we found some new dirs then write it into the list file that is
# included in configure.ac and tell the user. The fact that we have
# modified a file included in configure.ac will cause make to rerun
# autoconf and configure.
touch "$confdir/makefile_list"
if ! diff -q "$confdir/new_makefile_list" "$confdir/makefile_list" > /dev/null 2>&1; 
then
    echo "New/removed directories detected and $confdir/makefile_list updated,"
    echo "./configure will be rerun automatically by make."
    mv "$confdir/new_makefile_list" "$confdir/makefile_list"
fi


# If this is a new build or a forced rebuild then we need to explicitly run
# all the autotools magic now.
if $generate_config_files == "true"; then
    # Run all the autotools and just do the right things to generate
    # configure, Makefile.in and all the dependency relationships.
    autoreconf --install --force
fi


if [[ $only_generate_config_files == "true" ]]; then
    # Done, exit with success code
    exit 0
fi


# Set up configure options
#============================================================================
# The folder of the configure options
configure_options_dirname=`dirname "$configure_options_file"`;

# The folder of the configure options
configure_options_basename=`basename "$configure_options_file"`;

# Temporary file to store flags
temporary_configure_options_file="$configure_options_dirname/.$configure_options_basename.$RANDOM";

# If the filename is taken then loop until we create one that isn't
while [ -f $temporary_configure_options_file ]
do
    # Make a new filename
    temporary_configure_options_file="$configure_options_folder/.$configure_options_basename.$RANDOM";
done

# Make a temporary copy
cp "$configure_options_file" "$temporary_configure_options_file"

# The additional flags needed to link against shared libraries. 
# PM: When building shared third-party libraries, we need to ensure that we
# can access their code irregardless of the load address. Such code is
# referred to as position independent code (PIC). To make sure the code is
# compiled as such we need to provide the -fPIC and -DPIC flags to our
# compiler options but we shouldn't expect the user to add them themselves
# Update the config file (if necessary) to include the -fPIC flag, and
# don't forget to include the -DPIC linker flag which ensures that
# we include code in #ifdef PIC /*...*/ #endif statements.
extra_flags_for_linking_against_shared_libraries="-fPIC -DPIC"

# Loop over the flags to include
for i_flag in $extra_flags_for_linking_against_shared_libraries
do
    # Insert the i-th flag in the set of flags to include
    InsertExtraFlags "$temporary_configure_options_file" "$i_flag"
done

# Read the options from the files and convert them into a single one-line string
new_configure_options=$(ProcessOptionsFile < "$temporary_configure_options_file")
old_configure_options=$(ProcessOptionsFile < config/configure_options/current)

# If configure options have changed then we need to reconfigure
if [[ "$new_configure_options" != "$old_configure_options" ||
	  "$generate_config_files" == "true" ]]; then

    # Check that the options are in the correct order
    configure_options_are_ok="$(CheckOptions $temporary_configure_options_file)"
    
    # Kill the temporary file
    rm -f "$temporary_configure_options_file"
    
    # Check if the options are okay
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
        
    # Slight problem here: if we change the options and add a new
    # driver at the same time then configure will end up being re-run twice.
    # Don't think there's anything we can do about it

    echo "Using configure options:"
    echo "$new_configure_options"

    # Update current options, unless the files are the same
    if [[ "$(AbsPath $configure_options_file)" != "$(AbsPath config/configure_options/current)" ]];
    then
        cp "$configure_options_file" "config/configure_options/current"
    fi

    # Finally run configure itself to convert "Makefile.in"s into "Makefile"s
    echo
    echo "Running ./configure --prefix $build_dir $new_configure_options $extra_configure_options " 
    echo
    /bin/sh -c "./configure --prefix $build_dir $new_configure_options $extra_configure_options " 

    # Test that the mpi commands work with these configure options
    # (automatically passes if no variable MPI_RUN_COMMAND in makefile).
    # This needs to go after configure so that we can use the generated
    # Makefile to (robustly) get the run and compile commands.
    set +e
    ./bin/check_mpi_command.sh Makefile
    set -e
fi

# If the temporary file still exists
if [ -f "$temporary_configure_options_file" ]
then  
    # Kill it
    rm -f "$temporary_configure_options_file"
fi

# make is smart enough to automatically rerun automake, configure etc. if
# they are needed for other reasons (e.g if we have added new dirs to one
# of the dir lists or modified a Makefile.am).



# Build!
# ============================================================


# Make all libraries
echo " "
echo "Running make with make_options:"
echo " "
echo $make_options
echo " " 

make $make_options

echo "done make from within noninteractive autogen.sh"


# Install the libraries (in build directory specified above)
echo " "
echo "Running make install with make_options"
echo " "
echo $make_options
echo " " 

make $make_options install

echo "done make install from within noninteractive autogen.sh"


