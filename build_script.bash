#! /bin/bash


# Crash if any sub command crashes
set -o errexit

# Crash if any unset variables are used
set -o nounset


# Default: verbose mode is off
verbose=off

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
    # Enable verbose mode
    --verbose) verbose=true ;; 
    -h|--help)
      echo "Usage: $0 [--verbose]"
      exit 0
      ;;
    *) 
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
  shift
done


# Get the directory that autogen.sh is in (stolen from stackoverflow:
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
# ), this is the oomph-lib root directory. Doesn't follow symlinks to the
# script itself, should be robust for anything else. If you move autogen.sh
# this will need to change a little.
oomph_root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$oomph_root"


# Customisation options
#----------------------

# Get all files in the directory
FILES=(config/cmake_build_options/*)

if [ ! -e config/cmake_build_options/customised_options.json ]; then
    echo " "
    echo "We're currently building with the default options. Would you like to customise these?"
    echo " "
    echo "Here are various customised option files. (Note that you can add your own by"
    echo "adding a suitable json file to the cmake_build_options directory):"
    echo " "
else
    echo "We're currently building with the customised options specified "
    echo "in the file"
    echo " "
    echo "       config/cmake_build_options/customised_options.json"
    echo " "
    echo "which contains: "
    echo " "
    cat config/cmake_build_options/customised_options.json
    echo " "
    echo "Note that these overwrite the relevant default settings specified"
    echo "in the file "
    echo " " 
    echo "       config/cmake_build_options/default_options.json"
    echo " " 
    echo "Here are various alternative customised option files. (Note that you can add your own by"
    echo "adding a suitable json file to the config/cmake_build_options directory):"
    echo " "
fi


choice_made=n
while [ $choice_made != "y" ]; do
    # List files with enumeration
    echo "Customised options in config/cmake_build_options:"
    for i in "${!FILES[@]}"; do
        echo "$((i+1)): $(basename "${FILES[$i]}")"
    done
    echo " "
    
    # Prompt user for selection
    read -p "Which one would you like?: " CHOICE
    
    # Validate input
    if [[ ! "$CHOICE" =~ ^[0-9]+$ ]] || (( CHOICE < 1 || CHOICE > ${#FILES[@]} )); then
        echo "Invalid choice. Exiting."
        exit 1
    fi
    
    # Show selected file
    SELECTED_FILE="${FILES[$((CHOICE-1))]}"
    echo "Here are the customised options from the selected file ($SELECTED_FILE)"
    echo " " 
    cat "$SELECTED_FILE"
    echo " "
    read -p "Is this OK? [y,n]" choice_made
    if [[  "$choice_made" != "y" && "$choice_made" != "n" ]]; then
        echo "Illegal response; doing it again."
        choice_made=n
    fi

done
cp "$SELECTED_FILE"  config/cmake_build_options/customised_options.json




#######################################################################################
# Get settings from json files: 
#######################################################################################

# Assign value/key from default options to shell variables of the same name. 
eval $(jq -r '
  to_entries[] | 
  .value | 
  to_entries | 
  map("\(.key)=" + @sh "\(.value)") | .[]
' config/cmake_build_options/default_options.json)


# hierher list options and loop until one is accepted

# Assign value/key from customised options to shell variables of the same name.
if [ -e config/cmake_build_options/customised_options.json ]; then
  eval $(jq -r '
   to_entries[] | 
   .value | 
   to_entries | 
   map("\(.key)=" + @sh "\(.value)") | .[]
'  config/cmake_build_options/customised_options.json)
else
  echo "No customised options specified in customised_options.json, so sticking with defaults"
fi


# Tell us what you got:
echo " " 
echo "Build settings:"
echo "---------------"
echo " " 
echo "General settings:"
echo "-----------------"
echo "- CMAKE_BUILD_TYPE                          : "$CMAKE_BUILD_TYPE
echo "- OOMPH_ENABLE_MPI                          : "$OOMPH_ENABLE_MPI
echo "- OOMPH_USE_OPENBLAS_FROM                   : "$OOMPH_USE_OPENBLAS_FROM " <--- hierher puneet: do we absolutely have to insist on open blas or can it be any blas?"
echo "- OOMPH_USE_BOOST_FROM                      : "$OOMPH_USE_BOOST_FROM
echo "- OOMPH_ENABLE_SANITISER_SUPPORT            : "$OOMPH_ENABLE_SANITISER_SUPPORT


echo " " 
echo "Third party library settings:"
echo "------------------------------"
if [ $BUILD_THIRD_PARTY_LIBRARIES == "ON" ]; then
    echo "- BUILD_THIRD_PARTY_LIBRARIES               : "$BUILD_THIRD_PARTY_LIBRARIES
    echo "- OOMPH_BUILD_OPENBLAS                      : "$OOMPH_BUILD_OPENBLAS
    echo "- OOMPH_BUILD_SUPERLU                       : "$OOMPH_BUILD_SUPERLU
    echo "- OOMPH_BUILD_SUPERLU_DIST                  : "$OOMPH_BUILD_SUPERLU_DIST
    echo "- OOMPH_BUILD_CGAL                          : "$OOMPH_BUILD_CGAL
    echo "- OOMPH_BUILD_MUMPS                         : "$OOMPH_BUILD_MUMPS
    echo "- OOMPH_BUILD_HYPRE                         : "$OOMPH_BUILD_HYPRE
    echo "- OOMPH_BUILD_TRILINOS                      : "$OOMPH_BUILD_TRILINOS
    echo "- OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS   : "$OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS
    echo "- OOMPH_THIRD_PARTY_INSTALL_DIR             : "$OOMPH_THIRD_PARTY_INSTALL_DIR
else
    echo " Not building third-party libraries!"
fi
echo " " 
echo "oomph-lib settings:"
echo "-------------------"
echo "- OOMPH_DONT_SILENCE_USELESS_WARNINGS           : " $OOMPH_DONT_SILENCE_USELESS_WARNINGS
echo "- OOMPH_MPI_NUM_PROC                            : " $OOMPH_MPI_NUM_PROC
echo "- OOMPH_ENABLE_PARANOID                         : " $OOMPH_ENABLE_PARANOID
echo "- OOMPH_ENABLE_RANGE_CHECKING                   : " $OOMPH_ENABLE_RANGE_CHECKING
echo "- OOMPH_SUPPRESS_TRIANGLE_LIB                   : " $OOMPH_SUPPRESS_TRIANGLE_LIB
echo "- OOMPH_SUPPRESS_TETGEN_LIB                     : " $OOMPH_SUPPRESS_TETGEN_LIB
echo "- OOMPH_USE_OPENBLAS_FROM                       : " $OOMPH_USE_OPENBLAS_FROM
echo "- OOMPH_USE_GKLIB_FROM                          : " $OOMPH_USE_GKLIB_FROM
echo "- OOMPH_USE_METIS_FROM                          : " $OOMPH_USE_METIS_FROM
echo "- OOMPH_USE_SUPERLU_FROM                        : " $OOMPH_USE_SUPERLU_FROM
echo "- OOMPH_USE_PARMETIS_FROM                       : " $OOMPH_USE_PARMETIS_FROM
echo "- OOMPH_USE_SUPERLU_DIST_FROM                   : " $OOMPH_USE_SUPERLU_DIST_FROM
echo "- OOMPH_USE_BOOST_FROM                          : " $OOMPH_USE_BOOST_FROM
echo "- OOMPH_USE_CGAL_FROM                           : " $OOMPH_USE_CGAL_FROM
echo "- OOMPH_USE_MUMPS_FROM                          : " $OOMPH_USE_MUMPS_FROM
echo "- OOMPH_USE_TRILINOS_FROM                       : " $OOMPH_USE_TRILINOS_FROM
echo "- OOMPH_USE_HYPRE_FROM                          : " $OOMPH_USE_HYPRE_FROM



# OpenBlas and Boost are absolutely required, so check that
# we're either going to build it via our own third party
# build machinery or that the libraries have been specified
echo " "
echo "Checking spec of essential libraries, openblas and boost:"
echo " " 
third_party_library_list="OPENBLAS BOOST"
for third_party_library in `echo $third_party_library_list`; do
    full_var="OOMPH_USE_"$third_party_library"_FROM"

    # We're planning to use the third-party library built by us
    if [ `echo ${!full_var} | grep "<oomph_third_party_install_dir_if_it_exists>" | wc -c` != 0 ]; then
        # ... we we'd better make sure that we actually build it
        if [ $BUILD_THIRD_PARTY_LIBRARIES == "ON" ]; then
            echo "Third party library "$third_party_library" will be built by us."
        else
            echo "ERROR: "${full_var}" = "${!full_var}" but third-party "
            echo "       libraries will not be built by us! Please change"
            echo "       the setting of BUILD_THIRD_PARTY_LIBRARIES in in the file"
            echo " "
            echo "          config/cmake_build_options/customised_options.json"
            echo " "
            exit
        fi
    # otherwise the specified library must exist
    else
        echo "Third party library "$third_party_library" will not be built by us, so has to exist."
        if [ -e ${!full_var} ]; then
            echo "OK; great. The third party library "${full_var}" = "${!full_var} " exists!"
        else
            echo "ERROR: The third party library "${!full_var} " does not exist!"
            echo "       Please correct the setting of "${full_var}" in the file"
            echo " "
            echo "          config/cmake_build_options/customised_options.json"
            echo " "
            echo "       or omit it and enable us to build the library for you by setting"
            echo "       BUILD_THIRD_PARTY_LIBRARIES=\"ON\"."
            echo " "
            exit
        fi
    fi

done



#######################################################################################
# Third Party build
#######################################################################################
if [ $BUILD_THIRD_PARTY_LIBRARIES == "ON" ]; then
    echo " "
    echo "========================================================================== "
    echo " " 
    echo "Starting third party library build:"
    echo "-----------------------------------"

    cmake_flags="-DOOMPH_ENABLE_MPI=$OOMPH_ENABLE_MPI"
    cmake_flags=$cmake_flags" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE"
# not used in third party libraries hierher puneet: true?    cmake_flags=$cmake_flags" -DOOMPH_ENABLE_SANITISER_SUPPORT=$OOMPH_ENABLE_SANITISER_SUPPORT"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_OPENBLAS=$OOMPH_BUILD_OPENBLAS"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_SUPERLU=$OOMPH_BUILD_SUPERLU"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_CGAL=$OOMPH_BUILD_CGAL"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_MUMPS=$OOMPH_BUILD_MUMPS"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_HYPRE=$OOMPH_BUILD_HYPRE"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_TRILINOS=$OOMPH_BUILD_TRILINOS"
    cmake_flags=$cmake_flags" -DOOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS=$OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS"

    if [ $OOMPH_ENABLE_MPI == "ON" ]; then
        cmake_flags=$cmake_flags" -DOOMPH_BUILD_SUPERLU_DIST=$OOMPH_BUILD_SUPERLU_DIST"
    else
        cmake_flags=$cmake_flags" -DOOMPH_BUILD_SUPERLU_DIST=OFF"
    fi

    if [ $OOMPH_THIRD_PARTY_INSTALL_DIR == "<oomph_root>/external_distributions/install/" ]; then
        echo "Installing third-party libraries into "$oomph_root"/external_distributions/install/"
    else
        # hierher check this
        cmake_flags=$cmake_flags" -DOOMPH_THIRD_PARTY_INSTALL_DIR=$OOMPH_THIRD_PARTY_INSTALL_DIR"
        echo "Installing third-party libraries into "$OOMPH_THIRD_PARTY_INSTALL_DIR
    fi
    echo " "
    echo "cmake flags for third-party library build: "  $cmake_flags
    echo " "
    echo "Let the build commence!"
    echo " " 
    cd external_distributions
    if [ "$verbose" == true ]; then 
        echo " "
        echo "Note: You can hide the full log by running the script without --verbose"
        echo " "
        cmake -G Ninja -B build $cmake_flags
        cmake --build build    
    else
        echo " "
        echo "Note I'm writing full log of the configure/build/install process to "
        echo " "
        echo "        oomph_lib_third_party_library_build.log"
        echo " " 
        echo "You can make it visible by running the script with --verbose"
        echo " "
        cmake -G Ninja -B build $cmake_flags > oomph_lib_third_party_library_build.log
        cmake --build build >> oomph_lib_third_party_library_build.log
    fi
else
    echo " "
    echo "Bypassing third party library build:"
    echo " "
fi


#######################################################################################
# Oomph-lib build
#######################################################################################
cd "$oomph_root"

echo " "
echo "========================================================================== "
echo " " 
echo "Starting oomph-lib build:"
echo "-------------------------"


# Default flags
#--------------
cmake_flags="-DOOMPH_ENABLE_MPI=$OOMPH_ENABLE_MPI"
cmake_flags=$cmake_flags" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE"
cmake_flags=$cmake_flags" -DOOMPH_ENABLE_SANITISER_SUPPORT=$OOMPH_ENABLE_SANITISER_SUPPORT"
cmake_flags=$cmake_flags" -DOOMPH_DONT_SILENCE_USELESS_WARNINGS=$OOMPH_DONT_SILENCE_USELESS_WARNINGS" 
cmake_flags=$cmake_flags" -DOOMPH_MPI_NUM_PROC=$OOMPH_MPI_NUM_PROC" 
cmake_flags=$cmake_flags" -DOOMPH_ENABLE_PARANOID=$OOMPH_ENABLE_PARANOID" 
cmake_flags=$cmake_flags" -DOOMPH_ENABLE_RANGE_CHECKING=$OOMPH_ENABLE_RANGE_CHECKING" 
cmake_flags=$cmake_flags" -DOOMPH_SUPPRESS_TRIANGLE_LIB=$OOMPH_SUPPRESS_TRIANGLE_LIB" 
cmake_flags=$cmake_flags" -DOOMPH_SUPPRESS_TETGEN_LIB=$OOMPH_SUPPRESS_TETGEN_LIB"


# Optional libraries:
#--------------------

# Which libraries do we have:
# -- Either the default ones, specified in the json file
#    for the default settings, and then built by us in the external
#    distribution directory. We skip them if they don't exist.
# -- Or any third party ones, explicitly specified in the customised json
#    file; again skipped and warned about if they don't exist.
third_party_library_list="OPENBLAS GKLIB METIS SUPERLU PARMETIS BOOST CGAL MUMPS TRILINOS HYPRE"

# Add SuperLU_dist only if we're in mpi mode
if [ $OOMPH_ENABLE_MPI == "ON" ]; then
    third_party_library_list=$third_party_library_list" SUPERLU_DIST "
fi

# Check if the libraries exist, replacing the "<oomph_root>" placeholder
# for where we've actually installed them.
for third_party_library in `echo $third_party_library_list`; do
    full_var="OOMPH_USE_"$third_party_library"_FROM"
    # Replace any placeholders for our own built libraries; this will do nothing if we've specified some other existing library
    lib_dir_name=`echo ${!full_var} | sed "s|<oomph_third_party_install_dir_if_it_exists>|$oomph_root/external_distributions/install|g"`
    if [ -e $lib_dir_name ]; then
        cmake_flags=$cmake_flags" -DOOMPH_USE_"$third_party_library"_FROM="$lib_dir_name
        echo "Yay! Library "$lib_dir_name" does exist; we're using it for oomph-lib build."

    else
        echo "Library "$lib_dir_name" doesn't exist."
        echo "Not using library "$third_party_library" in oomph-lib build."
    fi
done


echo " "
echo "cmake flags for oomph-lib build: "  $cmake_flags
echo " "
echo "Let the build commence!"
echo " "

if [ "$verbose" == true ]; then 
   echo " "
   echo "Note: You can hide the full log by running the script without --verbose"
   echo " "
   cmake -G Ninja -B build $cmake_flags
   cmake --build build 
   cmake --install build 
else
   echo " "
   echo "Note I'm writing full log of the configure/build/install to "
   echo " "
   echo "     oomph_lib_build.log"
   echo " " 
   echo "You can make it visible by running the script with --verbose"
   echo " "
   cmake -G Ninja -B build $cmake_flags > oomph_lib_build.log 2>&1
   cmake --build build  >> oomph_lib_build.log 2>&1
   cmake --install build >> oomph_lib_build.log 2>&1
fi






echo " "
echo "done; doing a quick serial test:"
echo " "
cd $oomph_root/demo_drivers/poisson/one_d_poisson
cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE
cmake --build build
cd build/
ctest


if [ $OOMPH_ENABLE_MPI == "ON" ]; then
    echo " "
    echo "...and a parallel one:"
    echo " "
    cd $oomph_root/demo_drivers/mpi/solvers
    cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE
    cmake --build build
    cd build/
    ctest
fi
