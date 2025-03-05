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







#######################################################################################
# Get settings from json files: 
#######################################################################################

# Assign value/key from default options to shell variables of the same name. 
eval $(jq -r '
  to_entries[] | 
  .value | 
  to_entries | 
  map("\(.key)=" + @sh "\(.value)") | .[]
' cmake_build_options/default_options.json)


# hierher list options and loop until one is accepted

# Assign value/key from customised options to shell variables of the same name.
if [ -e cmake_build_options/customised_options.json ]; then
  eval $(jq -r '
   to_entries[] | 
   .value | 
   to_entries | 
   map("\(.key)=" + @sh "\(.value)") | .[]
'  cmake_build_options/customised_options.json)
else
  echo "No customised options specified in customised_options.json, so sticking with defaults"
fi


# Tell us what you got:
echo " " 
echo "Default settings:"
echo "-----------------"
echo " " 
echo "General settings:"
echo "-----------------"
echo "- OOMPH_ENABLE_MPI                          : "$OOMPH_ENABLE_MPI
echo "- CMAKE_BUILD_TYPE                          : "$CMAKE_BUILD_TYPE

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
    echo "- OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTING : "$OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTING
    echo "- OOMPH_THIRD_PARTY_INSTALL_DIR             : "$OOMPH_THIRD_PARTY_INSTALL_DIR
    echo "- OOMPH_USE_OPENBLAS_FROM                   : "$OOMPH_USE_OPENBLAS_FROM   " <-- hierher test this (this is needed by third parties if blas  isn't built, right?)"
    echo "- OOMPH_USE_BOOST_FROM                      : "$OOMPH_USE_BOOST_FROM      " <-- hierher test this (this is needed by third parties if boost isn't built, right?)"
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


#######################################################################################
# Third Party build
#######################################################################################
if [ $BUILD_THIRD_PARTY_LIBRARIES == "ONhierher" ]; then
    echo " "
    echo "========================================================================== "
    echo " " 
    echo "Starting third party library build:"
    echo "-----------------------------------"

    cmake_flags="-DOOMPH_ENABLE_MPI=$OOMPH_ENABLE_MPI"
    cmake_flags=$cmake_flags" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_OPENBLAS=$OOMPH_BUILD_OPENBLAS"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_SUPERLU=$OOMPH_BUILD_SUPERLU"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_CGAL=$OOMPH_BUILD_CGAL"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_MUMPS=$OOMPH_BUILD_MUMPS"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_HYPRE=$OOMPH_BUILD_HYPRE"
    cmake_flags=$cmake_flags" -DOOMPH_BUILD_TRILINOS=$OOMPH_BUILD_TRILINOS"
    cmake_flags=$cmake_flags" -DOOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTING=$OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTING"

    if [ $OOMPH_ENABLE_MPI == "ON" ]; then
        cmake_flags=$cmake_flags" -DOOMPH_BUILD_SUPERLU_DIST=$OOMPH_BUILD_SUPERLU_DIST"
    else
        cmake_flags=$cmake_flags" -DOOMPH_BUILD_SUPERLU_DIST=OFF"
    fi

    echo "Installing third party libraries in $OOMPH_THIRD_PARTY_INSTALL_DIR"
    if [ $OOMPH_THIRD_PARTY_INSTALL_DIR == "<oomph_root>/external_distributions/install/" ]; then
        echo "no need to set directory..."
    else
       cmake_flags=$cmake_flags" -DOOMPH_THIRD_PARTY_INSTALL_DIR=$OOMPH_THIRD_PARTY_INSTALL_DIR"
    fi
    echo "cmake flags: "  $cmake_flags
    
    cd external_distributions
    cmake -G Ninja -B build $cmake_flags
    cmake --build build

    # hierher test these
    #-DOOMPH_USE_OPENBLAS_FROM=$OOMPH_USE_OPENBLAS_FROM \ 
    #-DOOMPH_USE_BOOST_FROM=$OOMPH_USE_BOOST_FROM 
    
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
    #echo "no  ! : "${full_var}"
    #echo "yes ! : "${!full_var}"
    lib_dir_name=`echo ${!full_var} | sed "s|<oomph_root>|$oomph_root|g"`
    if [ -e $lib_dir_name ]; then
        cmake_flags=$cmake_flags" -DOOMPH_USE_"$third_party_library"_FROM="$lib_dir_name
    else
        echo "Library "$lib_dir_name" doesn't exist."
        echo "Not using library "$third_party_library" in oomph-lib build."
    fi
done

    
echo "cmake flags: "  $cmake_flags


# hierher
echo "hierher check install dir; the above probably won't quite work if we're specifying another install dir: "$OOMPH_THIRD_PARTY_INSTALL_DIR
    
exit
cmake_flags=$cmake_flags" -DOOMPH_USE_GKLIB_FROM=$OOMPH_USE_GKLIB_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_METIS_FROM=$OOMPH_USE_METIS_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_SUPERLU_FROM=$OOMPH_USE_SUPERLU_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_PARMETIS_FROM=$OOMPH_USE_PARMETIS_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_BOOST_FROM=$OOMPH_USE_BOOST_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_CGAL_FROM=$OOMPH_USE_CGAL_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_MUMPS_FROM=$OOMPH_USE_MUMPS_FROM" 
cmake_flags=$cmake_flags" -DOOMPH_USE_TRILINOS_FROM=$OOMPH_USE_TRILINOS_FROM"
cmake_flags=$cmake_flags" -DOOMPH_USE_HYPRE_FROM=$OOMPH_USE_HYPRE_FROM" 

if [ $OOMPH_ENABLE_MPI == "ON" ]; then
    cmake_flags=$cmake_flags" -DOOMPH_USE_SUPERLU_DIST_FROM=$OOMPH_USE_SUPERLU_DIST_FROM"
fi

# hierher what if we dont have the (other) third party libraries?





cmake -G Ninja -B build $cmake_flags
cmake --build build
cmake --install build



exit


#############################################################################

# General build flags
#--------------------

# Build directory
build_dir="build"

# Build type
build_type="Debug"

# Configure flags
general_cmake_config_flags="-G Ninja -B $build_dir -DCMAKE_BUILD_TYPE=\"$build_type\""


echo "General build flag: "
echo " "
echo "   "$general_cmake_config_flags
echo " "


# Third party library flags:
#---------------------------
build_third_party_libraries=1

third_party_libraries_cmake_config_flags=$general_cmake_config_flags" -DOOMPH_ENABLE_MPI=ON"

if [ $build_third_party_libraries -eq 1 ]; then
    echo "Third party library build flag: "
    echo " "
    echo "   "$third_party_libraries_cmake_config_flags
    echo " "
fi


##############################################################################



# Build external libraries
if [ $build_third_party_libraries -eq 1 ]; then
    echo " "
    echo "Building third party libraries"
    echo " "    
    cd external_distributions
    cmake $third_party_libraries_cmake_config_flags
    cmake --build $build_dir
    echo " "
    echo "Done building third party libraries"
    echo " "    
fi



