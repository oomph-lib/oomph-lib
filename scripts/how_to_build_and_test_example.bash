#! /bin/bash

###############################################################
# Script to test/illustrate all possible (?) combinations
###############################################################


# Don't allow uninialised variables (cf. cd $doesn't_exit; rm -rf *}
set -o nounset

# echo commands; may be able to extra them for md file
# set -x

#==========================================
# Global settings: What do you want me to do?
#==========================================

# clone repos?
oomphlib_branch="development"
clone_oomphlib=1
clone_stand_alone_code=1

# What do we demonstrate?
demo_default_cmake=1
demo_default_script=1
demo_non_default_install_cmake=1
demo_non_default_install_script=1
demo_developer_cmake=1
demo_developer_script=1

# Back to the roots!
playground_home_dir=`pwd`

# If installing libraries outside oomph-lib, this is where we put them:
non_default_tpl_install_dir=$playground_home_dir/tpl_installation
non_default_oomph_lib_install_dir=$playground_home_dir/oomph_lib_installation



#==========================================
# helper function to wipe build/install dirs
#==========================================
cleanup()
{
    where_i_was=`pwd`
    echo "Cleaning up in "$playground_home_dir
    if [ -e $playground_home_dir/oomph-lib ]; then
        cd $playground_home_dir/oomph-lib
        rm -rf build install external_distributions/build external_distributions/install
    fi
    if [ -e $playground_home_dir/stand_alone_oomph-lib_user_code ]; then
        cd $playground_home_dir/stand_alone_oomph-lib_user_code
        if [ -e build ]; then
            rm -rf build
        fi
    fi
    cd $where_i_was
}




#==========================================
# Check out the two repos
#==========================================

# Clone oomph-lib?
cd $playground_home_dir
if [ $clone_oomphlib == 1 ]; then
    git clone https://github.com/oomph-lib/oomph-lib.git
    cd oomph-lib
    echo " "
    echo "NOTE: I'm switching oomph-lib to branch  " $oomphlib_branch
    echo " " 
    git checkout $oomphlib_branch
else
    echo " "
    echo "not checking out oomph-lib"
    echo " "
fi

# Clone stand-alone code?
cd $playground_home_dir
if [ $clone_stand_alone_code == 1 ]; then
    git clone https://github.com/oomph-lib/stand_alone_oomph-lib_user_code.git
else
    echo " "
    echo "not checking out stand alone code"
    echo " "
fi




#==========================================
# Start building
#==========================================


###########################################################
# Default installation; raw cmake
###########################################################
if [ $demo_default_cmake == 1 ]; then
    
    
    echo "Doing default cmake demo"
    cd $playground_home_dir
    cleanup
    
    cd $playground_home_dir/oomph-lib
    
    
    #-----------------------------------------------
    # Third party libraries
    #-----------------------------------------------
    
    # Go to the third-party distributions directory
    cd external_distributions
    
    # Configure the third party library build
    cmake -G Ninja -B build
    
    # Build and install them (this takes a while). Note that no separate
    # install step is required here. By default the libraries are
    # installed locally in the external_distributions/install directory
    cmake --build build
    
    
    #-----------------------------------------------
    # oomph-lib
    #-----------------------------------------------
    cd $playground_home_dir/oomph-lib
    
    # Configure and generate the build system. -G specifies the build
    # system generator; -B specifies the build directory (here "build").
    # The final argument specifies the location of the third-party
    # libraries (built in the previous step; the command shown here
    # assumes that the third party libraries were installed in the
    # default location external_distributions/install. Update this if
    # you installed them somewhere else).
    cmake -G Ninja -B build $(cat external_distributions/install/cmake_flags_for_oomph_lib.txt)
    
    # Build the oomph-lib libraries (i.e. compile the sources
    # and turn them into libraries); specify the directory
    # that was created at the configure stage (here "build").
    # This takes a while...
    cmake --build build
    
    # Install, again specifying the build directory created above
    # (here "build")
    cmake --install build
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=`pwd`/install
    
    # Check if this is likely to be correct, i.e. does this directory
    # contain the oomphlibConfig.cmake file?
    legal=`find $oomph_lib_install_dir -name "oomphlibConfig.cmake" | wc -l`
    if [ $legal == 1 ]; then
        echo "Found exactly one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed; looking good!"
    else
        echo "Found less or more than one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed. Bailing!"
        exit 1
    fi
    
    # Test it
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code: Pass MY_FLAG macro to C++
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DMY_FLAG=ON
    cd build
    ninja
    # Note: the code outputs which flags it was compiled with, so you
    # can check if the flags provided during the oomph-lib build were
    # passed through to the stand-alone driver code.
    ./one_d_poisson

else

    echo "Not doing default cmake demo"

fi






###########################################################
# Default installation; script
###########################################################
if [ $demo_default_script == 1 ]; then
    
    
    echo "Doing default script demo"
    cd $playground_home_dir
    cleanup
    
    cd $playground_home_dir/oomph-lib
    ./oomph_build.py 
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=`pwd`/install
    
    # Check if this is likely to be correct, i.e. does the directory
    # contain the oomphlibConfig.cmake file?
    legal=`find $oomph_lib_install_dir -name "oomphlibConfig.cmake" | wc -l`
    if [ $legal == 1 ]; then
        echo "Found exactly one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed; looking good!"
    else
        echo "Found less or more than one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed. Bailing!"
        exit 1
    fi
    
    # Test it
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code: Pass MY_FLAG macro to C++
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DMY_FLAG=ON
    cd build
    ninja
    # Note: the code outputs which flags it was compiled with, so you
    # can check if the flags provided during the oomph-lib build were
    # passed through to the stand-alone driver code.
    ./one_d_poisson

else

    echo "Not doing default cmake demo"

fi






###########################################################
# Installation in specified directory; raw cmake
###########################################################
if [ $demo_non_default_install_cmake == 1 ]; then
        
    echo "Doing non-default install location cmake demo"
    cd $playground_home_dir
    cleanup
    
    # Wipe existing installations
    rm -rf $non_default_tpl_install_dir
    rm -rf $non_default_oomph_lib_install_dir

    cd $playground_home_dir/oomph-lib

    #-----------------------------------------------
    # Third party libraries
    #-----------------------------------------------
    
    # Go to the third-party distributions directory
    cd external_distributions
    
    # Configure the third party library build; specify install directory
    # for the libraries
    cmake -G Ninja -B build -DCMAKE_INSTALL_PREFIX=$non_default_tpl_install_dir
    
    # Build and install them (this takes a while). Note that no separate
    # install step is required here.
    cmake --build build
    
    
    #-----------------------------------------------
    # oomph-lib
    #-----------------------------------------------
    cd $playground_home_dir/oomph-lib
    
    # Configure and generate the build system. -G specifies the build
    # system generator; -B specifies the build directory (here "build").
    # Get configuration details for third-party libraries from the
    # cmake_flags_for_oomph_lib.txt file in their install directory.
    # Also specify non-default install directory for oomph-lib.
    cmake -G Ninja -B build $(cat $non_default_tpl_install_dir/cmake_flags_for_oomph_lib.txt) -DCMAKE_INSTALL_PREFIX=$non_default_oomph_lib_install_dir
    
    # Build the oomph-lib libraries (i.e. compile the sources
    # and turn them into libraries); specify the directory
    # that was created at the configure stage (here "build").
    # This takes a while...
    cmake --build build
    
    # Install, again specifying the build directory created above
    # (here "build")
    cmake --install build
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$non_default_oomph_lib_install_dir
    
    # Check if this is likely to be correct, i.e. does this directory
    # contain the oomphlibConfig.cmake file?
    legal=`find $oomph_lib_install_dir -name "oomphlibConfig.cmake" | wc -l`
    if [ $legal == 1 ]; then
        echo "Found exactly one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed; looking good!"
    else
        echo "Found less or more than one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed. Bailing!"
        exit 1
    fi
    
    # Test it. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code: Pass MY_FLAG macro to C++
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DMY_FLAG=ON
    cd build
    ninja
    # Note: the code outputs which flags it was compiled with, so you
    # can check if the flags provided during the oomph-lib build were
    # passed through to the stand-alone driver code.
    ./one_d_poisson

else

    echo "Not doing non-default install location cmake demo"

fi




###########################################################
# Installation in specified directory; script
###########################################################
if [ $demo_non_default_install_script == 1 ]; then
    
    
    echo "Doing non-default install location script demo"
    cd $playground_home_dir
    cleanup
    
    # Wipe existing installations
    rm -rf $non_default_tpl_install_dir
    rm -rf $non_default_oomph_lib_install_dir

    # Install tpl and oomph-lib using script; specify install directories
    # for third-party (external) and oomph-lib libraries
    cd $playground_home_dir/oomph-lib
    ./oomph_build.py --ext-CMAKE_INSTALL_PREFIX=$non_default_tpl_install_dir  --oomph-CMAKE_INSTALL_PREFIX=$non_default_oomph_lib_install_dir 
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$non_default_oomph_lib_install_dir
    
    # Check if this is likely to be correct, i.e. does this directory
    # contain the oomphlibConfig.cmake file?
    legal=`find $oomph_lib_install_dir -name "oomphlibConfig.cmake" | wc -l`
    if [ $legal == 1 ]; then
        echo "Found exactly one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed; looking good!"
    else
        echo "Found less or more than one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed. Bailing!"
        exit 1
    fi
    
    # Test it. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code: Pass MY_FLAG macro to C++
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DMY_FLAG=ON
    cd build
    ninja
    # Note: the code outputs which flags it was compiled with, so you
    # can check if the flags provided during the oomph-lib build were
    # passed through to the stand-alone driver code.
    ./one_d_poisson

else
    
    echo "Not doing non-default install location script demo"
    
fi






###########################################################
# Installation in specified directory; raw cmake; developer
# settings. Proper handling of compiler flag for oomph-lib
# gets passed through to downstream stand-alone driver codes
###########################################################
if [ $demo_developer_cmake == 1 ]; then
        
    echo "Doing non-default install location cmake developer demo"
    cd $playground_home_dir
    cleanup
    
    # Wipe existing installations
    rm -rf $non_default_tpl_install_dir
    rm -rf $non_default_oomph_lib_install_dir

    cd $playground_home_dir/oomph-lib

    #-----------------------------------------------
    # Third party libraries
    #-----------------------------------------------
    
    # Go to the third-party distributions directory
    cd external_distributions
    
    # Configure the third party library build. Specify their install directory.
    # Need to specify mpi here!
    cmake -G Ninja -B build -DCMAKE_INSTALL_PREFIX=$non_default_tpl_install_dir \
          -DOOMPH_ENABLE_MPI=ON 
    
    # Build and install them (this takes a while). Note that no separate
    # install step is required here.
    cmake --build build
    
    
    #-----------------------------------------------
    # oomph-lib
    #-----------------------------------------------
    cd $playground_home_dir/oomph-lib
    
    # Configure and generate the build system. -G specifies the build
    # system generator; -B specifies the build directory (here "build").
    # Get configuration details for third-party libraries from the
    # cmake_flags_for_oomph_lib.txt file stored in their install directory.
    # MPI use comes across from third party libraries. We're also compiling
    # in debug mode, with paranoia, range checking and symbolic links for
    # the oomph-lib headers.
    cmake -G Ninja -B build $(cat $non_default_tpl_install_dir/cmake_flags_for_oomph_lib.txt) \
          -DCMAKE_INSTALL_PREFIX=$non_default_oomph_lib_install_dir \
          -DCMAKE_BUILD_TYPE=Debug \
          -DOOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON \
          -DOOMPH_ENABLE_PARANOID=ON \
          -DOOMPH_ENABLE_RANGE_CHECKING=ON \
          -DOOMPH_EXTRA_COMPILE_DEFINES="OOMPH_SPECIAL_FLAG OOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    
    # Build the oomph-lib libraries (i.e. compile the sources
    # and turn them into libraries); specify the directory
    # that was created at the configure stage (here "build").
    # This takes a while...
    cmake --build build
    
    # Install, again specifying the build directory created above
    # (here "build")
    cmake --install build
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$non_default_oomph_lib_install_dir
    
    # Check if this is likely to be correct, i.e. does this directory
    # contain the oomphlibConfig.cmake file?
    legal=`find $oomph_lib_install_dir -name "oomphlibConfig.cmake" | wc -l`
    if [ $legal == 1 ]; then
        echo "Found exactly one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed; looking good!"
    else
        echo "Found less or more than one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed. Bailing!"
        exit 1
    fi
    
    # Test it. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir
    cd build
    ninja
    ctest
    
    
    # ..and an mpi example
    cd $playground_home_dir/oomph-lib
    cd demo_drivers/mpi/solvers
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir    
    cd build
    ninja
    ctest

    
    #-----------------------------------------------
    # Stand-alone code: Pass MY_FLAG macro to C++.
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DMY_FLAG=ON
    cd build
    ninja
    # Note: the code outputs which flags it was compiled with, so you
    # can check if the flags provided during the oomph-lib build were
    # passed through to the stand-alone driver code.
    ./one_d_poisson
    echo " "
    echo "Note how this code saw the OOMPH_SPECIAL_FLAG and "
    echo "OOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION flags because"
    echo "they were specified correctly and transparently during the "
    echo "oomph-lib build!"
    echo " "
else

    echo "Not doing non-default install location cmake developer demo"

fi


###########################################################
# Installation in specified directory; script; developer
# settings. Proper handling of compiler flag for oomph-lib
# gets passed through to downstream stand-alone driver codes
###########################################################
if [ $demo_developer_script == 1 ]; then
    
    
    echo "Doing developer script demo"
    cd $playground_home_dir
    cleanup
    
    # Wipe existing installations
    rm -rf $non_default_tpl_install_dir
    rm -rf $non_default_oomph_lib_install_dir

    # Install tpl and oomph-lib using script; specify install directories
    # mpi usage, and use of symbolic links for the headers. The
    # --oomph-OOMPH_EXTRA_COMPILE_DEFINES is the proper way to specify
    # C++ macros for oomph-lib. 
    cd $playground_home_dir/oomph-lib
    ./oomph_build.py --ext-CMAKE_INSTALL_PREFIX=$non_default_tpl_install_dir \
                     --oomph-CMAKE_INSTALL_PREFIX=$non_default_oomph_lib_install_dir \
                     --OOMPH_ENABLE_MPI=ON \
                     --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON \
                     --oomph-OOMPH_EXTRA_COMPILE_DEFINES="OOMPH_SPECIAL_FLAG OOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"


    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$non_default_oomph_lib_install_dir
    
    # Check if this is likely to be correct, i.e. does this directory
    # contain the oomphlibConfig.cmake file?
    legal=`find $oomph_lib_install_dir -name "oomphlibConfig.cmake" | wc -l`
    if [ $legal == 1 ]; then
        echo "Found exactly one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed; looking good!"
    else
        echo "Found less or more than one instance of oomphlibConfig.cmake in the directory"
        echo "      " $oomph_lib_install_dir
        echo "where oomph-lib is supposed to be installed. Bailing!"
        exit 1
    fi
    
    # Test it. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    echo "serial self-test:"
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir
    cd build
    ninja
    ctest
    
    # ..and an mpi example
    echo "parallel self-test:"
    cd $playground_home_dir/oomph-lib
    cd demo_drivers/mpi/solvers
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir 
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code: Pass MY_FLAG macro to C++
    #-----------------------------------------------
    echo "stand-alone code:"
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DMY_FLAG=ON
    cd build
    ninja
    # Note: the code outputs which flags it was compiled with, so you
    # can check if the flags provided during the oomph-lib build were
    # passed through to the stand-alone driver code.
    ./one_d_poisson
    echo " "
    echo "Note how this code saw the OOMPH_SPECIAL_FLAG and "
    echo "OOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION flags because"
    echo "they were specified correctly and transparently during the "
    echo "oomph-lib build!"
    echo " "

else
    
    echo "Not doing developer script demo"
    
fi
