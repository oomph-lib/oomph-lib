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
oomphlib_branch="cmake-beta"
clone_oomphlib=0
clone_stand_alone_code=0

# What do we demonstrate?
demo_default_cmake=1
demo_default_script=1
demo_non_default_install_cmake=1
demo_non_default_install_script=1
demo_developer_naughty_cmake=1
demo_developer_proper_cmake=1
demo_developer_script=1

# Back to the roots!
playground_home_dir=`pwd`

# If installing libraries outside oomph-lib, this is where we put them:
tpl_install_dir=$playground_home_dir/tpl_installation
oomph_lib_install_dir=$playground_home_dir/oomph_lib_installation



#==========================================
# helper function to wipe build/install dirs
#==========================================
cleanup()
{
    where_i_was=`pwd`
    echo "Cleaning up in "$playground_home_dir
    cd $playground_home_dir/oomph-lib
    rm -rf build install external_distributions/build external_distributions/install
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    rm -rf build
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
    cleanup
    
    cd $playground_home_dir/oomph-lib
    
    
    #-----------------------------------------------
    # Third party
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
    # libraries (built in the  previous step; the command shown here
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
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG"
    cd build
    ninja
    ./one_d_poisson

else

    echo "Not doing default cmake demo"

fi






###########################################################
# Default installation; script
###########################################################
if [ $demo_default_script == 1 ]; then
    
    
    echo "Doing default script demo"
    cleanup
    
    cd $playground_home_dir/oomph-lib
    ./oomph_build.py 
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=`pwd`/install
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG"
    cd build
    ninja
    ./one_d_poisson

else

    echo "Not doing default cmake demo"

fi






###########################################################
# Installation in specified diretory; raw cmake
###########################################################
if [ $demo_non_default_install_cmake == 1 ]; then
        
    echo "Doing non-default install location cmake demo"
    cleanup
    
    # Wipe existing installations
    rm -rf $tpl_install_dir
    rm -rf $oomph_lib_install_dir

    cd $playground_home_dir/oomph-lib

    #-----------------------------------------------
    # Third party
    #-----------------------------------------------
    
    # Go to the third-party distributions directory
    cd external_distributions
    
    # Configure the third party library build
    cmake -G Ninja -B build -DCMAKE_INSTALL_PREFIX=$tpl_install_dir
    
    # Build and install them (this takes a while). Note that no separate
    # install step is required here.
    cmake --build build
    
    
    #-----------------------------------------------
    # oomph-lib
    #-----------------------------------------------
    cd $playground_home_dir/oomph-lib
    
    # Configure and generate the build system. -G specifies the build
    # system generator; -B specifies the build directory (here "build").
    # The final argument specifies the location of the third-party
    # libraries (built in the  previous step); the command shown here
    # assumes that the third party libraries were installed in the
    # default location external_distributions/install. Update this if
    # you installed them somewhere else).
    cmake -G Ninja -B build $(cat $tpl_install_dir/cmake_flags_for_oomph_lib.txt) -DCMAKE_INSTALL_PREFIX=$oomph_lib_install_dir
    
    # Build the oomph-lib libraries (i.e. compile the sources
    # and turn them into libraries); specify the directory
    # that was created at the configure stage (here "build").
    # This takes a while...
    cmake --build build
    
    # Install, again specifying the build directory created above
    # (here "build")
    cmake --install build
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$oomph_lib_install_dir
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG"
    cd build
    ninja
    ./one_d_poisson

else

    echo "Not doing non-default install location cmake demo"

fi




###########################################################
# Installation in specified diretory; script
###########################################################
if [ $demo_non_default_install_script == 1 ]; then
    
    
    echo "Doing non-default install location script demo"
    cleanup
    
    # Wipe existing installations
    rm -rf $tpl_install_dir
    rm -rf $oomph_lib_install_dir

    # Install tpl and oomph-lib using script
    cd $playground_home_dir/oomph-lib
    ./oomph_build.py --ext-CMAKE_INSTALL_PREFIX=$tpl_install_dir  --oomph-CMAKE_INSTALL_PREFIX=$oomph_lib_install_dir 
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$oomph_lib_install_dir
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir
    cd build
    ninja
    ctest
    
    
    #-----------------------------------------------
    # Stand-alone code
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG"
    cd build
    ninja
    ./one_d_poisson

else
    
    echo "Not doing non-default install location script demo"
    
fi



###########################################################
# Installation in specified diretory; raw cmake; developer
# settings. Brute forcing a compiler flag into everything!
###########################################################
if [ $demo_developer_naughty_cmake == 1 ]; then
        
    echo "Doing non-default install location cmake developer demo"
    cleanup
    
    # Wipe existing installations
    rm -rf $tpl_install_dir
    rm -rf $oomph_lib_install_dir

    cd $playground_home_dir/oomph-lib

    #-----------------------------------------------
    # Third party
    #-----------------------------------------------
    
    # Go to the third-party distributions directory
    cd external_distributions
    
    # Configure the third party library build
    # Need to specify mpi here too!
    cmake -G Ninja -B build -DCMAKE_INSTALL_PREFIX=$tpl_install_dir \
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
    # The final argument specifies the location of the third-party
    # libraries (built in the  previous step); the command shown here
    # assumes that the third party libraries were installed in the
    # default location external_distributions/install. Update this if
    # you installed them somewhere else).
    # MPI use comes across from third party libraries.
    # We're brute forcing compiler flags into the oomph-lib compilation
    # but doing it this way (with the -DCMAKE_CXX_FLAGS flag) they're not
    # passed through to any downstream codes (user drivers, say), so have
    # to be redefined there to avoid the potential for nasty crashes.
    cmake -G Ninja -B build $(cat $tpl_install_dir/cmake_flags_for_oomph_lib.txt) \
          -DCMAKE_INSTALL_PREFIX=$oomph_lib_install_dir \
          -DCMAKE_BUILD_TYPE=Debug \
          -DOOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON \
          -DOOMPH_ENABLE_PARANOID=ON \
          -DOOMPH_ENABLE_RANGE_CHECKING=ON \
          -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"

    # Build the oomph-lib libraries (i.e. compile the sources
    # and turn them into libraries); specify the directory
    # that was created at the configure stage (here "build").
    # This takes a while...
    cmake --build build
    
    # Install, again specifying the build directory created above
    # (here "build")
    cmake --install build
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$oomph_lib_install_dir
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    # Note that the additional c++ compiler flags have to be-defined here; they
    # don't come across from oomph-lib!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    cd build
    ninja
    ctest
    
    
    # ..and an mpi example
    cd $playground_home_dir/oomph-lib
    cd demo_drivers/mpi/solvers
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    cd build
    ninja
    ctest

    
    #-----------------------------------------------
    # Stand-alone code; Note that if we remove
    # the cxx flag here, the macro doesn't get
    # carried across and there's potential for
    # nasty seg faults because the library
    # and the driver code see different versions
    # of the header file!
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG -DOOMPH_SPECIAL_FLAG  -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    cd build
    ninja
    ./one_d_poisson

else

    echo "Not doing non-default install location cmake developer demo"

fi




###########################################################
# Installation in specified diretory; raw cmake; developer
# settings. Proper handling of compiler flag for oomph-lib
# gets passed through to downstream stand-alone driver codes
###########################################################
if [ $demo_developer_proper_cmake == 1 ]; then
        
    echo "Doing non-default install location cmake developer demo"
    cleanup
    
    # Wipe existing installations
    rm -rf $tpl_install_dir
    rm -rf $oomph_lib_install_dir

    cd $playground_home_dir/oomph-lib

    #-----------------------------------------------
    # Third party
    #-----------------------------------------------
    
    # Go to the third-party distributions directory
    cd external_distributions
    
    # Configure the third party library build.
    # Need to specify mpi here!
    cmake -G Ninja -B build -DCMAKE_INSTALL_PREFIX=$tpl_install_dir \
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
    # The final argument specifies the location of the third-party
    # libraries (built in the  previous step); the command shown here
    # assumes that the third party libraries were installed in the
    # default location external_distributions/install. Update this if
    # you installed them somewhere else)
    # MPI use comes across from third party libraries
    cmake -G Ninja -B build $(cat $tpl_install_dir/cmake_flags_for_oomph_lib.txt) \
          -DCMAKE_INSTALL_PREFIX=$oomph_lib_install_dir \
          -DCMAKE_BUILD_TYPE=Debug \
          -DOOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON \
          -DOOMPH_ENABLE_PARANOID=ON \
          -DOOMPH_ENABLE_RANGE_CHECKING=ON \
          -DOOMPH_EXTRA_COMPILE_DEFINES="OOMPH_SPECIAL_FLAG=1 OOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION=1"

    # Note that the final option is equivalent to the brute force alternative
    #
    #   -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    #
    # but this would only get the flags into oomph-lib itself and wouldn't
    # pass them to downstream to driver codes; potential for massive problems!

    
    # Build the oomph-lib libraries (i.e. compile the sources
    # and turn them into libraries); specify the directory
    # that was created at the configure stage (here "build").
    # This takes a while...
    cmake --build build
    
    # Install, again specifying the build directory created above
    # (here "build")
    cmake --install build
    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$oomph_lib_install_dir
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out. Note that we now need to specify where oomph-lib was installed;
    # It's not in the default oomph-lib location!
    cd demo_drivers/poisson/one_d_poisson
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir

    # Note that if the flags were brute forced above, they have to be redefined here
    # they're not coming across automatically, so we'd have to add
    #
    #   -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    #
    # again

    
    cd build
    ninja
    ctest
    
    
    # ..and an mpi example
    cd $playground_home_dir/oomph-lib
    cd demo_drivers/mpi/solvers
    rm -rf build
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir

    # Note that if the flags were brute forced above, they have to be redefined here
    # they're not coming across automatically, so we'd have to add
    #
    #   -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    #
    # again
    
    cd build
    ninja
    ctest

    
    #-----------------------------------------------
    # Stand-alone code; Note that if we remove
    # the cxx flag here, the macro doesn't get
    # carried across and there's potential for
    # nasty seg faults because the library
    # and the driver code see different versions
    # of the header file!
    #-----------------------------------------------
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG"
    
    # Note that if the flags were brute forced above, they have to be redefined here
    # they're not coming across automatically, so we'd have to add
    #
    #   -DCMAKE_CXX_FLAGS="-DOOMPH_SPECIAL_FLAG -DOOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION"
    #
    # again

    cd build
    ninja
    ./one_d_poisson

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
    cleanup
    
    # Wipe existing installations
    rm -rf $tpl_install_dir
    rm -rf $oomph_lib_install_dir

    # Install tpl and oomph-lib using script
    cd $playground_home_dir/oomph-lib
    ./oomph_build.py --ext-CMAKE_INSTALL_PREFIX=$tpl_install_dir \
                     --oomph-CMAKE_INSTALL_PREFIX=$oomph_lib_install_dir \
                     --OOMPH_ENABLE_MPI=ON \
                     --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON \
                     --oomph-OOMPH_EXTRA_COMPILE_DEFINES="OOMPH_SPECIAL_FLAG=1 OOMPH_3_5_BRICK_FOR_MESHING_ONLY_NO_INTEGRATION=1"


    
    # Where have we installed oomph-lib?
    oomph_lib_install_dir=$oomph_lib_install_dir
    
    # check if it's likely to be correct, i.e. does it contain the oomphlibConfig.cmake file?
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
    
    # Check it out. Note that we now need to specify where oomph-lib was installed;
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
    # Stand-alone code
    #-----------------------------------------------
    echo "stand-alone code:"
    cd $playground_home_dir/stand_alone_oomph-lib_user_code
    cmake -G Ninja -B build -Doomphlib_ROOT=$oomph_lib_install_dir \
          -DCMAKE_CXX_FLAGS="-DMY_FLAG"
    cd build
    ninja
    ./one_d_poisson

else
    
    echo "Not doing developer script demo"
    
fi
