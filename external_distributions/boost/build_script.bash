#! /bin/bash


if [ $# -ne 3 ]; then
 echo "Need three args: BOOST_LIB BOOST_INCLUDE BOOST_DIR"
 exit 1
fi

BOOST_LIB=$1 
BOOST_INCLUDE=$2 
BOOST_DIR=$3

echo " "
echo " "
echo "================================================================"
echo " "
echo " "
echo "NOTE: oomph-lib boost build currently ignores oomph-lib build"
echo "      flags because boost's build system is quite different"
echo "      from ours. Next time somebody gets bored, we should try"
echo "      to port CC, CXX, CFLAGS and CXXFLAGS across..."
echo " "
echo " "
echo "================================================================"
echo " "
echo " "

# update this if we move to a new version
tar_file=boost_1_65_1.tar.gz

# clean up previous installation
dir=`basename $tar_file .tar.gz`
echo "dir: " $dir
ls -l
#if [ -e $dir ]; then
#    if [ ! -d $dir ]; then
#        echo "Directory " $dir " is not a directory"
#        exit
#    fi
#fi
#rm -rf $dir

install_dir=$BOOST_DIR 
echo "install dir: " $install_dir
lib_dir=$BOOST_LIB 
echo "lib dir: " $lib_dir
include_dir=$BOOST_INCLUDE 
echo "include dir: " $include_dir

# Check if installed version of boost already exists
if [ -e $include_dir ]; then
   printf "Was about to build boost but it already appears to have been \n"
   printf "installed because the directory:\n\n"
   printf "      " $include_dir" \n\n"
   printf "already exists. I'm not installing boost again.\n\n\n"
   exit
else
    echo "boost install dir doesn't exist yet; will be created during installation"
fi

# Get ready for new install
tar xf $tar_file
cd $dir
./bootstrap.sh --prefix=$install_dir --libdir=$libdir --includedir=$include_dir --without-libraries=python
./b2 install link=static --toolset=darwin --prefix=$install_dir --libdir=$libdir --includedir=$include_dir



