#! /bin/bash

echo " "
echo " "
echo "================================================================"
echo " "
echo " "
echo "NOTE: oomph-lib cgal build currently ignores oomph-lib build"
echo "      flags because cgal's build system is quite different"
echo "      from ours. Next time somebody gets bored, we should try"
echo "      to port CC, CXX, CFLAGS and CXXFLAGS across..."
echo " "
echo " "
echo "================================================================"
echo " "
echo " "

echo $*

CGAL_PERMANENT_INSTALL_DIRECTORY=""
if [ $# -eq 5 ]; then
    echo  "bla"
    CGAL_PERMANENT_INSTALL_DIRECTORY=`pwd`
elif [ $# -eq 6 ]; then
    CGAL_PERMANENT_INSTALL_DIRECTORY=$6
else
    echo "Need five or six args: GMP_INCLUDE GMP_LIB MPFR_INCLUDE MPFR_LIB BOOST_INCLUDE (CGAL_PERMANENT_INSTALLATION_DIR)"
    exit 1
fi


# update this if we move to a new version
tar_file=CGAL-4.11.tar.xz

gmp_include=$1                      # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/gmp/gmp_default_installation/include
gmp_lib=$2                          # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/gmp/gmp_default_installation/lib
if [ `uname` == "Darwin" ]
then
    gmp_actual_library=$2/libgmp.dylib     # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/gmp/gmp_default_installation/lib/libgmp.dylib
else
gmp_actual_library=$2/libgmp.so     # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/gmp/gmp_default_installation/lib/libgmp.so
fi
mpfr_include=$3                     # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/mpfr/mpfr_default_installation/include
mpfr_lib=$4                         # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/mpfr/mpfr_default_installation/lib
if [ `uname` == "Darwin" ]
then
    mpfr_actual_library=$4/libmpfr.dylib  # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/mpfr/mpfr_default_installation/lib/libmpfr.dylib
else
mpfr_actual_library=$4/libmpfr.so   # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/mpfr/mpfr_default_installation/lib/libmpfr.so
fi
boost_include=$5                    # /home/mheil/version_for_merging_in_louis_stuff/external_distributions/boost/boost_default_installation/include
                                    # boost_lib=/home/mheil/version_for_merging_in_louis_stuff/external_distributions/boost/boost_default_installation/lib



if [ ! -e $gmp_actual_library ]; then
    echo " "
    echo "ERROR! gmp library doesn't exist as: "
    echo "    "$gmp_actual_library 
    echo " "
    exit
else
    echo " "
    echo "gmp library exists as: "
    echo "    "$gmp_actual_library 
    echo " "
fi

if [ ! -e $mpfr_actual_library ]; then
    echo " "
    echo "ERROR! mpfr library doesn't exist as: "
    echo "    "$mpfr_actual_library 
    echo " "
    exit
else
    echo " "
    echo "mpfr library exists as: "
    echo "    "$mpfr_actual_library 
    echo " "
fi


install_dir=$CGAL_PERMANENT_INSTALL_DIRECTORY/"cgal_default_installation"
echo "install dir: " $install_dir
lib_dir=$install_dir/lib
echo "lib dir: " $lib_dir
include_dir=$install_dir/include
echo "include dir: " $include_dir


if [ -e $install_dir ]; then
    echo "cgal install dir already exists -- deleting it!"
    rm -rf $install_dir
else
    echo "cgal install dir doesn't exist yet; will be created during installation"
fi

dir=`basename $tar_file .tar.xz`
echo "dir: " $dir
if [ -e $dir ]; then
    if [ ! -d $dir ]; then
        echo "Directory " $dir " is not a directory"
        exit
    fi
fi
#ls -l
#echo "About to delete: "$dir
#rm -rf $dir
#ls -l
#pwd

# Get ready for new install
tar xf $tar_file

# Hack to comment out rounding check
dir_above_tar_file=`pwd`
offensive_file=`find . -name 'Interval_nt.h'`
offensive_dir=`dirname $offensive_file`
cd $offensive_dir
mv Interval_nt.h Interval_nt.h.orig
sed 's/CGAL_assertion_msg(-CGAL_IA_MUL/\/\/CGAL_assertion_msg(-CGAL_IA_MUL/g' Interval_nt.h.orig | sed 's/\"Wrong rounding:/\/\/\"Wrong rounding:/g' | sed 's/CGAL_assertion_msg(-CGAL_IA_DIV/\/\/CGAL_assertion_msg(-CGAL_IA_DIV/g' > Interval_nt.h
echo " " 
echo "oomph-lib has commented out rounding check"
echo " " 
cd $dir_above_tar_file

cd $dir
src_dir=`pwd`
mkdir $install_dir
cd $install_dir
echo "about to start build in: "`pwd`
build_opts="-DWITH_GMP=true \
            -DGMP_INCLUDE_DIR=$gmp_include \
            -DGMP_LIBRARIES=$gmp_actual_library \
            -DWITH_MPFR=true \
            -DMPFR_INCLUDE_DIR=$mpfr_include \
            -DMPFR_LIBRARIES=$mpfr_actual_library \
            -DBOOST_INCLUDEDIR=$boost_include \
            $src_dir"
echo $build_opts
cmake $build_opts
make

# for some reason the header files don't come across...
echo "I'm here: " 
pwd
echo " "
echo "src_dir: $src_dir"
echo "copying: $src_dir/include/CGAL"
echo "into   : $include_dir/CGAL"
echo " "
cp -r $src_dir/include/CGAL/* $include_dir/CGAL



