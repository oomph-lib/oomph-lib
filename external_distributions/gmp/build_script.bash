#! /bin/bash

#set -o xtrace

# update this if we move to a new version
tar_file=gmp-6.1.2.tar.xz

GMP_DIR=$1
CXX=$2
CC=$3

echo $conf_flags > .tmp_flag_converter.txt
source .tmp_flag_converter.txt
echo "Building with: "
echo "CXX     : " $CXX
echo "CC      : " $CC
rm .tmp_flag_converter.txt


# clean up previous installation
dir=`basename $tar_file .tar.xz`

# Get ready for new install
install_dir=$GMP_DIR 
echo "install dir: " $install_dir
if [ -e $install_dir ]; then
    echo "gmp install dir already exists -- not doing anything!"
    exit
else
    echo "gmp install dir doesn't exist yet; will be created during installation"
fi

tar xf $tar_file
cd $dir

CC=$CC CXX=$CXX ./configure --prefix=$install_dir
make
make install

