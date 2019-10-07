#! /bin/bash

#set -o xtrace

# update this if we move to a new version
tar_file=gmp-6.1.2.tar.xz

GMP_LIB=$1 
GMP_INCLUDE=$2 
GMP_DIR=$3

# Strip out oomph-lib cxx compiler flags
conf_flags=""
count=0
for word in "$@"; do
    let count=$count+1
    if [ $count -gt 3 ]; then
        if [ "$word" != "-DPARANOID" -a "$word" != "-DRANGE_CHECKING" ]; then
            conf_flags=$conf_flags" "$word
        fi
    fi
done

echo $conf_flags > .tmp_flag_converter.txt
source .tmp_flag_converter.txt
echo "Building with: "
echo "CXX     : " $CXX
echo "CXXFLAGS: " $CXXFLAGS
echo "CC      : " $CC
echo "CFLAGS  : " $CFLAGS
rm .tmp_flag_converter.txt


# clean up previous installation
dir=`basename $tar_file .tar.xz`
#echo "dir: " $dir
#if [ -e $dir ]; then
#    if [ ! -d $dir ]; then
#        echo "Directory " $dir " is not a directory"
#        exit
#    fi
#fi
#rm -rf $dir

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
echo "------------"
echo "configure  flags: " `echo $conf_flags`
echo "------------"
junk=`echo $conf_flags | sed "s/'//g"`
export CXX=$CXX
export CXXFLAGS="-g -Wall"
export CFLAGS="-g -Wall"
export CC=$CC

./configure --prefix=$install_dir
make
make install

