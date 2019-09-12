#! /bin/bash


# update this if we move to a new version
tar_file=mpfr-3.1.6.tar.gz


MPFR_LIB=$1 
MPFR_INCLUDE=$2 
MPFR_DIR=$3
gmp_include=$4 
gmp_lib=$5



# Strip out oomph-lib cxx compiler flags
conf_flags=""
count=0
for word in "$@"; do
    echo $word
    let count=$count+1
    if [ $count -gt 5 ]; then
        if [ "$word" != "-DPARANOID" -a "$word" != "-DRANGE_CHECKING" ]; then
            conf_flags=$conf_flags" "$word
        fi
    fi
done
echo "conf flags: " $conf_flags
echo "gmp_include: " $gmp_include
echo "gmp_lib: " $gmp_lib

echo $conf_flags > .tmp_flag_converter.txt
source .tmp_flag_converter.txt
echo "Building with: "
echo "CC      : " $CC
echo "CFLAGS  : " $CFLAGS
rm .tmp_flag_converter.txt


install_dir=$MPFR_DIR 
echo "install dir: " $install_dir
include_dir=$install_dir/include

# Check if installed version of mpfr already exists
if [ -e $include_dir ]; then
   printf "Was about to build mpfr but it already appears to have been \n"
   printf "installed because the directory:\n\n"
   printf "      " $include_dir" \n\n"
   printf "already exists. I'm not installing mpfr again.\n\n\n"
   exit
else
    echo "mpfr install dir doesn't exist yet; will be created during installation"
fi

# Clean up previous install
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

# Get ready for new install
tar xf $tar_file
cd $dir
./configure --prefix=$install_dir --with-gmp-include=$gmp_include --with-gmp-lib=$gmp_lib; make; make install



