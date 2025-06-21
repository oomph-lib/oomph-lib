#! /bin/bash

# Move all external libraries and distributions to permanent location
if [ $# -ne 1 ]; then
    echo " "
    echo "Please specify the absolute path to the directory in which "
    echo "you want the external libraries and distributions to be"
    echo "installed (e.g. /home/mheil/local)."
    echo " " 
    exit 1
fi


install_dir=$1

echo "Trying to install the external distributions/libraries in "
echo " "
echo "     "$install_dir
echo " "

#Does it even exist?
if [ ! -e $install_dir ]; then
    echo "Installation directory"    
    echo " "
    echo "     "$install_dir
    echo " "
    echo "doesn't exist. Bailing out."
    exit 1
else
    echo "Installation directory"    
    echo " "
    echo "     "$install_dir
    echo " "
    echo " exists!"
fi


oomph_dir=`pwd`

# Make blas library
if [ ! -e external_src/oomph_blas ]; then
    echo " "
    echo "Can't find the external_src/oomph_blas directory."
    echo "Are you running this script from the oomph-lib home directory?"
    echo "Bailing out."
    echo " "
else
    if [ -e $install_dir/lib/blas ]; then
        echo " "
        echo "Blas lib directory:"
        echo " "
        echo "   "$install_dir/lib/blas
        echo " "
        echo "already exists."
    else
        echo " "
        echo "Creating blas lib directory:"
        echo " "
        echo "   "$install_dir/lib/blas
        echo " "
        mkdir -p $install_dir/lib/blas
    fi
    cd external_src/oomph_blas
    echo "Rebuilding object files (if necessary)"
    make
    echo "Turning them into library"
    ar rs blas.a *.o
    echo "Moving"
    mv blas.a $install_dir/lib/blas/
    echo "Done!"
fi
cd $oomph_dir

# Make lapack library
if [ ! -e external_src/oomph_flapack ]; then
    echo " "
    echo "Can't find the external_src/oomph_flapack directory."
    echo "Are you running this script from the oomph-lib home directory?"
    echo "Bailing out."
    echo " "
else
    if [ -e $install_dir/lib/lapack ]; then
        echo " "
        echo "Lapack lib directory:"
        echo " "
        echo "   "$install_dir/lib/lapack
        echo " "
        echo "already exists."
    else
        echo " "
        echo "Creating lapack lib directory:"
        echo " "
        echo "   "$install_dir/lib/lapack
        echo " "
        mkdir -p $install_dir/lib/lapack
    fi
    cd external_src/oomph_flapack
    echo "Rebuilding object files (if necessary)"
    make
    echo "Turning them into library"
    ar rs lapack.a *.o
    echo "Moving"
    mv lapack.a $install_dir/lib/lapack/
    echo "Done!"
fi
cd $oomph_dir


# Install directories that can be moved in their entirety
moveable_dir_list="external_distributions/boost/boost_default_installation external_distributions/hypre/hypre_default_installation external_distributions/cgal/cgal_default_installation external_distributions/trilinos/trilinos_default_installation external_distributions/gmp/gmp_default_installation external_distributions/mpfr/mpfr_default_installation external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation"
for moveable_dir in `echo $moveable_dir_list`; do
    if [ -e $moveable_dir ]; then
        target_dir=$install_dir"/"`basename $moveable_dir`
        if [ -e $target_dir ]; then
            echo " "
            echo "Target directory"
            echo " "
            echo "    "$target_dir
            echo " " 
            echo "already exists; bailing out!"
            echo " " 
            exit
        else
            echo " "
            echo "Target directory"
            echo " "
            echo "    "$target_dir
            echo " "
            echo "doesn't exist yet -- can move!"
        fi
        echo " "
        echo "Moving default installation"
        echo " "
        echo "     "$moveable_dir
        echo " " 
        echo "to permanent location "
        echo " " 
        echo "     "$target_dir
        echo " "
        cp -r $moveable_dir $target_dir
        echo " "
        echo "Done!"
        echo " " 
    else
        echo " "
        echo "Default installation"
        echo " "
        echo "     "$moveable_dir
        echo " " 
        echo "doesn't exist. Not moving it!"
        echo " "
    fi
done

echo " "
echo " "
echo "Success!"
echo " "
echo "Now add the following flags to your configure script: "
echo " "
echo " "
if [ -e $install_dir/mumps_and_scalapack_default_installation/lib ]; then echo "--with-blacs=$install_dir/mumps_and_scalapack_default_installation/lib"; fi
if [ -e $install_dir/mumps_and_scalapack_default_installation/lib/libscalapack.a ]; then echo "--with-scalapack=$install_dir/mumps_and_scalapack_default_installation/lib/libscalapack.a"; fi
if [ -e $install_dir/mumps_and_scalapack_default_installation/lib/libpord.a ]; then echo "--with-pord=$install_dir/mumps_and_scalapack_default_installation/lib/libpord.a"; fi
if [ -e $install_dir/mumps_and_scalapack_default_installation ]; then echo "--with-mumps=$install_dir/mumps_and_scalapack_default_installation"; fi
if [ -e $install_dir/hypre_default_installation ]; then echo "--with-hypre=$install_dir/hypre_default_installation"; fi
if [ -e $install_dir/trilinos_default_installation ]; then echo "--with-trilinos=$install_dir/trilinos_default_installation"; fi
if [ -e $install_dir/boost_default_installation ]; then echo "--with-boost=$install_dir/boost_default_installation"; fi
if [ -e $install_dir/gmp_default_installation ]; then echo "--with-gmp=$install_dir/gmp_default_installation"; fi
if [ -e $install_dir/mpfr_default_installation ]; then echo "--with-mpfr=$install_dir/mpfr_default_installation"; fi
if [ -e $install_dir/cgal_default_installation ]; then echo "--with-cgal=$install_dir/cgal_default_installation"; fi
if [ -e $install_dir/lib/blas/blas.a ]; then echo "--with-blas=$install_dir/lib/blas/blas.a"; fi
if [ -e $install_dir/lib/lapack/lapack.a ]; then echo "--with-lapack=$install_dir/lib/lapack/lapack.a"; fi
echo " " 


exit



--with-blacs=/home/mheil/local/mumps_and_scalapack_default_installation/lib
--with-scalapack=/home/mheil/local/mumps_and_scalapack_default_installation/lib/libscalapack.a
--with-pord=/home/mheil/local/mumps_and_scalapack_default_installation/lib/libpord.a 
--with-mumps=/home/mheil/local/mumps_and_scalapack_default_installation
--with-blas=/home/mheil/local/lib/blas/blas.a
--with-lapack=/home/mheil/local/lib/lapack/lapack.a
--with-hypre=/home/mheil/local/hypre_default_installation_mpi
--with-trilinos=/home/mheil/local/trilinos_default_installation_mpi
--with-boost=/home/mheil/local/boost_default_installation
--with-gmp=/home/mheil/local/gmp_default_installation
--with-mpfr=/home/mheil/local/mpfr_default_installation
--with-cgal=/home/mheil/local/cgal_default_installation


exit


# Get tar files for external distributions from oomph-lib webpage
# Needs internet access and assumes that it's run from the oomph-lib
# home directory.

echo " "
echo "I'm trying to get the tar files for external distributions "
echo "from oomph-lib webpage."
echo " "
echo "CHECK 1: Do we have interweb access?"
echo " "

# Do we have the interweb?
# from https://stackoverflow.com/questions/929368/how-to-test-an-internet-connection-with-bash
url=https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/
wget -q --spider $url
if [ $? -eq 0 ]; then
    echo "Yay!"
    echo " "
else
    echo "Nay! Bummer. Can't connect to "
    echo " "
    echo "    "$url
    echo " "
    echo "Make sure sure your computer is connected to the web."
    echo "Bailing out"
    echo " "
    exit
fi

echo " "
echo "CHECK 2: Is this script run from the oomph-lib home directory?"
echo " "

# Check if we can find the installation directories:
dir_list='external_distributions/trilinos external_distributions/hypre external_distributions/mumps_and_scalapack'
for dir in `echo $dir_list`; do
    
    if [ ! -e $dir ]; then
        echo "Something is wrong! I can't find the directory "
        echo " "
        echo "     "$dir
        echo " "
        echo "Are you sure you're running this script from the oomph-lib"
        echo "home directory? Bailing out."
        echo " " 
        exit
    else
        echo "Yay! Found "$dir
    fi

done
echo " "

home_dir=`pwd`


echo " "
echo "Getting trilinos tar file: "
echo " " 
cd external_distributions/trilinos
tar_file=trilinos-11.8.1-Source.tar.gz
if [ -e $tar_file ]; then
    echo "        tar file "$tar_file" already exists; not downloading it again!"
else
    wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$tar_file
    if [ -e $tar_file ]; then
        echo "Success!"
    else
        echo "That doesn't seem to have worked; at least the tar file "
        echo " "
        echo "       "$tar_file
        echo " "
        echo "isn't where it should be..."
        echo " " 
        exit
    fi
fi
cd $home_dir

echo " "
echo "Getting hypre tar file: "
echo " " 
cd external_distributions/hypre
tar_file=hypre-2.0.0.tar.gz
if [ -e $tar_file ]; then
    echo "        tar file "$tar_file" already exists; not downloading it again!"
else
    wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$tar_file
    if [ -e $tar_file ]; then
        echo "Success!"
    else
        echo "That doesn't seem to have worked; at least the tar file "
        echo " "
        echo "       "$tar_file
        echo " "
        echo "isn't where it should be..."
        echo " " 
        exit
    fi
fi
cd $home_dir


echo " "
echo "Getting mumps/scalapack tar files: "
echo " " 
cd external_distributions/hypre
tar_file_list="scalapack_installer.tgz MUMPS_4.10.0.tar.gz"
for tar_file in `echo $tar_file_list`; do
    if [ -e $tar_file ]; then
        echo "        tar file "$tar_file" already exists; not downloading it again!"
    else
        wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$tar_file
        if [ -e $tar_file ]; then
            echo "Success!"
        else
            echo "That doesn't seem to have worked; at least the tar file "
            echo " "
            echo "       "$tar_file
            echo " "
            echo "isn't where it should be..."
            echo " " 
            exit
        fi
    fi
done
cd $home_dir

echo " " 
echo "Done; if you now (re-)install oomph-lib, we will build the above external"
echo "distributions from source and make them available to oomph-lib."
echo "However, note that the mumps/scalapack build requires the"
echo "specification of the "
echo " "
echo "   --with-mpi-include-directory"
echo " "
echo "configure flag where the argument should be the directory containing"
echo "the mpi.h header file of your mpi installation. On my machine this"
echo "is /usr/lib/openmpi/include, so the full flag should be"
echo " "
echo "   --with-mpi-include-directory=/usr/lib/openmpi/include"
echo " "
echo "You can use the locate command "
echo " "
echo "     locate mpi.h"
echo " "
echo "to, guess what, locate that directory!"
echo " "
echo " " 
