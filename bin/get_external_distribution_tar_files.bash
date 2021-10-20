

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
cd external_distributions/mumps_and_scalapack
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
echo "Following the build process, you can then use the script"
echo " " 
echo "    bin/move_external_libraries_and_distributions_to_permanent_location.bash"
echo " "
echo "to move the libraries to a permanent location (so they don't have"
echo "to be re-built every time!"
echo " "
