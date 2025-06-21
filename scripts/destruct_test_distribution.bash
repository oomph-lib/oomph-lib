#! /bin/bash

#======================================================
# Destruct test distribution by running full self-tests
# with (nearly) all configure options. Some options
# are hard-coded (e.g. the names of tar files for the
# external distributions and the location of the mpi
# header files (needed for mumps build) but are checked.
# Also marked with "[update]".
#
# A few more (names of mpi compilers; base url for
# location of tar files for external distributions,
# currently (and for the foreseeable future!)
#
#   https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/
#
# should probably be handled the same way but aren't
# yet. Works fine on leylandii/vummath.
#======================================================
if [ $# -ne 1 ]; then
 echo "Usage: destruct_test_distribution.bash NAME_OF_TAR_FILE"
 exit 1
fi

# Remember where we are
current_dir=`pwd`


full_command="destruct_test/paranoia_0_mpi_0_external_dist_0/oomph-lib-*/bin/find_errors_and_warnings_in_build_log.bash destruct_test/*/*/test_build.log"
filtered_command="destruct_test/paranoia_0_mpi_0_external_dist_0/oomph-lib-*/bin/find_errors_and_warnings_in_build_log.bash destruct_test/*/*/test_build.log  | grep -v deliberately_broken | grep -v with_warning | grep -v 'code with warning' "
test_logs_location="destruct_test/*/*/test_build.log"

echo " "
echo "---------------------------------------------------------------"
echo " "
echo "Destruct test script is currently customised for running on "
echo "School of Maths linux system. Some variables need to be updated"
echo "on other machines. A certain amount of self-checking is done"
echo "to make sure the variables work. Check the \"Customised variables\""
echo "section in the script if anything doesn't work out." 
echo "Names of compilers are currently hard-coded (for gcc/open-mpi)."
echo " "
echo "Full destruct test will require about 215G of diskspace and take "
echo "a long time to run. Tests will run in newly created sub-directory "
echo "destruct_test (and it won't run if that directory already exists)."
echo "On-screen output is redirected into "
echo " "
echo "          $test_logs_location" # destruct_test/*/*/test_build.log"
echo " "
echo "and can be checked for errors using"
echo " "
echo "    $full_command" #"   destruct_test/paranoia_0_mpi_0_external_dist_0/oomph-lib-*/bin/find_errors_and_warnings_in_build_log.bash destruct_test/*/*/test_build.log "
echo " "
echo "This will obviously complain about the buggy codes that were "
echo "deliberately introduced to create warnings/errors (to test"
echo "that the destruct test machinery itself works properly...). To omit"
echo "these use "
echo " "
echo "    $filtered_command" #"   destruct_test/paranoia_0_mpi_0_external_dist_0/oomph-lib-*/bin/find_errors_and_warnings_in_build_log.bash destruct_test/*/*/test_build.log  | grep -v deliberately_broken | grep -v with_warning | grep -v 'code with warning' "
echo " "
echo "instead. This should not show ANY errors/warnings."
echo "  "
echo "Please make a note of this now and then hit any key to continue"
read -n 1

echo "...here we go..."


# Break to bypass actual destruct test
if [ 1 == 1 ]; then

##########################################################################
##########################################################################
##########################################################################

#======================================================
# Customised variables
#======================================================

#------------------------------------------------------
# Directory that contains mpi.h [update]
#------------------------------------------------------
directory_with_mpi_dot_h=/usr/include/openmpi-x86_64/   # ok for leylandii and vummath
if [ ! -e ${directory_with_mpi_dot_h}mpi.h ]
    then
    echo " "
    echo "ERROR"
    echo " "
    echo "Didn't find "$directory_with_mpi_dot_h"mpi.h"
    echo " "
    echo "which is needed for mumps build). Please update "
    echo "variable directory_with_mpi_dot_h in script".
    echo " "
    exit
else
    echo " "
    echo "Managed to find"
    echo " "
    echo "    "$directory_with_mpi_dot_h"mpi.h"
    echo 
fi

#------------------------------------------------------
# Name of trilinos tar file at 
# https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/
# [update]
#------------------------------------------------------
trilinos_tar_file=trilinos-11.8.1-Source.tar.gz

#------------------------------------------------------
# Name of hypre tar file at 
# https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/
# [update]
#------------------------------------------------------
hypre_tar_file=hypre-2.0.0.tar.gz

#------------------------------------------------------
# Name of mumps tar file at 
# https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/
# [update]
#------------------------------------------------------
mumps_tar_file=MUMPS_4.10.0.tar.gz

#------------------------------------------------------
# Name of scalapack tar file at 
# https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/
# [update]
#------------------------------------------------------
scalapack_tar_file=scalapack_installer.tgz


echo " " 
echo "CHECKING EXTERNAL DISTRIBUTIONS"
echo "==============================="
echo " " 

full_url=https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$trilinos_tar_file 
echo "Tar file: "
echo " "
echo "     "$full_url
echo " "
# This is from http://stackoverflow.com/questions/2924422/how-do-i-determine-if-a-web-page-exists-with-shell-scripting
if curl --output /dev/null --silent --head --fail "$full_url" 
then
    echo "...found!"
else
    echo "...not found!. Pleaes update variable trilinos_tar_file in script."
    exit
fi
echo " "

full_url=https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$hypre_tar_file 
echo "Tar file: "
echo " "
echo "     "$full_url
echo " "
# This is from http://stackoverflow.com/questions/2924422/how-do-i-determine-if-a-web-page-exists-with-shell-scripting
if curl --output /dev/null --silent --head --fail "$full_url" 
then
    echo "...found!"
else
    echo "...not found!. Pleaes update variable hypre_tar_file in script."
    exit
fi
echo " "

full_url=https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$mumps_tar_file 
echo "Tar file: "
echo " "
echo "     "$full_url
echo " "
# This is from http://stackoverflow.com/questions/2924422/how-do-i-determine-if-a-web-page-exists-with-shell-scripting
if curl --output /dev/null --silent --head --fail "$full_url" 
then
    echo "...found!"
else
    echo "...not found!. Pleaes update variable mumps_tar_file in script."
    exit
fi
echo " "

full_url=https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/$scalapack_tar_file 
echo "Tar file: "
echo " "
echo "     "$full_url
echo " "
# This is from http://stackoverflow.com/questions/2924422/how-do-i-determine-if-a-web-page-exists-with-shell-scripting
if curl --output /dev/null --silent --head --fail "$full_url" 
then
    echo "...found!"
else
    echo "...not found!. Pleaes update variable scalapack_tar_file in script."
    exit
fi
echo " "


##########################################################################
##########################################################################
##########################################################################


#------------------------------------------------------
# Tar file that contains the distribution
#------------------------------------------------------
tar_file=$1

#--------------------------------------------------
# Directory that contains the unpacked distribution
#--------------------------------------------------
unpacked_dist_dir=`basename $tar_file .tar.gz`


#------------------------------------------------------
# Directory for destruct tests
#------------------------------------------------------
if [ -e destruct_test   ] 
then
   echo "Please delete directory destruct_test and try again"
   exit
fi
mkdir destruct_test
cp $tar_file destruct_test
cd destruct_test


#------------------------------------------------------
# MPI
#------------------------------------------------------
for do_mpi in 0 1; do
    
    
#------------------------------------------------------
# Paranoia
#------------------------------------------------------
    for do_paranoia in 0 1 2; do 
        
        
        
#------------------------------------------------------
# external distributions
#------------------------------------------------------
        for do_ext_dist in 0 1; do
            
            
#------------------------------------------------------
# Create build script
#------------------------------------------------------
            echo "#! /bin/bash" > build_script.bash
            chmod a+x build_script.bash
            echo "tar xfz "$tar_file >> build_script.bash  
            echo "cd "$unpacked_dist_dir >> build_script.bash 

            if [ $do_ext_dist -eq 1 ]
                then

                # TRILINOS [update version number]
                #---------------------------------
                echo "cd external_distributions/trilinos/  " >> build_script.bash 
                echo " wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/"$trilinos_tar_file >> build_script.bash 
                echo "if [ ! -e $trilinos_tar_file ]; then  " >> build_script.bash 
                echo "   echo \"Download of trilinos tar file failed\"  " >> build_script.bash 
                echo "   exit  " >> build_script.bash 
                echo "else  " >> build_script.bash   
                echo "    echo \"trilinos tar file succesfully downloaded\" " >> build_script.bash 
                echo "fi" >> build_script.bash 
                echo "cd ../.." >> build_script.bash 
                
                # HYPRE [updateversion number]
                #-----------------------------
                echo "cd external_distributions/hypre" >> build_script.bash 
                echo "wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/"$hypre_tar_file >> build_script.bash 
                echo "if [ ! -e $hypre_tar_file ]; then" >> build_script.bash 
                echo "    echo \"Download of hypre tar file failed\" " >> build_script.bash 
                echo "    exit" >> build_script.bash 
                echo "else" >> build_script.bash 
                echo "    echo \"hypre tar file succesfully downloaded\" " >> build_script.bash 
                echo "fi" >> build_script.bash 
                echo "cd ../.." >> build_script.bash 
                
                # MUMPS/SCALAPACK [update version number]
                #----------------------------------------
                echo "cd external_distributions/mumps_and_scalapack" >> build_script.bash 
                echo "wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/"$scalapack_tar_file >> build_script.bash 
                echo "wget https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles/"$mumps_tar_file >> build_script.bash 
                echo "if [ ! -e $mumps_tar_file ]; then" >> build_script.bash 
                echo "    echo \"Download of mumps tar file failed\" " >> build_script.bash 
                echo "    exit" >> build_script.bash 
                echo "else" >> build_script.bash 
                echo "    echo \"mumps tar file succesfully downloaded\" " >> build_script.bash 
                echo "fi" >> build_script.bash 
                echo "if [ ! -e $scalapack_tar_file ]; then" >> build_script.bash 
                echo "   echo \"Download of scalapack installer tar file failed\" " >> build_script.bash 
                echo "   exit" >> build_script.bash 
                echo "else" >> build_script.bash 
                echo "     echo \"scalapack installer tar file succesfully downloaded\" " >> build_script.bash 
                echo "fi" >> build_script.bash 
                echo "cd ../.." >> build_script.bash 
            fi


            echo "rm config/configure_options/current" >> build_script.bash 
            echo "echo \"--enable-symbolic-links-for-headers\" >> config/configure_options/current" >> build_script.bash 
            echo "echo \"--enable-multiple_teuchos_libraries\" >> config/configure_options/current" >> build_script.bash 
            echo "echo \"--enable-suppress-doc\" >> config/configure_options/current" >> build_script.bash  # hierher suppressing doc

            if [ $do_mpi -eq 1 ]
            then
                echo "echo \"--enable-MPI\" >> config/configure_options/current" >> build_script.bash 
                echo "echo \"--with-mpi-self-tests=\\\"mpirun -np 2\\\"\" >> config/configure_options/current" >> build_script.bash 
                echo "echo \"--with-mpi-self-tests-variablenp=\\\"mpirun -np OOMPHNP\\\"\" >> config/configure_options/current" >> build_script.bash 
                if [ $do_ext_dist -eq 1 ]
                then
                    bla=
                    echo "echo \"--with-mpi-include-directory=\\\""${directory_with_mpi_dot_h}"\\\"\"  >> config/configure_options/current" >> build_script.bash 
                fi
            fi
            case "$do_paranoia" in
                
                "0")
                    echo "echo \"CXXFLAGS=\\\"-Wall -O3\\\"\" >> config/configure_options/current " >> build_script.bash ;;
                "1")
                    echo "echo \"CXXFLAGS=\\\"-Wall -g -DPARANOID\\\"\" >> config/configure_options/current " >> build_script.bash ;;
                "2")
                    echo "echo \"CXXFLAGS=\\\"-Wall -g -DPARANOID -DRANGE_CHECKING \\\"\" >> config/configure_options/current " >> build_script.bash ;;
            esac
            echo "echo \"CFLAGS=\\\"-O3\\\"\" >> config/configure_options/current" >> build_script.bash  
            echo "echo \"FFLAGS=\\\"-O3\\\"\" >> config/configure_options/current" >> build_script.bash 
            echo "echo \"FFLAGS_NO_OPT=\\\"-O0\\\"\" >> config/configure_options/current" >> build_script.bash  

            if [ $do_mpi -eq 1 ]
                then
                echo "echo \"CXX=mpic++\" >> config/configure_options/current" >> build_script.bash 
                echo "echo \"CC=mpicc\" >> config/configure_options/current" >> build_script.bash 
                echo "echo \"F77=mpif77\" >> config/configure_options/current" >> build_script.bash 
                echo "echo \"FC=mpif90\" >> config/configure_options/current" >> build_script.bash 
                echo "echo \"LD=mpif77\" >> config/configure_options/current" >> build_script.bash 
            fi
            # Build and run self tests
            echo "./non_interactive_autogen.sh &> test_build.log  " >> build_script.bash            
            echo "make -k check &>> test_build.log  " >> build_script.bash       
            
#------------------------------------------------------
# Make the test directory and run tests
#------------------------------------------------------
            test_dir="paranoia_"$do_paranoia"_mpi_"$do_mpi"_external_dist_"$do_ext_dist
            mkdir $test_dir
            cp $tar_file $test_dir
            mv build_script.bash $test_dir
            cd $test_dir
            echo "ABOUT TO RUN BUILD SCRIPT IN: "`pwd`
            ./build_script.bash &
            cd ..
            
        done
    done
done


#------------------------------------------------------
# Wait until all background jobs have finished
#------------------------------------------------------
echo "Waiting for destruct tests to finish..."
wait


# hierher bypass actual destruct test
fi

#------------------------------------------------------
# Diagnose and store result in file for attachment
#------------------------------------------------------
cd $current_dir
redirected_full_command="$full_command > destruct_test/full_destruct_test_results.txt"
echo "$redirected_full_command" > tmp_run.bash
chmod a+x tmp_run.bash
./tmp_run.bash
rm tmp_run.bash

redirected_filtered_command="$filtered_command > destruct_test/filtered_destruct_test_results.txt"
echo "$redirected_filtered_command" > tmp_run.bash
chmod a+x tmp_run.bash
./tmp_run.bash
rm tmp_run.bash


message="Destruct test done.\n\nOn-screen output is redirected into\n\n         $test_logs_location \n\nand can be checked for errors using\n\n   $full_command
\n\nThis will obviously complain about the buggy codes that were\ndeliberately introduced to create warnings/errors (to test\nthat the destruct test machinery itself works properly...).\nTo omit these use\n\n  $filtered_command \n\ninstead. This should not show ANY errors/warnings.\n\n"
echo -e $message | mail -s "destruct test done" -a destruct_test/full_destruct_test_results.txt -a destruct_test/filtered_destruct_test_results.txt $USER
echo -e $message
echo " "
echo "done"
echo " "
