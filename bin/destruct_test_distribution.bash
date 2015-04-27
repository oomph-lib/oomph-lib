#! /bin/bash

#======================================================
# Destruct test distribution by running full self-tests
# with (nearly) all configure options.
#======================================================


echo "LIKELY TO BE BROKEN AFTER DAVID/MATTHIAS'S REWRITE; NEED TO UPDATE"
echo "THE RESPONSES"
exit

#------------------------------------------------------
# Tar file that contains the distribution [Adjust]
#------------------------------------------------------
tar_file=oomph-lib-0.90.tar.gz

#-------------------------------------------------------------
# Directory that contains the unpacked distribution  [Adjust]
#-------------------------------------------------------------
unpacked_dist_dir=oomph-lib-0.90


#------------------------------------------------------
# Location of blas and lapack [Adjust]
#------------------------------------------------------
blas_lapack="--with-blas=/home/mheil/local/lib/blas/blas.a --with-lapack=/home/mheil/local/lib/lapack/lapack.a "

#------------------------------------------------------
# Directory for destruct tests
#------------------------------------------------------
dir=destruct_test
#rm -rf $dir 
if [ -e $dir   ] 
then
   echo "Please delete directory $dir and try again"
   exit
fi
mkdir $dir
cp $tar_file $dir
cd $dir

#------------------------------------------------------
# MPI
#------------------------------------------------------
for do_mpi in 0 1; do
    
    if [ $do_mpi -eq 1 ]
        then
        mpi="--enable-MPI --with-mpi-self-tests=\"mpirun -np 2\" " 
    else
        mpi=" " 
    fi
    
#------------------------------------------------------
# Paranoia
#------------------------------------------------------
    for do_paranoia in 0 1 2; do 
        
        case "$do_paranoia" in
            
            "0")
            paranoia=" -O6 ";;
            
            "1")
            paranoia=" -g -DPARANOID ";;
            
            "2")
            paranoia=" -g -DPARANOID -DRANGE_CHECKING " ;;
            
        esac
        
        
#------------------------------------------------------
# external distributions
#------------------------------------------------------
        for do_ext_dist in 0 1; do
                        
            if [ $do_ext_dist -eq 1 ]
                then
                if [ $do_mpi -eq 1 ]
                    then
                    ext_dist="--with-hypre=/home/mheil/local/hypre_default_installation_mpi --with-trilinos=/home/mheil/local/trilinos_default_installation_mpi "
                else
                    ext_dist="--with-hypre=/home/mheil/local/hypre_default_installation_serial --with-trilinos=/home/mheil/local/trilinos_default_installation_serial "
                fi
            else
                ext_dist=" " 
            fi
            
            
#------------------------------------------------------
# Now build the configure options
#------------------------------------------------------
            configure_options="--enable-symbolic-links-for-headers "
            configure_options=$configure_options`echo $blas_lapack`" "
            configure_options=$configure_options`echo $ext_dist`" "
            configure_options=$configure_options`echo $mpi`" "
            configure_options=$configure_options" CXXFLAGS=\"-Wall "`echo $paranoia`"\" CFLAGS=\"-O6\" FFLAGS=\"-O6\""
            if [ $do_mpi -eq 1 ]
                then
                configure_options=$configure_options" CXX=mpic++ CC=mpicc F77=mpif77 LD=mpif77"
            fi

            echo " " 
            echo "CONFIGURE OPTIONS: " $configure_options
            echo " " 

            
#------------------------------------------------------
# Create build script
#------------------------------------------------------
            echo "#! /bin/bash" > build_script.bash
            chmod a+x build_script.bash
            echo "tar xvfz "$tar_file >> build_script.bash  
            echo "cd "$unpacked_dist_dir >> build_script.bash            
            echo "./autogen.sh < ../replies.txt > test_build.log & " >> build_script.bash

            # Replies to queries in autogen.sh:
            echo "y" > replies.txt # confirm build dir
            echo "y" >> replies.txt # doc size
            echo "n" >> replies.txt # config options are not ok
            echo "0" >> replies.txt # enter config options via command line
            echo $configure_options >> replies.txt # here they come...
            echo "y" >> replies.txt # config options are now ok
            echo "y" >> replies.txt # run the self tests
            echo " " >> replies.txt # hit return to start the build process


            
            
#------------------------------------------------------
# Make the test directory and run tests
#------------------------------------------------------
            test_dir="paranoia_"$do_paranoia"_mpi_"$do_mpi"_external_dist_"$do_ext_dist
            mkdir $test_dir
            cp $tar_file $test_dir
            mv build_script.bash $test_dir
            mv replies.txt $test_dir
            cd $test_dir
            ./build_script.bash &
            cd ..
            
        done
    done
done

