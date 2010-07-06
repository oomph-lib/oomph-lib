#! /bin/bash

#-------------------------------------------------------
# Script to strip down oomph-lib disribution to bare
# essentials to allow test-driving of installation
# process.
#-------------------------------------------------------
home_dir=`pwd`


#====================================================================
# A few helper functions
#====================================================================
# A little function 'borrowed' from the tecplot installation script...
OptionPrompt() 
{ 
 printf "%s " "$1" 
}

# Another little function 'borrowed' from the tecplot installation script...
OptionRead()
{
 read Opt
 if test "$Opt" = "" ; then
  Opt=$1
 fi
 echo $Opt
}
echo " " 
echo "======================================================== " 
echo " " 
echo "WARNING: This script strips out most of the distribution "
echo "         and retains only the bare minimum (the build"
echo "         machinery, src/generic and src/poisson) required"
echo "         to test-drive the installation procedure."
echo "         Files that are not required will be DELETED!"
echo " "
echo "======================================================== " 
echo " " 
echo " "
echo "We are currently in: " $home_dir
echo " " 
echo " "
echo " Do you want to strip down this distribution [y/n -- default: n]"
reply=`OptionRead`
if test "$reply" != "y" -a "$reply" != "Y" ; then 
    echo "not doing it"
    exit 0
else
    echo "doing it"
    exit 0
fi

# Delete command
del="echo Deleting:  `ls -d`"

# Generic top level stuff
#========================

# Update makefiles
cp $home_dir/config/stripped_down_files/Makefile.am.top_level Makefile.am
cp $home_dir/config/stripped_down_files/Makefile.am.bin bin/Makefile.am

# Get rid of user src and drivers
$del user_src
$del user_drivers

# Kill private directory
$del private

# Get rid of entire doc directory
$del doc




# Self tests: Only keep analyse stage
#====================================
cd self_test
mv analyse_self_tests ../tmp_analyse_self_tests
$del *
mv ../tmp_analyse_self_tests analyse_self_test
cp $home_dir/config/stripped_down_files/Makefile.am.self_test Makefile.am





# Demo drivers
#=============
cd $home_dir
cd demo_drivers


# Kill everything apart from Poisson and mpi
#-------------------------------------------
mv poisson ../tmp_poisson
mv mpi ../tmp_mpi
$del *
mv ../tmp_poisson poisson
mv ../tmp_mpi mpi
cp $home_dir/config/stripped_down_files/Makefile.am.demo_drivers Makefile.am


# Within Poisson get rid of everything apart from adaptive flux bc
#-----------------------------------------------------------------
cd poisson
mv two_d_poisson_flux_bc_adapt ../tmp_two_d_poisson_flux_bc_adapt
$del *
mv ../tmp_two_d_poisson_flux_bc_adapt two_d_poisson_flux_bc_adapt 
cp $home_dir/config/stripped_down_files/Makefile.am.demo_drivers.poisson Makefile.am
cd ..


# Within mpi get rid of everything apart from adaptive flux bc
#-------------------------------------------------------------
cd mpi
mv distribution ../tmp_distribution
$del *
mv ../tmp_distribution distribution
cp $home_dir/config/stripped_down_files/Makefile.am.demo_drivers.mpi Makefile.am

cd distribution
mv two_d_poisson_flux_bc_adapt ../tmp_two_d_poisson_flux_bc_adapt
$del * 
mv ../tmp_two_d_poisson_flux_bc_adapt two_d_poisson_flux_bc_adapt
cp $home_dir/config/stripped_down_files/Makefile.am.demo_drivers.mpi.distribution Makefile.am




# Sources
#========
cd $home_dir/src


# Keep only Poisson, generic and meshes
#--------------------------------------
mv generic ../tmp_generic
mv poisson ../tmp_poisson
mv meshes ../tmp_meshes
$del *
mv ../tmp_generic generic
mv ../tmp_poisson poisson
mv ../tmp_meshes  meshes
cp $home_dir/config/stripped_down_files/Makefile.am.src Makefile.am


# Keep only rectangular quad meshes
#----------------------------------
cd meshes
mkdir ../tmp_meshes
cp rectangular_quadmesh.template.h simple_rectangular_quadmesh.template.h ../tmp_meshes
cp rectangular_quadmesh.template.cc simple_rectangular_quadmesh.template.cc ../tmp_meshes
$del *.h *.cc
mv ../tmp_meshes/* .
rmdir ../tmp_meshes
cp $home_dir/config/stripped_down_files/mesh_names.list .
cp $home_dir/config/stripped_down_files/Makefile.am.meshes Makefile.am





# Update configure scripts
#=========================
cd $home_dir
cd config/configure.ac_scripts
cp $home_dir/config/stripped_down_files/user_drivers.dir_list .
cp $home_dir/config/stripped_down_files/user_src.dir_list .
cp $home_dir/config/stripped_down_files/core.dir_list .
cp $home_dir/config/stripped_down_files/doc.dir_list .



