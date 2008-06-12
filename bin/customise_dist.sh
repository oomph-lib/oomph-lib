#! /bin/sh

#--------------------------------------------------------------
# Shell script to customise distribution -- gives 
# user the option to wipe/re-install selected sub-directories before
# tar-ing up. Input: Original directory and dist directory.
#--------------------------------------------------------------



# Number of command line arguments must be two
if (test $# != 3); then 
    echo " "
    echo "This script needs three command line arguments: The original"
    echo "oomph-lib home directory path, the dist directory path and"
    echo "the relative path from the original to the dist directory. "
    echo " "
    echo "exiting... "
    exit
else 
    echo " " 
    orig_dir=$1
    echo "Orig directory         : " $orig_dir
    dist_dir=$2
    echo "Dist directory         : " $dist_dir
    rel_dist_dir=$3
    echo "Relative dist directory: " $rel_dist_dir

fi


echo " "
echo "================================================================="
echo " "
echo "                  Customising the distribution"
echo " "
echo "================================================================="
echo " "
echo " "

# All this needs to run in the dist directory!
cd $dist_dir


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
echo "====================================================================="
echo " "
echo "There are several default customations (if that makes sense...)"
echo "which package up either everything or selected parts of the  "
echo "library. You can either choose these or perform the customisation " 
echo "manually."
echo " "
echo "The options are:"
echo " " 
echo " 0: Manual customisation [default]"
echo " 1: Include the lot... mainly useful for developers"
echo " 2: Distribution including doc" 
echo " 3: Distribution without doc" 
echo " "
echo "====================================================================="
echo " " 
echo "Choose customisation flag [0,1,2,3 -- default: 0] "
reply=`OptionRead`

if (test \( $reply -eq 1 \) -o \( $reply -eq 2 \) -o \( $reply -eq 3 \) ) ; then 
    customisation_flag=$reply
else
    customisation_flag=0
fi
echo " " 



# Default flags for all default customisations
#---------------------------------------------
if !(test $customisation_flag -eq 0 ) ; then 
    keep_nondist_figures="n"
    keep_svn="n"
    wipe_hsl="y"
    wipe_arpack="y"
    wipe_user_drivers="y"
    wipe_user_src="y"
    include_private_directories="n"
    if (test $customisation_flag -eq 2 ) ; then 
        wipe_doc="n"
    fi
    if (test $customisation_flag -eq 3 ) ; then 
        wipe_doc="y"
    fi
fi
     

# Default flags for "package up everything" 
#------------------------------------------
# [NOTE: "t" means "terminate" customisation]
#--------------------------------------------
echo "Flag: " $customisation_flag
if (test $customisation_flag -eq 1 ) ; then 
    echo "Packaging up the lot..."
    keep_nondist_figures="y"
    keep_svn="y"
    wipe_hsl="t"
    wipe_arpack="t"
    wipe_user_drivers="t"
    wipe_user_src="t"
    wipe_doc="n"
    include_private_directories="t"
fi


if (test $customisation_flag -eq 0) ; then 


echo " "
echo "====================================================================="
echo " "
echo "Most doc directories contain subdirectories called  nondist_figures"
echo "that contain tecplot macros/layout files etc. that are used "
echo "(and remain useful!) during the development/maintenance of the" 
echo "documentation."
echo " "
echo "Developers may wish to retain these directories in the distribution."
echo " "
echo "====================================================================="


# Add nondist_figures directories?
#---------------------------------
echo " "
echo " "
echo "Do you wish to retain the nondist_figures directories"
echo " in the distribution? [y/n - default: n]"
echo " "
keep_nondist_figures=`OptionRead`

fi

if test "$keep_nondist_figures" = "y" -o "$keep_nondist_figures" = "Y" ; then 
   echo "Including nondist_figures directories."
   (cd $orig_dir && pwd && find . -path ./$rel_dist_dir -prune -o -name 'nondist_figures' -exec cp -r {} $dist_dir/{} \; )
else
   echo "Not including nondist_figures directories."
fi



#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////




if (test $customisation_flag -eq 0) ; then 

echo " "
echo "================================================================="
echo " "
echo "oomph-lib uses the revision control system subversion"
echo " "
echo "       http://subversion.tigris.org/ "
echo " "
echo "The revision control directories .svn are of no use " 
echo "in the general distribution and by default they will now" 
echo "be wiped."
echo " "
echo "Developers may wish to retain the subversion directories."
echo " "
echo "================================================================="


# Retain subversion directories?
#-------------------------------
echo " "
echo " "
echo "Do you wish to retain the subversion directories in "
echo "the distribution? [y/n -- default: n]"
echo " "
keep_svn=`OptionRead`

fi


if (test "$keep_svn" = "y" -o "$keep_svn" = "Y") ; then 
   echo "Retaining subversion directories:"
   echo "Orig directory         : " $orig_dir
   echo "Dist directory         : " $dist_dir
   echo "Relative dist directory: " $rel_dist_dir
   (cd $orig_dir && find . -path ./$rel_dist_dir -prune -o -name '.svn' -exec cp -rf {}  $dist_dir/{} \; )

else
   echo "Removing subversion directories in " $dist_dir
   rm -rf `find $dist_dir -name .svn`
fi




#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////


if (test $customisation_flag -eq 0) ; then 


echo " "
echo "================================================================="
echo " "
echo "      We'll now run through a number of directories/files that you might"
echo "      not want to include in the distribution. You will need "
echo "      to explicitly type \"y\" to remove the indicated directory/file."
echo "      \"t\" terminates the selection processes -- all remaining"
echo "      directories will then automatically get included."
echo " "
echo "================================================================="


fi

# Wipe hsl sources?
#------------------

# Does the hsl file exist?
if test -e `pwd`/external_src/oomph_hsl/frontal.f; then
    # If it exists and has zero length, kill it as it was only
    # created (with touch) as a dummy during make dist
    if test -s `pwd`/external_src/oomph_hsl/frontal.f; then

        if (test $customisation_flag -eq 0) ; then 
            echo " "
            echo " "
            echo "hsl sources in frontal.f: "
            echo "       Do you want to wipe? [y/n/t -- default: n]"
            echo " "
            wipe_hsl=`OptionRead`
        fi

        if test "$wipe_hsl" = "t"  ; then 
            echo " "
            echo "=========================================================="
            echo " "
            echo "Terminating customisation of distribution."
            echo "Remaining directories/files get included automatically"
            echo " "
            echo "==========================================================" 
            echo " " 
            exit
        fi
        if test "$wipe_hsl" == "y"  ; then 
            echo "Wiping " `pwd`/external_src/oomph_hsl/frontal.f 
            rm -f `pwd`/external_src/oomph_hsl/frontal.f 
        else
            echo "Not wiping " `pwd`/external_src/oomph_hsl/frontal.f  
        fi
    else
        # Wipe the dummy 
        echo "Wiping the empty dummy file " `pwd`/external_src/oomph_hsl/frontal.f  
        rm -f `pwd`/external_src/oomph_hsl/frontal.f 
    fi
fi






#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////



# Wipe ARPACK sources?
#---------------------

# Does the ARPACK file exist?
if test -e `pwd`/external_src/oomph_arpack/all_arpack_sources.f; then
    # If it exists and has zero length, kill it as it was only
    # created (with touch) as a dummy during make dist
    if test -s `pwd`/external_src/oomph_arpack/all_arpack_sources.f; then

        if (test $customisation_flag -eq 0) ; then 
            echo " "
            echo " "
            echo "ARPACK sources in all_arpack_sources.f: "
            echo "    Do you want to wipe? [y/n/t -- default: n]"
            echo " "
            wipe_arpack=`OptionRead`
        fi
        if test "$wipe_arpack" = "t"  ; then 
            echo " "
            echo "=========================================================="
            echo " "
            echo "Terminating customisation of distribution."
            echo "Remaining directories/files get included automatically"
            echo " "
            echo "==========================================================" 
            echo " " 
            exit
        fi
        if test "$wipe_arpack" == "y"  ; then 
            echo "Wiping " `pwd`/external_src/oomph_arpack/all_arpack_sources.f 
            rm -f `pwd`/external_src/oomph_arpack/all_arpack_sources.f 
        else
            echo "Not wiping " `pwd`/external_src/oomph_arpack/all_arpack_sources.f  
        fi
    else
        # Wipe the dummy 
        echo "Wiping the empty dummy file " `pwd`/external_src/oomph_arpack/all_arpack_sources.f  
        rm -f `pwd`/external_src/oomph_arpack/all_arpack_sources.f 
    fi
fi



#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////




# Wipe doc directory?
#--------------------
if (test $customisation_flag -eq 0) ; then 
    echo " "
    echo " "  
    echo "doc directory: "
    echo "        Do you want to wipe? [y/n/t -- default: n]"
    echo " "
    wipe_doc=`OptionRead`
fi
if test "$wipe_doc" = "t"  ; then 
  echo " "
  echo "================================================================="
  echo " "
  echo "Terminating customisation of distribution."
  echo "Remaining directories/files get included automatically"
  echo " "
  echo "=================================================================" 
  echo " " 
  exit
fi
if test "$wipe_doc" == "y"  ; then 
   echo "Wiping " `pwd`/doc 
   rm -rf `pwd`/doc
   rm -f `pwd`/config/configure.ac_scripts/doc.dir_list
   echo "done"
else
   echo "Not wiping " `pwd`/doc 
fi



#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////


# Wipe private directories?
#--------------------------
if (test -d $orig_dir/private); then
if (test $customisation_flag -eq 0) ; then 
    echo " "
    echo " "  
    echo "private directories: "
    echo "        Do you want to include? [y/n/t -- default: n]"
    echo " "
    include_private_directories=`OptionRead`
fi
else
  include_private_directories="n"
fi
if test "$include_private_directories" = "t"; then
  echo " "
  echo "================================================================="
  echo " "
  echo "Terminating customisation of distribution."
  echo "Remaining directories/files get included automatically"
  echo " "
  echo "=================================================================" 
  echo " " 
  (cd $orig_dir && pwd && cp -r private $dist_dir)
  if test "$keep_svn" = "n"; then
   rm -rf `find $dist_dir/private -name .svn`
  fi
  exit
fi
if test "$include_private_directories" = "y"; then
   echo " "
   (cd $orig_dir && pwd && cp -r private $dist_dir)
   if test "$keep_svn" = "n"; then
    rm -rf `find $dist_dir/private -name .svn`
   fi
   echo "done"
else
   echo "Not including " `pwd`/private
   rm -f `pwd`/config/configure.ac_scripts/private_user_drivers.dir_list
   rm -f `pwd`/config/configure.ac_scripts/private_user_src.dir_list
   rm -f `pwd`/config/configure.ac_scripts/private.dir_list
fi

#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////


# Wipe non-demo user_drivers directory?
#--------------------------------------
if (test $customisation_flag -eq 0) ; then 
    echo " "
    echo " "
    echo "user_drivers directories: "
    echo "       Do you want to wipe the non-demo ones? [y/n/t --default: n]"
    echo " "
    wipe_user_drivers=`OptionRead`
fi
if test "$wipe_user_drivers" = "t"  ; then 
  echo " "
  echo "================================================================="
  echo " "
  echo "Terminating customisation of distribution."
  echo "Remaining directories/files get included automatically"
  echo " "
  echo "=================================================================" 
  echo " "
  #Regenerate the config files in case other directories have been wiped   
  `pwd`/bin/regenerate_config_files.sh `pwd` 
  exit
fi

if test "$wipe_user_drivers" == "y"  ; then 
   echo "Wiping non-demo directories in  `pwd`/user_drivers" 
   RETURN_DIR=`pwd`
   echo $RETURN_DIR
   cd `pwd`/user_drivers
   # Get all files in user driver directory
   FILE_AND_DIR_LIST=`ls -d *`
   # Temporary storage for the directories
   rm -rf oomphs_own_tmp_directory
   mkdir oomphs_own_tmp_directory
   # Move the directories into temporary storage
   echo $FILE_AND_DIR_LIST
   for file in $FILE_AND_DIR_LIST; do
       if (test -d $file); then
           #echo $file " is a directory"
           mv $file oomphs_own_tmp_directory
       #else
       #    echo $file " is NOT a directory"
       fi
   done
   # Move the demo user driver directories back from temporary storage
   # and clean up
   mv oomphs_own_tmp_directory/jack_cool .
   mv oomphs_own_tmp_directory/joe_cool .
   rm -rf  oomphs_own_tmp_directory
   cd $RETURN_DIR
   # Finally: Overwrite the list of user driver directories in
   # config/configure.ac_scripts/user_drivers.dir_list
   # by default which only contains the two demo ones
   cp -f  `pwd`/config/configure.ac_scripts/user_drivers.dir_list.default \
          `pwd`/config/configure.ac_scripts/user_drivers.dir_list 

else
   echo "Not wiping " `pwd`/user_drivers 
fi


# Wipe non-demo user_src directories?
#------------------------------------
if (test $customisation_flag -eq 0) ; then 
    echo " "
    echo " "
    echo "user_src directory: "
    echo "       Do you want to wipe the non-demo ones? [y/n/t]"
    echo " "
    wipe_user_src=`OptionRead`
fi
if test "$wipe_user_src" = "t"  ; then 
  echo " "
  echo "================================================================="
  echo " "
  echo "Terminating customisation of distribution."
  echo "Remaining directories/files get included automatically"
  echo " "
  echo "=================================================================" 
  echo " "
  #Regenerate the config files in case other directories have been wiped   
  `pwd`/bin/regenerate_config_files.sh `pwd`
  exit
fi
if test "$wipe_user_src" == "y"  ; then 

   echo "Wiping non-demo directories in  `pwd`/user_src" 
   RETURN_DIR=`pwd`
   echo $RETURN_DIR
   cd `pwd`/user_src
   # Get all files in user src directory
   FILE_AND_DIR_LIST=`ls -d *`
   # Temporary storage for the directories
   rm -rf oomphs_own_tmp_directory
   mkdir oomphs_own_tmp_directory
   # Move the directories into temporary storage
   echo $FILE_AND_DIR_LIST
   for file in $FILE_AND_DIR_LIST; do
       if (test -d $file); then
           #echo $file " is a directory"
           mv $file oomphs_own_tmp_directory
       #else
       #    echo $file " is NOT a directory"
       fi
   done
   # Move the demo user src directories back from temporary storage
   # and clean up
   mv oomphs_own_tmp_directory/jack_cool .
   rm -rf  oomphs_own_tmp_directory
   cd $RETURN_DIR
   # Finally: Overwrite the list of user src directories in
   # config/configure.ac_scripts/user_src.dir_list
   # by default which only contains the demo one
   cp -f  `pwd`/config/configure.ac_scripts/user_src.dir_list.default \
          `pwd`/config/configure.ac_scripts/user_src.dir_list 
   echo "done"

else
   echo "Not wiping " `pwd`/user_src
fi



#End the customisation procedure by regenerating the config files
`pwd`/bin/regenerate_config_files.sh `pwd` 


