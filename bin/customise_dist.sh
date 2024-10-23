#! /bin/sh

set -o errexit
set -o nounset

#------------------------------------------------------------------------------
# Shell script to customise distribution -- gives user the option to
# wipe/re-install selected sub-directories before tar-ing up.
# -----------------------------------------------------------------------------


# Number of command line arguments must be two
if (test $# != 3); then
    echo " "
    echo "This script needs three command line arguments: The original"
    echo "oomph-lib home directory path, the dist directory path and"
    echo "the relative path from the original to the dist directory. "
    echo " "
    echo "exiting... "
    exit 5
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

RegenerateConfigFiles()
{
    set -o errexit
    set -o nounset

    oomph_dist_dir="$1"

    # Use autogen.sh to generate an updated configure.ac etc. with the
    # correct list of Makefiles.
    cd "$oomph_dist_dir"
    ./non_interactive_autogen.sh -o -r

    # automake has problems if the packaged distribution doesn't contain a
    # doc dir, but the doc dir won't be copied if --enable-supress-doc is
    # set in the configure options. So we create the dir here.
    mkdir -p "doc"
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
echo " 2: Distribution including doc, including validata"
echo " 3: Distribution without   doc, including validata"
echo " 4: Distribution including doc, without   validata"
echo " 5: Distribution without   doc, without   validata"
echo " "
echo "====================================================================="
echo " "
echo "Choose customisation flag [0,1,2,3,4,5 -- default: 0] "
customisation_flag=`OptionRead 0`

if (test ! \( \( $customisation_flag -eq 0 \) -o \( $customisation_flag -eq 1 \) -o \( $customisation_flag -eq 2 \) -o \( $customisation_flag -eq 3 \) -o \( $customisation_flag -eq 4 \) -o \( $customisation_flag -eq 5 \) \) ) ; then
    echo "unrecognised option $customisation_flag"
    exit 5
fi


# Default flags for all default customisations
#---------------------------------------------
if !(test $customisation_flag -eq 0 ) ; then
    keep_nondist_figures="n"
    keep_svn="n"
    wipe_hsl="y"
    wipe_user_drivers="y"
    wipe_user_src="y"
    include_private_directories="n"
    if (test $customisation_flag -eq 2 ) ; then
        wipe_doc="n"
        wipe_validata="n"
    fi
    if (test $customisation_flag -eq 3 ) ; then
        wipe_doc="y"
        wipe_validata="n"
    fi
    if (test $customisation_flag -eq 4 ) ; then
        wipe_doc="n"
        wipe_validata="y"
    fi
    if (test $customisation_flag -eq 5 ) ; then
        wipe_doc="y"
        wipe_validata="y"
    fi
fi


# Default flags for "package up everything"
#------------------------------------------
echo "Flag: " $customisation_flag
if (test $customisation_flag -eq 1 ) ; then
    echo "Packaging up the lot..."
    keep_nondist_figures="y"
    keep_svn="y"
    wipe_hsl="n"
    wipe_user_drivers="n"
    wipe_user_src="n"
    wipe_doc="n"
    wipe_validata="n"
    include_private_directories="y"
fi

#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////


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
    keep_nondist_figures=`OptionRead n`

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
    echo "====================================================================="
    echo " "
    echo "All demo driver directories contain subdirectories called validata"
    echo "that contain reference data that is used to check the correctness"
    echo "of the computed results during oomph-lib's self-test procedures."
    echo " "
    echo "To reduce the size of the distribution you may wish to exclude this "
    echo "data (though you won't be able to run sensible self-tests)."
    echo " "
    echo "====================================================================="


    # Wipe contents of validata directories?
    #---------------------------------------
    echo " "
    echo " "
    echo "Do you wish to wipe the contents of the validata "
    echo "directories? [y/n - default: n]"
    echo " "
    wipe_validata=`OptionRead n`

fi

if test "$wipe_validata" = "y" -o "$wipe_validata" = "Y" ; then
    echo "Wiping the contents of the validata directories."
    (cd $dist_dir; for dir in `find . -name validata`; do cd $dir; rm -rf * ; cd $dist_dir; done)
else
    echo "Not wiping contents of the validata directory."
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
    keep_svn=`OptionRead n`

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
            echo "       Do you want to wipe? [y/n -- default: n]"
            echo " "
            wipe_hsl=`OptionRead n`
        fi
        if test "$wipe_hsl" = "y"  ; then
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




#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////




# Wipe doc directory?
#--------------------
if (test $customisation_flag -eq 0) ; then
    echo " "
    echo " "
    echo "doc directory: "
    echo "        Do you want to wipe? [y/n -- default: n]"
    echo " "
    wipe_doc=`OptionRead n`
fi
if test "$wipe_doc" = "y"  ; then
    echo "Wiping " `pwd`/doc
    rm -rf `pwd`/doc
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
        echo "        Do you want to include? [y/n -- default: n]"
        echo " "
        include_private_directories=`OptionRead n`
    fi
else
    include_private_directories="n"
fi
if test "$include_private_directories" = "y"; then
    echo " "
    if test "$keep_svn" = "n"; then
        rm -rf `find $dist_dir/private -name .svn`
    fi
    echo "done"
else
    echo "Not including " `pwd`/private
    rm -rf $dist_dir/private
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
    echo "       Do you want to wipe the non-demo ones? [y/n --default: n]"
    echo " " 
    echo "[Note that if directories are not (manually) included in the"
    echo " user_drivers/Makefile.am SUBDIRS variable then they will automatically"
    echo " not be included in the tar file.]"
    wipe_user_drivers=`OptionRead n`
fi


if test "$wipe_user_drivers" = "y"  ; then
    echo "Wiping non-demo directories in  `pwd`/user_drivers"
    RETURN_DIR=`pwd`
    echo $RETURN_DIR
    cd `pwd`/user_drivers
    # Get all files in user driver directory
    FILE_AND_DIR_LIST=`ls -d *`
    # Temporary storage for the directories
    touch oomphs_own_tmp_directory
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
    # Move the demo user driver directories back from temporary storage and
    # clean up
    mv oomphs_own_tmp_directory/jack_cool .
    mv oomphs_own_tmp_directory/joe_cool .
    rm -rf  oomphs_own_tmp_directory
    cd $RETURN_DIR
else
    echo "Not wiping " `pwd`/user_drivers
fi


# Wipe non-demo user_src directories?
#------------------------------------
if (test $customisation_flag -eq 0) ; then
    echo " "
    echo " "
    echo "user_src directory: "
    echo "       Do you want to wipe the non-demo ones? [y/n --default: n]"
    echo " "
    echo "[Note that if directories are not (manually) included in the"
    echo " user_src/Makefile.am SUBDIRS variable then they will automatically"
    echo " not be included in the tar file.]"
    wipe_user_src=`OptionRead n`
fi
if test "$wipe_user_src" = "y"  ; then

    echo "Wiping non-demo directories in  `pwd`/user_src"
    RETURN_DIR=`pwd`
    echo $RETURN_DIR
    cd `pwd`/user_src
    # Get all files in user src directory
    FILE_AND_DIR_LIST=`ls -d *`
    # Temporary storage for the directories
    touch oomphs_own_tmp_directory
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
    echo "done"

else
    echo "Not wiping " `pwd`/user_src
fi



# Done with the customisation, now do some clean up
# ============================================================

# Remove the temporary empty frontal.f that
# might have been created in the original tree If this is not done the
# configure scripts thinks that we do have the hsl sources and
# builds the wrong thing. This is only a problem is one does a make dist
# and then reruns the autogen stuff, but that could happen.
frontal_file=$orig_dir/external_src/oomph_hsl/frontal.f
if ( [ -e $frontal_file ] && [ ! -s $frontal_file ] ); then
    rm $frontal_file;
fi;


# Automake (Makefile.am in oomph root) requires that the doc dir at least
# exists, so create it
mkdir -p "$dist_dir/doc"


# End the customisation procedure by regenerating the config files
RegenerateConfigFiles "$dist_dir"

# ...but replace the symbolic links by the files themselves
mkdir tmp_junk
find . -type l -exec cp {} tmp_junk \;
find . -type l -exec rm -f {} \;
mv tmp_junk/* . || true # allow this to fail if there are no files in tmp_junk/*
rm -rf tmp_junk
