#! /bin/bash

set -o errexit


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
 echo $Opt
}

#This little function takes the input, removes anything following a #
#deletes blanks lines and then replaces all newlines by spaces
ProcessOptionsFile()
{
echo `cat $1 | sed 's/#.*$//' | sed '/^$/d' | tr '\012' ' '`
}



#Check that the "--" options preceed other configure options.
CheckOptions()
{
   awk '
   BEGIN{encountered_first_minus_minus=0
         encountered_first_non_minus_minus_after_first_minus_minus=0}
   NF {# pattern NF ignores blank lines since it expands to 0 for empty lines!
   # Ignore any comments (first entry in row starts with "#")
   if (substr($1,1,1)!="#")
    { 
     # Does the first entry in the line start with "--"?
     if (substr($1,1,2)=="--"){encountered_first_minus_minus=1}

     # Have we encountered the first "--" entry already?
     if (encountered_first_minus_minus==1)
       {
        # Does the first entry in the line not start with "--"?
        if (substr($1,1,2)!="--")
         {
          encountered_first_non_minus_minus_after_first_minus_minus=1
         }
       }
     # Should if this is followed by another "--" entry!
     if ((encountered_first_minus_minus==1)&&
         (encountered_first_non_minus_minus_after_first_minus_minus==1))
      {
       if (substr($1,1,2)=="--")
        {
         ok=0
         print "ERROR: The entry\n\n" $0 "\n"
         print "is in the wrong place. All the \"--\" prefixed options should go first!\n\n"
        }
      }
    }
    }' `echo $1`


#old echo `echo -n $@ | sed 's/^/ /' | sed -n '/[ ].[^-].* --/p'`
#echo `echo $@ | sed 's/^/ /' | sed -n '/[ ].[^-].* --/p'`
# echo `printf "$@" | sed 's/^/ /' | sed -n '/[ ].[^-].* --/p'`
}


#This function returns a match (non-null string) if the input string 
#contains a long option (starting with --) after a short option 
#(no --)
#CheckOptions()
#{
##old echo `echo -n $@ | sed 's/^/ /' | sed -n '/[ ].[^-].* --/p'`
#echo `echo $@ | sed 's/^/ /' | sed -n '/[ ].[^-].* --/p'`
# echo `printf "$@" | sed 's/^/ /' | sed -n '/[ ].[^-].* --/p'`
#}

#This little function echo's the usage information
EchoUsage()
{
        echo "Usage: "
        echo "------ "
        echo " "
        echo "[without flags]: Normal \"./configure; make; make install; make check -k\" sequence."
        echo " "
        echo " --rebuild     : Complete re-configure, followed by normal build sequence."
        echo " "
        echo "--jobs[=N]     :  Run N make jobs simultaneously."
        echo "                  Useful for speeding up the build on multi-core processors." 
        exit
}




#====================================================================
# Start Q/A session
#====================================================================

echo " "
echo "============================================================= "
echo "              oomph-lib installation script" 
echo "============================================================= "
echo " "

# Do you want to rebuild from scratch?
#-------------------------------------
#If so specify --rebuild as command line argument. Default is 
# to just do the normal ".configure, make, make install, make check -k" sequence.

#Bail out if more than two command line arguments
if (test $# -gt 2); then 
 EchoUsage 
fi   

#Process the command line options
raw_build=false;
make_options=" ";
while (test $# -gt 0)
do
   case "$1" in
     #Set the rebuild flag
     --rebuild) 
      echo "             [Doing complete rebuild from scratch.]"
      raw_build=true;;
     #Set the jobs flag
     --jobs*)
      make_options="$1";;
     #Anything else bail out     
      *)  
       EchoUsage;;
   esac
   shift
done

if (test "$raw_build" = "false"); then
   echo "                     [Doing normal build.]"
fi   


# Set the script to crash if any un set variables are used (we put this after
# the options processsing since some command line arguments may legitimately not
# exist).
set -o nounset

# Read out root install directory
#--------------------------------
MY_HOME_WD=`pwd`



# If this is a rebuild: Check for helper scripts 
#-----------------------------------------------
if $raw_build; then

  SCRIPT_LIST=`echo config.guess config.sub depcomp install-sh ltmain.sh missing aclocal.m4 mkinstalldirs `
  SCRIPTS_EXIST="no"
  for script in $SCRIPT_LIST
   do
    if (test -e $script); then
        SCRIPTS_EXIST="yes"
    fi 
  done
  if test "$SCRIPTS_EXIST" = "yes" ; then 
    echo " "
    echo "You may wipe the symbolic links to the autoconf/automake helper scripts"
    echo " "
    for script in $SCRIPT_LIST
      do
        if (test -e $script); then
          echo "   " $script
        fi 
      done
    echo " "
    echo "[This is recommended if you have moved the sources to a different"
    echo " machine without packaging them up with make dist. The symbolic "
    echo " links tend to be machine-specific so it's best to force "
    echo " autoconf/automake to rebuild them on the new machine]."
    echo " "
    echo " Do you want to wipe the helper scripts [y/n -- default: n]"
    reply=`OptionRead`
    if test "$reply" = "y" -o "$reply" = "Y" ; then 
        echo " "
        echo "As a backup: Here are the old symbolic links:"
        echo " "
        for script in $SCRIPT_LIST
          do
          if (test -L $script); then
              ls -L $script
              ls -l $script > old_symbolic_links.txt
          fi
        done
        echo " "
        echo "We have stored this information in old_symbolic_links.txt"
        echo " "
        echo "Wiping them..."
        rm -f  $SCRIPT_LIST
        echo "Done"
    fi   
else
    echo " "
    echo "[No autoconf/automake helper scripts to be wiped...]"
    echo " "
fi
fi


# Autodetect folders in user_drivers
#-----------------------------------

# Backup old file (use -f so it doesn't give an error if the file doesn't exist)
touch ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list
mv -f ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list.backup

# Get a list of locations of files named Makefile.am, modify a little and write to user_drivers.dir_list.
find ${MY_HOME_WD}/user_drivers -type f -name "Makefile.am" \
    | grep -v "^${MY_HOME_WD}/user_drivers/Makefile.am" \
    | sed 's:/Makefile.am$::' \
    | sed "s:^${MY_HOME_WD}/::" \
    > ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list

# The grep and sed commands above do the following: 1) Remove the line that
# corresponds to the Makefile.am in user_drivers itself. 2) Remove
# "/Makefile.am" from each line leaving only the directory (dirname doesn't work
# with pipes). 3) Remove the start of the path from each line leaving only the
# location relative to the oomph-lib root directory.

echo
echo "User driver folders included are:"
cat ${MY_HOME_WD}/config/configure.ac_scripts/user_drivers.dir_list
echo


# Set the build directory (for lib,include), relative to root
#------------------------------------------------------------
build_sub_dir=build

# Suggested build directory
#--------------------------
build_dir=$MY_HOME_WD/$build_sub_dir

# Check build directory with user
#--------------------------------
echo " "
echo " "
echo "I'm going to install the distribution (the lib and include directories)"
echo "in:"
echo " "
echo "    " $build_dir
echo " "
echo " "
OptionPrompt " Is this OK? [y/n -- default: n]"
reply=`OptionRead`
if test "$reply" != "y" -a "$reply" != "Y" ; then 
   OptionPrompt "Specify build directory [e.g. /home/joe_user/build] :"
   build_dir=`OptionRead`
else
    echo "It's ok"
fi


# Summary of build info
#----------------------
echo " "
echo "============================================================= "
echo " "
echo "Build directory is: " 
echo " " 
echo "     " $build_dir
echo " " 
echo "--> The include directory will be in: "
echo " " 
echo "    " $build_dir"/include" 
echo " " 
echo "--> The lib directory will be in: "
echo " " 
echo "    " $build_dir"/lib" 
echo " "
echo "etc.       " 
echo " "
echo "============================================================= "
echo " "


# Create configure options file
#------------------------------

# If "current" configure options file does not exist then copy in the
# default one:
if (test ! -f config/configure_options/current); then
    cp config/configure_options/default config/configure_options/current
fi

# Process configure options from the file config/configure_options/current
# Ignore any line that starts with "#";
# Add continuation slash at the end of each
# line. Check that all options (starting with "--") 
# come first. 

#Continue asking if the options are OK until approved
accept_configure_options=0
full_list="false"
list_changed="false"
while (test $accept_configure_options -eq 0)
do


#Read the options from the file and convert them into a single one-line string
configure_options=`ProcessOptionsFile config/configure_options/current`


#Check that the options are in the correct order
configure_options_are_ok=`CheckOptions config/configure_options/current`
if test "$configure_options_are_ok" != ""; then

  echo " "
  echo "==============================================================="
  echo "Error message from autogen.sh:"
  echo " " 
  echo $configure_options_are_ok
  echo " " 
  reply="n"
  OptionRead #This is just a pause

#If the options are in the correct order, ask whether they are OK
else

  echo " "
  echo "Configure options are: "
  echo 
  echo $configure_options
  echo 
  if test $list_changed = "false"; then
   OptionPrompt "Is this OK? [y/n -- default: y]"
   reply=`OptionRead`
  else
   reply="n"
   list_changed="false"
  fi
fi


#If it's not OK, then read in alternative options from a file, or
#specify on command line
private_configure_option_files=""
if test "$reply" = "n" -o "$reply" = "N"; then
 
  #Remove the current symbolic link (or file)
  #rm -f config/configure_options/current   

  # Link in the private ones:
  return_dir_before_link_in_private=`pwd`
  # Kill stray symlinks
  cd config/configure_options
  find . -type l -exec rm {} \; 
  if test  "$full_list" = "true"; then
   cd private_configure_options
   private_configure_option_files=`ls `
   cd ..
   for file in $private_configure_option_files; do ln -s private_configure_options/$file ; done
  fi
  cd $return_dir_before_link_in_private

  # Ooops: Non-portable gnu extension to ls
  #configure_option_files=`ls --ignore=private_configure_options config/configure_options`

  # Thanks for this fix, Andy!
  configure_option_files=`ls config/configure_options | grep -v  private_configure_options` 

  echo " "
  echo "======================================================================"
  echo 
  echo "Choose an alternative configuration file "
  #Loop over files and display a menu
  count=0
  for file in $configure_option_files
   do
    #Increase the counter
    count=`expr $count + 1`
    echo $count ": " $file
   done #End of loop over files in config/configure_options
 echo

  echo "Enter the Desired configuration file [1-"$count"]"
  if test $full_list = "false"; then
   echo "Enter -1 to show an extended list of options"
  else 
   echo "Enter -1 to show a short list of options"
  fi
  echo "Enter 0 to specify the options on the command line"
  #Read in the Desired File
  file_number=`OptionRead`

  #If options are to be read from the command line then store the#
  #options in the file config/configure_options/current
  if (test $file_number -eq 0); then
   echo 
   echo "Enter options"
   configure_options=`OptionRead`  
   echo $configure_options > config/configure_options/current

  #Otherwise copy the desired options file to config/configure_options/current
  elif (test $file_number -eq -1); then
   list_changed="true"
   if test $full_list = "true"; then
    full_list="false"
   else
    full_list="true"
   fi
  else   
   #Reset the counter
   count=0
   #Loop over the files until the counter equals the chosen file_number
   for file in $configure_option_files
     do
     #Increase the counter
     count=`expr $count + 1`
     if (test $count -eq $file_number); then
        cp -f config/configure_options/$file config/configure_options/current
        break
     fi
   done #End of loop over files
   fi #End of create symbolic link code

#If the configuration is OK, accept it
else
 echo " " 
 echo "Configure options have been accepted."
 accept_configure_options=1
fi

done #End of while loop over customisation of configure options


# Undo links to private configure options
return_dir_before_link_in_private=`pwd`
cd config/configure_options
for file in $private_configure_option_files; do rm -f $file; done
cd $return_dir_before_link_in_private

echo " "
echo "==================================================================="
echo " "
echo " "
echo "End of customisation -- the actual build process is about to start."
echo "This may take a while... No user intervention is required during"
echo "the build procedure, so go and take the dog for a walk..."
echo " " 
OptionPrompt "Hit enter to continue."
tmp=`OptionRead`





#====================================================================
# Start actual build process
#====================================================================


# If we are doing a raw build or if ./configure does not yet exist then generate
# all config files needed.
#--------------------------------------------------------
if [ $raw_build -o ! -e ./configure ]; then
 $MY_HOME_WD/bin/regenerate_config_files.sh $MY_HOME_WD
fi

# Run configure command
echo " "
echo "Running ./configure --prefix $build_dir $configure_options"
echo " " 
/bin/sh -c "./configure --prefix $build_dir $configure_options"

echo " " 
echo " " 
echo "done running ./configure"
echo " " 
echo " " 

# Test that the mpi commands work (automatically passes if no variable
# MPI_RUN_COMMAND in makefile). This needs to go after configure so that we
# can use the generated Makefile to (robustly) get the run and compile
# commands.
$MY_HOME_WD/bin/check_mpi_command.sh $MY_HOME_WD/Makefile


# Make all libraries
#-------------------
echo " "
echo "Running make $make_options" 
make $make_options
echo "done"


# Install the libraries (in build directory specified above)
#-----------------------------------------------------------
echo " "
echo "running make $make_options install"
make $make_options install
echo "done" 

echo " "
echo "autogen.sh has finished! If you can't spot any error messages" 
echo "above this, oomph-lib should now be ready to use... " 
echo " " 
echo "If you encounter any problems, please study the installation" 
echo "instructions and the FAQ before contacting the developers. " 
echo " "
echo "To run self tests use \"make check -k\" or ./bin/parallel_self_test.py"
