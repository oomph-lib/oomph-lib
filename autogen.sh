#! /bin/sh

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
   {
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
    }' config/configure_options/current


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


# Wipe previous builds
#---------------------
if (test -d  $build_dir); then 
    echo " "
    echo "Note: Build directory " $build_dir " exists."
    echo " "
    OptionPrompt "Do you want to wipe it [y/n -- default: n]"
    reply=`OptionRead`
    if test "$reply" = "y" -o "$reply" = "Y" ; then 

       echo " "
       echo "Sorry to be over-cautious here, but we'd better double-check"
       echo "before we delete some of your precious files..."
       echo " "
       echo "The contents of " $build_dir " are:"
       echo " "
       ls -l  $build_dir
       echo " "
       echo "What you should see above are the include and lib directories, "
       echo "nothing else! "
       echo " "
       OptionPrompt "Are you still sure you want to wipe it [y/n -- default: n]"
       reply2=`OptionRead`
       if test "$reply2" = "y" -o "$reply2" = "Y" ; then 
          echo " "
          echo "Wiping it..."
          rm -f -r $build_dir
          echo "Done"
       fi
    fi
fi





#Doc size of build?
#------------------
echo " " 
echo " " 
OptionPrompt "Do you want to document the size of the installation ? [y/n -- default: n]"
want_doc_size=`OptionRead`
if test "$want_doc_size" = "y" -o "$want_doc_size" = "Y" ; then 
    doc_size="y"
    echo " " 
    echo "Size of the distribution at various stages of the"
    echo "build process will be documented in size.html"
    echo " " 
else
    doc_size="n"
    echo " " 
    echo "Not documenting size of build"
fi

#Create header for size.html file
#--------------------------------
if test "$doc_size" = "y"; then 
    echo " "
    echo "Assessing size of the distribution before the"
    echo "installation -- this can take a little while... "
    echo " " 
    echo "<html>" > size.html
    echo "<h1>" >> size.html
    echo "Build sizes for oomph-lib build in " `pwd` " on " `date` >> size.html
    echo "</h1>" >> size.html
    echo "<TABLE BORDER=1>" >> size.html
    echo "<TR>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "tar file" >> size.html
    echo "</TD>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "Documentation?" >> size.html
    echo "</TD>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "Validata?" >> size.html
    echo "</TD>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "Size of gzipped tar file" >> size.html
    echo "</TD>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "Size of unpacked distribution before build" >> size.html
    echo "</TD>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "Size after build" >> size.html
    echo "</TD>" >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "Size after self-tests" >> size.html
    echo "</TD>" >> size.html
    echo "</TR>" >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo "<TD align=\"center\">" >> size.html
    echo "<a href=\"oomph-lib-0.90.tar.gz\">oomph-lib-0.90.tar.gz</A>" >> size.html
    echo "</TD>" >> size.html


    # Indicate if we have documentation available
    #--------------------------------------------
    echo "<TD align=\"center\">" >> size.html
    if (test -e doc/doc.txt); then 
      echo "yes" >> size.html
    else 
      echo "no" >> size.html
    fi
    echo "</TD>" >> size.html


    # Indicate if we have validata available
    #---------------------------------------
    echo "<TD align=\"center\">" >> size.html
    if (test -e demo_drivers/poisson/two_d_poisson_flux_bc_adapt/validata/results.dat.gz); then 
      echo "yes" >> size.html
    else 
      echo "no" >> size.html
    fi
    echo "</TD>" >> size.html


    # Dummy output for size of gzipped tar file
    #------------------------------------------
    echo "<TD align=\"center\">" >> size.html
    echo "***G" >> size.html
    echo "</TD>" >> size.html


    # Size before build
    #------------------
    echo "<TD align=\"center\">" >> size.html
    echo `du -s -h .` >> size.html
    echo "</TD>" >> size.html

    echo " "
    echo "...done."
    echo " " 

fi


# Create configure options file
#------------------------------

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
echo "=============================================================" 
echo " " 


# Automatic run of self tests? 
#-----------------------------
echo " " 
echo "It is possible to execute the library's extensive self-test"
echo "when the build is complete. Depending on the speed of your"
echo "computer, this may take a long time as more than 450 test "
echo "codes are run and validated."
echo " " 
echo "Do you want to build the demo codes and run the self-tests at the end "
OptionPrompt "of the build procedure? [y/n -- default: n]"
want_self_tests=`OptionRead`

if test "$want_self_tests" = "y" -o "$want_self_tests" = "Y" ; then 
  echo " "
  echo "\"make check -k\" will be run at the end of the build procedure"
  echo "to build/self-test the demo codes."
  echo " "
  echo "The results of the self-tests will be written to the"
  echo "file validation.log in the top-level oomph-lib directory."
  echo " "
  if [ -e validation.log ]; then 
    echo " " 
    OptionPrompt "The file validation.log exists. Is it OK to wipe it? [y/n -- default: n]"
    reply=`OptionRead`
    if test "$reply" != "y" -a "$reply" != "Y" ; then 
       echo " "
       echo "In that case I am disabling the automatic self-test that you requested earlier."
       want_self_tests="n"
       echo " "
       OptionPrompt "Hit enter to acknowledge"
       tmp=`OptionRead`
    else
       rm -f validation.log
       echo "I have wiped the validation.log file."
    fi
  fi
else
  echo " "
  echo "The self-test procedure will not be run automatically."
  echo "You may initiate the self-tests manually by typing \"make check -k\" "
  echo "in the top-level oomph-lib directory when the build process "
  echo "is complete."
  echo " "
fi


# If mpi self-tests are requested, check if we can compile
#---------------------------------------------------------
# and run a sample mpi code
#---------------------------
if test "$want_self_tests" = "y" -o "$want_self_tests" = "Y" ; then 

  # Extract c++ compilation command
  cxx_compile_command=`echo $configure_options | \
    awk '{ \
    start=match($0,"CXX="); \
    if (start!=0) \
      {  \
        start+=4; rest=substr($0 ,start); \
        end=match(rest," ")-1; \
      } \
    } \
    END{if (start!=0) \
    { \
      print  substr(rest,1,end) \
    } \
  else \
    { \
      print "FAILED cxx" \
    } \
  }'` 


  mpi_run_command=`echo $configure_options | \
    awk '{ \
    start=match($0,"--with-mpi-self-tests"); \
    if (start!=0) \
      {  \
        start+=23; rest=substr($0 ,start); \
        end=match(rest,"\"")-1; \
      } \
    } \
    END{if (start!=0) \
    { \
      print  substr(rest,1,end) \
    } \
  else \
    { \
      print "FAILED" \
    } \
  }'` 

  if test "$mpi_run_command" != "FAILED" 
  then  
    echo " "
    echo "============================================================"
    echo " "
    echo "It appears that you want to excecute the mpi self-tests " 
    echo "by executing the mpi codes with the run command:"
    echo " " 
    echo "       " $mpi_run_command
    echo " "
    echo "I'm now going to check if we can compile and run a basic"
    echo "stand-alone mpi code."
    echo " " 

    if test "$cxx_compile_command" != "FAILED cxx" 
    then  
      echo "OK, let's try to compile the basic mpi test"
      echo "code with the compile command: " 
      echo " "            
      full_command=`echo $cxx_compile_command -o bin/minimal_mpi_test bin/minimal_mpi_test.cc` 
      echo "      "$full_command 
      echo " " 
      rm -f bin/minimal_mpi_test
      `echo $full_command` 
      if [ ! -e bin/minimal_mpi_test ]; then
        echo " " 
        echo "================================================="
        echo " "
        echo "WARNING (ISSUED BY OOMPH-LIB):"
        echo "------------------------------"
        echo " "
        echo "Compilation of bin/minimal_mpi_test.cc failed."
        echo "Are you sure your c++ compiler can compile mpi code?"
        echo " "
        echo "Note: This does not necessarily indicate a problem. "
        echo "      autogen.sh tries to extract the c++ compiler"
        echo "      from the configure options assuming that it is"
        echo "      specified in the form CXX=mpic++, say (no quotes,"
        echo "      no spaces)."
        echo " "
        echo "I will continue regardless but you shouldn't be"
        echo "surprised if large numbers of mpi self-test fail..."
        echo " "
        echo "================================================="
        echo " "
      else
        echo "Done: Executable was produced -- that's good!"
        echo " " 
        echo "Now let's run the minimal mpi test to see if mpi is up and running:" 
        echo "I'm going to run the code with the command: "
        full_command=`echo $mpi_run_command ` 
        full_command=$full_command" bin/minimal_mpi_test "
        rm -f bin/minimal_mpi_test.out
        echo " " 
        echo "      " $full_command
        echo " " 
        `echo $full_command`  > bin/minimal_mpi_test.out
        result=`grep 'This worked'  bin/minimal_mpi_test.out | wc | awk '{print $1}'`
        if [ "$result" -ne "2" ]
        then
          echo " " 
          echo "================================================="
          echo " "
          echo "WARNING (ISSUED BY OOMPH-LIB):"
          echo "------------------------------"
          echo " "
          echo "The mpi test code bin/minimal_mpi_test was not run"
          echo "successfully." 
          echo " " 
          echo "You may want to check the following:" 
          echo "-- Are you sure your mpi demons have been started?"
          echo "   E.g. under lam you have to use the lamboot "
          echo "   command to get mpi up and running...."
          echo "-- Is the mpi run command you specified via "
          echo "   the --with-mpi-self-tests flag in the"
          echo "   configure options valid?"
          echo " "
          echo "I will continue regardless but you shouldn't be"
          echo "surprised if large numbers of mpi self-test fail..."
          echo " "
          echo "================================================="
          echo " "
        else
          echo " " 
          echo "Done: mpi test code was executed succesfully. Good stuff."
        fi
      fi
    else
      echo " "
      echo "================================================="
      echo " "
      echo "WARNING (ISSUED BY OOMPH-LIB):"
      echo "------------------------------"
      echo " "
      echo "Sorry I got myself confused when parsing " 
      echo "the configure options and could not find"
      echo "the specification of the c++ mpi compiler"
      echo "via the CXX flag. Please make sure this is"
      echo "specified in configure/configure_optinios/current"
      echo "in the form: CXX=mpic++, say, (no quotes, no spaces)."
      echo " " 
      echo "NOTE: This is not necessarily a problem, but will"
      echo "keep me from checking if mpi is up and running"
      echo "before starting the self-tests."
      echo " "
      echo "================================================="
      echo " "        
    fi
    echo " "
  else
    echo " "
    echo "No mpi self-tests were requested"
    echo " "
  fi  

  # Testing mpi run command on variable number of processors.
  mpi_np_run_command=`echo $configure_options | \
    awk '{ \
    start=match($0,"--with-mpi-self-tests-variablenp"); \
    if (start!=0) \
      {  \
        start+=34; rest=substr($0 ,start); \
        end=match(rest,"\"")-1; \
      } \
    } \
    END{if (start!=0) \
    { \
      print  substr(rest,1,end) \
    } \
  else \
    { \
      print "FAILED" \
    } \
  }'` 

  if test "$mpi_np_run_command" != "FAILED" 
  then  
    echo " "
    echo "============================================================"
    echo " "
    echo "It appears that you want to excecute the mpi self-tests " 
    echo " which involves a varying number of processors by executing the "
    echo " mpi codes with the run command:"
    echo " " 
    echo "       " $mpi_np_run_command
    echo " "

    # Check if OOMPHNP is in the mpi run command.
    OOMPHNP_CHECK=`echo $mpi_np_run_command | \
      awk '{ \
      start=match($0,"OOMPHNP"); \
    } \
    END{if (start!=0) \
    { \
      print "PASSED" \
    } \
  else \
    { \
      print "FAILED" \
    } \
  }'` 

    if test "$OOMPHNP_CHECK" != "FAILED"
    then
      echo "I'm now going to check if we can compile and run a basic"
      echo "stand-alone mpi code on 1, 2, 3 and 4 cores."

      if test "$cxx_compile_command" != "FAILED cxx" 
      then  
        echo "OK, let's try to compile the basic mpi test"
        echo "code with the compile command: " 
        echo " "            
        full_command=`echo $cxx_compile_command -o bin/minimal_mpi_variablenp_test bin/minimal_mpi_variablenp_test.cc` 
        echo "      "$full_command 
        echo " " 


        rm -f bin/minimal_mpi_variablenp_test
        `echo $full_command` 
        if [ ! -e bin/minimal_mpi_variablenp_test ]; then
          echo " " 
          echo "================================================="
          echo " "
          echo "WARNING (ISSUED BY OOMPH-LIB):"
          echo "------------------------------"
          echo " "
          echo "Compilation of bin/minimal_mpi_test.cc failed."
          echo "Are you sure your c++ compiler can compile mpi code?"
          echo " "
          echo "Note: This does not necessarily indicate a problem. "
          echo "      autogen.sh tries to extract the c++ compiler"
          echo "      from the configure options assuming that it is"
          echo "      specified in the form CXX=mpic++, say (no quotes,"
          echo "      no spaces)."
          echo " "
          echo "I will continue regardless but you shouldn't be"
          echo "surprised if large numbers of mpi self-test fail..."
          echo " "
          echo "================================================="
          echo " "
        else
          echo "Done: Executable was produced -- that's good!"
          echo " " 
          echo "Now let's run the minimal mpi test to see if mpi is up and running:" 
          echo "I'm going to run the code with the command: "
          full_command=`echo $mpi_np_run_command `
          full_command=$full_command" bin/minimal_mpi_variablenp_test "
          echo " " 
          echo "      " $full_command
          echo "where OOMPHNP = 1, 2, 3, and 4." 
          echo " " 

          # one processor                                                                 
          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/1/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/2/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/3/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/4/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          rm -f bin/minimal_mpi_variablenp_test

          variablenp_outputs="mpi_seltest_np1rank0
          mpi_seltest_np2rank0 mpi_seltest_np2rank1 \
            mpi_seltest_np3rank0 mpi_seltest_np3rank1 mpi_seltest_np3rank2 \
            mpi_seltest_np4rank0 mpi_seltest_np4rank1 mpi_seltest_np4rank2 mpi_seltest_np4rank3"

          for OUTFILE in $variablenp_outputs 
          do
            result=`grep 'This worked'  $OUTFILE | wc | awk '{print $1}'`

            if [ "$result" -ne "1" ]
            then
              echo " " 
              echo "================================================="
              echo " "
              echo "WARNING (ISSUED BY OOMPH-LIB):"
              echo "------------------------------"
              echo " "
              echo "The mpi test code bin/minimal_mpi_test was not run"
              echo "successfully." 
              echo " " 
              echo "You may want to check the following:" 
              echo "-- Are you sure your mpi demons have been started?"
              echo "   E.g. under lam you have to use the lamboot "
              echo "   command to get mpi up and running...."
              echo "-- Is the mpi run command you specified via "
              echo "   the --with-mpi-self-tests-variablenp flag in the"
              echo "   configure options valid?"
              echo " "
              echo "I will continue regardless but you shouldn't be"
              echo "surprised if large numbers of mpi self-test fail..."
              echo " "
              echo "================================================="
              echo " "
            else
              echo " " 
              echo " The file $OUTFILE is correct, good stuff!"
            fi

            rm -f $OUTFILE
          done
        fi
      else
        echo " "
        echo "================================================="
        echo " "
        echo "WARNING (ISSUED BY OOMPH-LIB):"
        echo "------------------------------"
        echo " "
        echo "Sorry I got myself confused when parsing " 
        echo "the configure options and could not find"
        echo "the specification of the c++ mpi compiler"
        echo "via the CXX flag. Please make sure this is"
        echo "specified in configure/configure_optinios/current"
        echo "in the form: CXX=mpic++, say, (no quotes, no spaces)."
        echo " " 
        echo "NOTE: This is not necessarily a problem, but will"
        echo "keep me from checking if mpi is up and running"
        echo "before starting the self-tests."
        echo " "
        echo "================================================="
        echo " "        
      fi

    else
      echo "I cannot find OOMPHNP in the run command."
      echo "OOMPHNP will be replaced by the number of processors to run the code on."
      echo "For example, if your mpi run comamnd for two processors is "
      echo " "
      echo "  mpirun -np 2"
      echo " "
      echo "then the mpi run command for self tests on a variable number of"
      echo "processors must be"
      echo " "
      echo "  mpirun -np OOMPHNP"
      echo " "
      echo "I will continue regardless but you shouldn't be"
      echo "surprised if large numbers of mpi self-test fail..."
    fi # OOMPHNP_CHECK
    echo " "
  else
    echo " "
    echo "No mpi self-tests were requested"
    echo " "
  fi  

fi



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

# Now do the actual configure:
#------------------------------
# - prefix sets directory where lib and install directories get placed:
# - CXX =  C++ compiler (defaults to gcc/g++)
# - CC =   C compiler (defaults to gcc)
# - F77 =  F77 compiler (defaults to gcc/g77)
# 
# Options: "--enable-MPI" includes all MPI sources into the build
#


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




# Size after build
#-----------------
if test "$doc_size" = "y"; then 
    echo " "
    echo "Assessing size of the distribution after "
    echo "installation -- this can take a little while..."
    echo " " 
    echo "<TD align=\"center\">" >> size.html
    echo `du -s -h .` >> size.html
    echo "</TD>" >> size.html
    echo " "
    echo "...done."
    echo " " 
fi


# Make the demo codes and run the ones that are listed in the TESTS
#------------------------------------------------------------------
# declaration in the "Makefile.am"s
#-----------------------------------
if test "$want_self_tests" = "y" -o "$reply" = "Y" ; then 
   
  # We have to turn off "crash on errors" here because we don't want everything
  # to stop if a single self test fails.
  set -o errexit false

  # Use -k on make check so that make will keep going even if a demo driver
  # fails to build.
  echo " "
  echo "Running check to build/self-test the demo codes."
  echo "y" | make $make_options check -k
  echo "Done self test"

  # and now turn it back on (just in case...)
  set -o errexit

  # Size after self-tests
  #----------------------
  if test "$doc_size" = "y"; then 
      echo " "
      echo "Assessing size of the distribution after "
      echo "self tests -- this can take a little while... "
      echo " " 
      echo "<TD align=\"center\">" >> size.html
      echo `du -s -h .` >> size.html
      echo "</TD>" >> size.html
      echo " "
      echo "...done."
      echo " " 
  fi

else
  echo " "
  echo "The build process is complete. You may now "
  echo "initiate the self-tests by typing \"make check -k\". "
  echo " "


  # Size after self-tests
  #----------------------
  if test "$doc_size" = "y"; then 
      echo "<TD align=\"center\">" >> size.html
      echo "n/a" >> size.html
      echo "</TD>" >> size.html
  fi

fi


#Finish off size file
#--------------------
if test "$doc_size" = "y"; then 
    echo "</TR>" >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo " " >> size.html
    echo "</TABLE>" >> size.html
    echo "</html>" >> size.html
fi


echo " "
echo "autogen.sh has finished! If you can't spot any error messages" 
echo "above this, oomph-lib should now be ready to use... " 
echo " " 
echo "If you encounter any problems, please study the installation" 
echo "instructions and the FAQ before contacting the developers. " 
echo " " 

