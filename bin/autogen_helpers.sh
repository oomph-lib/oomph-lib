
# A few helper functions for autogen or other bash scripts

# Read an option using readline
OptionRead()
{
    read -e Opt
    echo $Opt
}


# Convert a string passed as argument to return true/false (ie a bool),
# return value 255 indicates error.
YNToBool()
{
    if [[ $1 == "y" || $1 == "Y" ]]; then
        return $(true) # 0 == true in bash
    elif [[ $1 == "n" || $1 == "N" ]]; then
        return $(false) # 1 is false in bash
    else
        echo "I don't understand \"$1\", y/n only please." 1>&2
        return 255
    fi
}


# Read a y/n answer from input, if the answer is not y/n then repeat until
# it is. Turns out this can be nasty when combined with set -o errexit: you
# can only call this inside an if statement otherwise bash thinks any false
# return value is an error!
YesNoRead()
{
    prompt="$1"
    default="$2"

    # Read an answer
    printf "%s y/n [default %s]: " "$1" "$2"
    read Opt

    # Convert to a bool
    bool=""
    if [[ $Opt == "" ]]; then
        YNToBool $default
        bool=$?
    else
        YNToBool $Opt
        bool=$?
    fi

    # If we didn't recognise it then try again
    if [[ $bool != 1 && $bool != 0 ]]; then
        YesNoRead $prompt $default
        return $?
    fi

    return $bool
}


# This little function takes the input, removes anything following a #
# deletes blanks lines and then replaces all newlines by spaces
ProcessOptionsFile()
{
    sed < $1 -e 's/#.*$//' -e '/^[ \t]*$/d' -e 's/\n/ /g'
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
    }' $(echo $1)
}


# This little function echo's the usage information
EchoUsage()
{
    echo "  h)"
    echo "      Show this usage info"
    echo "  "
    echo "  r)"
    echo "      Regenerate the configuration files before compliling,"
    echo "      this is needed if you have added new user drivers."
    echo "  "
    echo "  d)"
    echo "      Set the oomph-lib root directory (default: use pwd)."
    echo "  "
    echo "  c)"
    echo "      Set the configure options file to use (default: current)."
    echo "  "
    echo "  j)"
    echo "      Set the number of cores to use for compilation (default: 1)."
    echo "  "
    echo "  b)"
    echo "      Set the directory to put object files etc. in (default pwd/build)."
}
