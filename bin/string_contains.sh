#!/bin/sh

#---------------------------------------------------------------------------
# Taken from http://stackoverflow.com/questions/2829613/how-do-you-tell-if-a-string-contains-another-string-in-unix-shell-scripting
#
# A unix shell string contains function.
# This uses POSIX substring parameter expansion, 
# so it works in bash, dash, ksh...
#
# contains string, substring
#
# Returns 0 if the specified string contains the specified substring,
# otherwise returns 1.
#
# Example:
#
# # First source this file:
# . $OOMPH_ROOT_DIR/bin/string_contains.sh
#
# # Run the contains function:
# contains "my_long_string" "long" 
# 
# # We get the return value via $?
# constainsRetVal=$?
#
# # Now we can use this value, eg.
#
# if [ "$containsRetVal" -eq "0" ]
# then
#   echo "Long string contains the substring."
# else
#   echo "Substring is not in long string."
# fi
# 
#---------------------------------------------------------------------------

contains() {
    string="$1"
    substring="$2"
    if test "${string#*$substring}" != "$string"
    then
        return 0    # $substring is in $string
    else
        return 1    # $substring is not in $string
    fi
}
