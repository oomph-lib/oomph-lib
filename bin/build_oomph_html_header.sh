#! /bin/sh

#-------------------------------------------------------------------
# Shell script to build html header file for oomph-lib documentation.
# The root directory must be passed as the first argument.
#-------------------------------------------------------------------

# Test that we do have (only one) command line argument
if (test $# -ne 1); then
 echo "Path to oomph-lib root directory must be the only argument\
 to build_oomph_html_header.sh"
 exit 1
fi

#-------------------------------------------------------------------
# Build/execute sed script to replace the #### placeholder in the template
# file by relative path to the oomph-lib root directory:
#-------------------------------------------------------------------
echo -n 'sed "s/############/'  > tmp.sed
echo -n $1 | sed 's#/#\\\/#g' >> tmp.sed
echo -n '/g"'  >> tmp.sed
chmod a+x tmp.sed
./tmp.sed < $1/doc/oomph-lib_header.html.template > oomph-lib_header.html
rm tmp.sed




