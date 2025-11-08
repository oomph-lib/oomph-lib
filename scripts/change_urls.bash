#! /bin/bash

#-----------------------------------------------------------
# Search through all files and replace one url with another
# Ft. Puneet magic to escape awkward characters...
#-----------------------------------------------------------

orig_url=http://oomph-lib.maths.man.ac.uk/oomph-lib_external_distfiles
orig_url=http://www.maths.manchester.ac.uk/~oomphlib/oomph-lib_external_distfiles

new_url=https://personal.maths.manchester.ac.uk/oomphlib/oomph-lib_external_distfiles


# Helper function to handle the typically-messy escaping of characters
# Puneet's magic!
escape_sed_args() {
    echo "$1" | sed 's/[\/&]/\\&/g'
}

# Escape characters
orig_url=$(escape_sed_args "$orig_url")
new_url=$(escape_sed_args "$new_url")

sed_string="s/"$orig_url"/"$new_url"/g"
echo $sed_string



echo "find . \( ! -name .command.bash -a ! -name change_urls.bash \) -type f -exec grep -m 1 -H '"$orig_url"' {} \;" > .command.bash
chmod a+x .command.bash
full_list=`./.command.bash`
rm -f .command.bash

for full_line in `echo $full_list`; do
    file=`echo $full_line | awk 'BEGIN{FS=":"}{print $1}'`
    echo "doing "$file 
    cp $file $file.back
    sed $sed_string $file.back > $file
    rm -f $file.back
done
