#!/bin/bash
 
# Move to the demo drivers directory; edit (or get prompt from user)
# if this is not run from oomph-lib home directory
demo_drivers=`pwd`"/demo_drivers"
echo " "
echo "Running all self-tests that involve triangle meshes in: "
echo " " 
echo "    "$demo_drivers
echo " "
echo "[Note: This has to be the absolute address!] "
echo " "
cd $demo_drivers
 
# Name of the file
file="file_names.txt"
 
# Get the file names which use a triangle-based mesh
grep -R "TriangleMesh<" . | awk -F '[:]' '{print $1}' > ../$file
 
# Move back up
cd ..
 
# Remove the lines saying Binary file
grep -v "Binary file" $file > temp.txt
mv temp.txt $file
 
# Remove the last word from each line of the file. I.e. the .cc file name
awk -F '[.]' '{print $2}' $file > temp.txt
mv temp.txt $file
 
# Command to output text file
filelines=`cat $file`
 
# Loop over the lines in the file
for line in $filelines ; do
    # Grab the last word
    word=`echo $line | awk -F '[/]' 'NF>1{print $NF}'`
 
    # Remove this word from the end of the line and put it in a file
    path_old=${line%$word}
 
    # Update this
    echo ${path_old#"/"} >> temp.txt
done
 
# Replace the old file
mv temp.txt $file
 
# Take out duplicate lines
awk '!a[$0]++' $file > temp.txt
 
# Move it back
mv temp.txt $file
 
# Command to output text file
filelines=`cat $file`

# Now we've got the paths we need to run there and make sure they work
for line in $filelines; do

    # Jump to demo drivers
    cd $demo_drivers
     
    # Jump to file directory
    cd "$line";

    echo "Running in: "`pwd`
done

# Now we've got the paths we need to run there and make sure they work
for line in $filelines; do

    # Jump to demo drivers
    cd $demo_drivers
     
    # Jump to file directory
    cd "$line";

    # Make check
    make check -k
done
