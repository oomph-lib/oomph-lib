#! /bin/bash -f



echo "<?xml version=\"1.0\"?>" > combined.pvd
echo "<VTKFile type=\"Collection\" version=\"0.1\">" >> combined.pvd
echo "<Collection>" >> combined.pvd

count=0
for file in `echo $@`; do
    echo $file
    #oomph-convert $file
    #file_base_only=`basename $file .dat`
    #file_dir=`dirname $file`
    #file_base=`echo $file_dir"/"$file_base_only`
    #echo $file_base
    #echo "<DataSet timestep= "$count" part=\" 0 \" file=\""$file_base".vtu/>\"" >> combined.pvd
echo "<DataSet timestep= \" "$count" \" part=\" 0 \" file=\""$file"\"/>" >> combined.pvd
    let count=$count+1
done

echo "</Collection>" >> combined.pvd
echo "</VTKFile>" >> combined.pvd
