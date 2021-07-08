START
{
 print "// This file was generated automatically during the make process"
 print "// and it will be remade automatically"
}
{
 # Loop over all header file names
 for (i=2;i<=NF;i++) printf "#include<"$1"/"$i"> \n"
} 


