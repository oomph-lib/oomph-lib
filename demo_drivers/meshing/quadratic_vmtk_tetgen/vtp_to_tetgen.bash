#! /bin/bash

#===============================================================
# Shell script to generate tetgen fluid and solid meshes from 
# vmtk *.vtp files.
#
# Use edgelengthfactor of 0.6 and quality mesh for fine
#     edgelengthfactor of 0.8 and non-quality mesh for coarse
#===============================================================

# Stem -- adjust; Specifying "iliac" will work with iliac.vtp
#------------------------------------------------------------
stem=iliac

#stem=iliac_with_extensions

# wall thickness -- adjust
#-------------------------
wall_thickness=1.0

# Quality mesh? -- adjust
#------------------------
quality_mesh=0
    




# Build code that converts vmtk-generated fluid xda file to
# poly files for tetgen based mesh generation for fluid and solid
make xda_to_poly_fsi

# Create/empty working directory
if [ ! -e vmtk_working_dir ]; then
    mkdir vmtk_working_dir
fi
rm -f vmtk_working_dir/*


# Move important files across into working directory so we don't
# clutter up filespace
echo "========================================================================"
echo "I'm working in "`pwd` 
echo "========================================================================"
cp vmtk_files/*vtp vmtk_working_dir
cd vmtk_working_dir
ln -sf ../xda_to_poly_fsi

# Read the vtp file and turn into vtu: 0.6 for fine 0.8 for coarse
#-----------------------------------------------------------------
vmtksurfacereader -ifile `echo $stem`.vtp  \
--pipe vmtkcenterlines -endpoints 1  -seedselector carotidprofiles  \
--pipe vmtkdistancetocenterlines -useradius 1 -centerlineradius 1 \
--pipe vmtkmeshgenerator -ofile `echo $stem`.vtu \
-edgelengthfactor 0.8 \
-elementsizemode edgelengtharray -edgelengtharray DistanceToCenterlines


# NOTES:
#-------
# USE: -seedselector openprofiles TO IDENTIFY INLET/OUTLET IDS
#      INTERACTIVELY
#
#       --pipe vmtkmeshgenerator -ofile `echo $stem`.vtu -edgelengthfactor 0.6 \
#       TO CHANGE EDGELENGTH FACTOR FOR FINER/COARSER MESHES



# Convert to quadratic elements -- Note: Have to specify original vtp file too!
#------------------------------------------------------------------------------
mv `echo $stem`.vtu junk.vtu
vmtklineartoquadratic -ifile  junk.vtu \
   -ofile `echo $stem`.vtu \
   -rfile `echo $stem`.vtp  -entityidsarray CellEntityIds 
rm junk.vtu

# Convert to xda format
#----------------------
vmtkmeshwriter \
    -ifile `echo $stem`.vtu \
    -entityidsarray CellEntityIds \
    -ofile `echo $stem`.xda
    

# Create tetgen files for fluid and solid domains
#------------------------------------------------
rm -f .vtp_to_tetgen_replies.txt

# Stem of filename
echo $stem > .vtp_to_tetgen_replies.txt

# Wall thickness
echo $wall_thickness >> .vtp_to_tetgen_replies.txt

# Create meshes for fsi problem, i.e. label each
# interface facet as a separate boundary
echo "y" >> .vtp_to_tetgen_replies.txt
echo " " >> .vtp_to_tetgen_replies.txt

# Run conversion code
./xda_to_poly_fsi <  .vtp_to_tetgen_replies.txt

# Clean up 
rm -f .vtp_to_tetgen_replies.txt
    

# Build tetgen meshes
fluid_quality_flag=""
solid_quality_flag=""
if [ $quality_mesh -eq 1 ]; then
    fluid_quality_flag="-q1.2"
    solid_quality_flag="-q"
fi

tetgen `echo $fluid_quality_flag` fluid_`echo $stem`.poly
tetgen `echo $solid_quality_flag` solid_`echo $stem`.poly

cp solid_`echo $stem`.poly solid.poly
cp solid_`echo $stem`.1.ele  solid.1.ele
cp solid_`echo $stem`.1.node solid.1.node
cp solid_`echo $stem`.1.face solid.1.face

cp fluid_`echo $stem`.poly fluid.poly
cp fluid_`echo $stem`.1.ele  fluid.1.ele
cp fluid_`echo $stem`.1.node fluid.1.node
cp fluid_`echo $stem`.1.face fluid.1.face

cp `echo $stem`_boundary_enumeration.dat boundary_enumeration.dat
cp `echo $stem`_quadratic_fsi_boundary.dat quadratic_fsi_boundary.dat
cp `echo $stem`_quadratic_outer_solid_boundary.dat quadratic_outer_solid_boundary.dat

cd ..

