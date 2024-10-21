# This is discription on how to use gmsh in oomphlib
- use gmsh software to generate .msh file (please visit https://gmsh.info/)
- oomphlib gmsh reader renumber boundaries starting from zero e.g. if you have in .msh file 
$PhysicalNames
6
2 27 "XMIN"
2 28 "XMAX"
2 33 "ZMAX"
2 34 "ZMIN"
2 35 "YMIN"
2 36 "YMAX"
$EndPhysicalNames
then XMIN boundary value = 0, XMAX = 1, and ZMAX = 2 so on.

# Future improvements:
- the current reader does not support reading multiple regions( e.g. to assign different material density, consitiutive laws etc.)


# mesh_from_gmsh_solid.cc is a demo.
