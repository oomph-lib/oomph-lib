#! /bin/bash

#------------------------------------------------
# Shell script to build mesh of bifurcation
# from line potential
#------------------------------------------------

# Make the executables
make potential_for_bifurcation_mesh
make create_poly_for_bifurcation_from_tecplot

# Create potential in cubic box: This creates "potential.dat"
# hierher 
#./potential_for_bifurcation_mesh

# Plot "potential.dat" in tecplot and determine and extract iso-surface into 
# "potential_boundary.dat"
tecplot -mesa -b -p create_bifurcation_from_potential.mcr

# Display
xanim bifurcation_from_potential.avi

# Now convert extracted iso-surface from "potential_boundary.dat" into 
# "fsi_potential_bifurcation_fluid.poly"
./create_poly_for_bifurcation_from_tecplot

# Display poly file
#tetview fsi_potential_bifurcation_fluid.poly
#tetview fsi_potential_bifurcation_solid.poly

# Mesh it!
tetgen fsi_potential_bifurcation_fluid.poly
tetgen fsi_potential_bifurcation_solid.poly

# Display the mesh
#tetview fsi_potential_bifurcation_fluid.1.ele
#tetview fsi_potential_bifurcation_solid.1.ele

