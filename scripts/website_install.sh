#!/usr/bin/env bash 

build_dir="build/website"

# Ensure build directory exists
mkdir -p $build_dir

# Copy across all the stuff we need
cp -a bin config demo_drivers doc external_src self_test src user_drivers user_src -t $build_dir

# Delete some problematic sfiles
rm -f $build_dir/bin/fig2poly
rm -f $build_dir/bin/create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh

# Remove symlinks to files that don't get made for some reason, need to come back to this
rm -f $build_dir/doc/poisson/one_d_poisson/latex/refman.pdf
rm -f $build_dir/doc/unsteady_heat/two_d_unsteady_heat/latex/refman.pdf

# Move gitignore
cp website.gitignore $build_dir/.gitignore

# Make .nojekyll file to stop jekyll compiling
touch $build_dir/.nojekyll
