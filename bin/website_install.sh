#!/usr/bin/env bash 

build_dir="build/website"

# Ensure build directory exists
mkdir -p $build_dir

# Copy across all the stuff we need
cp -a bin config demo_drivers doc external_src self_test src user_src -t $build_dir

# Move gitignore
cp website.gitignore $build_dir/.gitignore

# Make .nojekyll file to stop jekyll compiling
touch $build_dir/.nojekyll
