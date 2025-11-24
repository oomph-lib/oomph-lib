#! /bin/bash

CGAL_ROOT_DIR=$1
LOG_DIR=$2
BUILD_OPTIONS="${@:3}"

echo "Entering CGAL root directory: ${CGAL_ROOT_DIR}"
cd ${CGAL_ROOT_DIR}

echo "Will log to directory: ${LOG_DIR}"

if [ -e oomph_test ]; then
    echo "Deleting oomph_test directory (resumably there from "
    echo "previous installation without wiping sources)"
    rm -rf oomph_test
fi

driver_src_dir=Spatial_searching/examples/Spatial_searching
scripts_src_dir=Scripts/scripts

mkdir oomph_test
cp ${driver_src_dir}/nearest_neighbor_searching.cpp oomph_test
cd oomph_test

if [ ! -e ../${scripts_src_dir}/cgal_create_CMakeLists ]; then
    echo "ERROR: The CGAL script "
    echo " "
    echo "   ${scripts_src_dir}/cgal_create_CMakeLists"
    echo " "
    echo "that we use to test our CGAL installation doesn't exist (any more?)"
    exit 1
fi

# Create CMakeLists file
../${scripts_src_dir}/cgal_create_CMakeLists -s nearest_neighbor_searching &>"${LOG_DIR}/build.log"

# Configure
echo "Building CGAL test with compiler_spec_string: "$compiler_spec_string
echo "Building CGAL test with build_opts: ${BUILD_OPTIONS}"
echo $compiler_spec_string" cmake -DCGAL_DIR=.. ${BUILD_OPTIONS} . 2>&1 >> ${LOG_DIR}/build.log" >.full_build_file
source .full_build_file

# Build
echo " "
make 2>&1 >>"${LOG_DIR}/build.log"
if [ ! -e ./nearest_neighbor_searching ]; then
    echo "ERROR: CGAL test code "
    echo " "
    echo "  "$(pwd)"/nearest_neighbor_searching"
    echo " "
    echo "which was copied from "
    echo" "
    echo "    ${driver_src_dir}"
    echo " "
    echo "failed to build. Check ${LOG_DIR}/build.log!"
    exit 1
else
    echo "Yay! CGAL test code "
    echo " "
    echo "  "$(pwd)"/nearest_neighbor_searching"
    echo " "
    echo "which was copied from ${driver_src_dir}"
    echo "was built!"
    output=$(./nearest_neighbor_searching)
    echo " "
    if [ "$output" != "0 0 0" ]; then
        echo "ERROR: CGAL failed: Output should be: \"0 0 0\" but is \"$output\""
        exit 1
    else
        echo "Yay: CGAL test passed: Output should be: \"0 0 0\" and is \"$output\""
    fi
    echo " "
fi

exit 0
