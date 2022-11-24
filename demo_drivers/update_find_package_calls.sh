#!/bin/bash

OOMPH_ROOT=$(pwd)/../
FILES=$(rg -l 'find_package\(oomphlib REQUIRED\)')

# find_package(oomphlib CONFIG REQUIRED PATHS "../install/")

for FILE in $FILES; do
    echo "${FILE}"
    PARENT_DIR=$(dirname ${FILE})
    RELATIVE_INSTALL_PATH=$(grealpath --relative-to="${PARENT_DIR}" "${OOMPH_ROOT}/install")
    sd -s "find_package(oomphlib REQUIRED)" "find_package(oomphlib CONFIG REQUIRED PATHS \"${RELATIVE_INSTALL_PATH}\")" "${FILE}"
done
