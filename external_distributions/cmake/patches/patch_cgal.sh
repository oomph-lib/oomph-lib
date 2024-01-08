#!/bin/bash

CGAL_ROOT_DIR=$1

echo "Patching files in CGAL_ROOT_DIR: ${CGAL_ROOT_DIR}"
echo ""
if [ ! -d "${CGAL_ROOT_DIR}" ]; then
    echo "Argh! This isn't a directory!"
    echo ""
    exit 1
fi

# Need to change all instances of 'BUILD_TESTING' --> 'CGAL_BUILD_TESTING'
if grep -qR '\bBUILD_TESTING\b' "${CGAL_ROOT_DIR}"; then
    echo "Found instances of 'BUILD_TESTING'! Will change them to 'CGAL_BUILD_TESTING'"
    echo ""

    # Update 'BUILD_TESTING' --> 'CGAL_BUILD_TESTING'
    if [ "$(uname)" == "Darwin" ]; then
        LC_ALL=C find "${CGAL_ROOT_DIR}" -type f -exec sed -i '' "s|BUILD_TESTING|CGAL_BUILD_TESTING|g" {} \;
    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
        LC_ALL=C find "${CGAL_ROOT_DIR}" -type f -exec sed -i "s|BUILD_TESTING|CGAL_BUILD_TESTING|g" {} \;
    fi

    # Make sure it worked
    if grep -qR '\bBUILD_TESTING\b' "${CGAL_ROOT_DIR}"; then
        echo "It looks like you failed to change all occurrences of BUILD_TESTING!"
        echo ""
        grep -R '\bBUILD_TESTING\b' "${CGAL_ROOT_DIR}"
    fi
fi

OLD_STRING='option(CGAL_BUILD_TESTING "Build the testing tree." OFF)'
NEW_STRING='option(CGAL_BUILD_TESTING "Build the testing tree." ${BUILD_TESTING})'
ESCAPED_OLD_STRING=$(printf '%s\n' "${OLD_STRING}" | sed -e 's/[\/&]/\\&/g')
ESCAPED_NEW_STRING=$(printf '%s\n' "${NEW_STRING}" | sed -e 's/[\/&]/\\&/g')

echo "Replacing '${ESCAPED_OLD_STRING}' with '${ESCAPED_NEW_STRING}'."
echo ""
sed -i '' "s/${ESCAPED_OLD_STRING}/${ESCAPED_NEW_STRING}/g" "${CGAL_ROOT_DIR}/CMakeLists.txt"

echo "Yay! Looks like it all worked!"
echo ""
