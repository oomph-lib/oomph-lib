#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# patch_cgal.sh  <CGAL_ROOT_DIR>
#
#  * Renames BUILD_TESTING --> CGAL_BUILD_TESTING everywhere.
#  * Adjusts CGAL_BUILD_TESTING option default.
#  * Rewrites Boost lookup so it uses CONFIG mode:
#       find_package( Boost 1.72 REQUIRED )
#          --> find_package( Boost 1.72 CONFIG REQUIRED )
# -----------------------------------------------------------------------------

CGAL_ROOT_DIR=$1

echo "Patching files in CGAL_ROOT_DIR: ${CGAL_ROOT_DIR}"
echo ""
if [ ! -d "${CGAL_ROOT_DIR}" ]; then
    echo "Argh! This isn't a directory!"
    echo ""
    exit 1
fi

# -----------------------------------------------------------------------------
# BUILD_TESTING --> CGAL_BUILD_TESTING
# -----------------------------------------------------------------------------
echo "Scanning for BUILD_TESTING"
if grep -qR '\bBUILD_TESTING\b' "${CGAL_ROOT_DIR}"; then
    echo ""
    echo "Found instances of 'BUILD_TESTING'! Will change them to 'CGAL_BUILD_TESTING'"

    # Update 'BUILD_TESTING' --> 'CGAL_BUILD_TESTING'
    if [ "$(uname)" == "Darwin" ]; then
        LC_ALL=C find "${CGAL_ROOT_DIR}" -type f -exec sed -i '' "s|BUILD_TESTING|CGAL_BUILD_TESTING|g" {} \;
    else
        LC_ALL=C find "${CGAL_ROOT_DIR}" -type f -exec sed -i "s|BUILD_TESTING|CGAL_BUILD_TESTING|g" {} \;
    fi

    # Make sure it worked
    echo ""
    if grep -qR '\bBUILD_TESTING\b' "${CGAL_ROOT_DIR}"; then
        echo "It looks like you failed to change all occurrences of BUILD_TESTING!"
        grep -R '\bBUILD_TESTING\b' "${CGAL_ROOT_DIR}"
    else
        echo "No BUILD_TESTING tokens found."
    fi
fi

# -----------------------------------------------------------------------------
# Force Boost to be found via CONFIG mode (CMake 1.72)
# File: cmake/modules/CGAL_SetupBoost.cmake
# -----------------------------------------------------------------------------
BOOST_FILE="${CGAL_ROOT_DIR}/cmake/modules/CGAL_SetupBoost.cmake"

echo ""
if [[ -f ${BOOST_FILE} ]]; then
    echo "Patching Boost find_package() call in ${BOOST_FILE}"
    if [ "$(uname)" == "Darwin" ]; then
        sed -E -i '' "s|find_package\( Boost 1\.72 REQUIRED \)|find_package( Boost 1.72 CONFIG REQUIRED )|g" "${BOOST_FILE}"
    else
        sed -E -i "s|find_package\( Boost 1\.72 REQUIRED \)|find_package( Boost 1.72 CONFIG REQUIRED )|g" "${BOOST_FILE}"
    fi
else
    echo "Expected Boost cmake module not found: ${BOOST_FILE}"
fi

echo ""
echo "Yay! Looks like it all worked!"
echo ""
