#!/usr/bin/env bash

# A tiny helper that downloads and installs a recent CMake (>= 3.24).
# It follows the exact steps documented in README.md.
#
# Usage:
#   ./scripts/install_cmake.sh                # install 3.24.4 under ~/.cmake-3.24.4
#   CMAKE_VERSION=3.25.3 ./scripts/install_cmake.sh   # choose a different version
#   CMAKE_DEST_DIR=$HOME/opt/cmake ./scripts/install_cmake.sh  # choose destination
#
# After it finishes, it prints an export line you can add to ~/.bashrc or ~/.zshrc.

set -euo pipefail

CMAKE_VERSION="${CMAKE_VERSION:-3.24.4}"
CMAKE_DEST_DIR="${CMAKE_DEST_DIR:-$HOME/.local/.cmake-${CMAKE_VERSION}}"

if [ -d ${CMAKE_DEST_DIR} ]; then
    echo "Desired installation directory for CMake, i.e."
    echo
    echo "      ${CMAKE_DEST_DIR}"
    echo
    echo "already exists!"
    exit 1
fi

# Resolve OS
OS="$(uname -s)"
TMP="$(mktemp)"

# Helper message
function final_instructions() {
    echo
    echo "CMake ${CMAKE_VERSION} installed to:"
    echo
    echo "      ${CMAKE_DEST_DIR}"
    echo
    echo "Add the following line to your shell start-up file to make it permanent:"
    echo
    echo "      export PATH=${1}:\$PATH"
}

case "$OS" in
Linux*)
    echo "Installing CMake ${CMAKE_VERSION} for Linux ..."
    mkdir -p "${CMAKE_DEST_DIR}"
    URL="https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh"
    curl -L "${URL}" -o "${TMP}"
    chmod +x "${TMP}"
    bash "${TMP}" --prefix="${CMAKE_DEST_DIR}" --exclude-subdir
    rm -f "${TMP}"
    final_instructions "${CMAKE_DEST_DIR}/bin"
    ;;

Darwin*)
    echo "Installing CMake ${CMAKE_VERSION} for macOS ..."
    URL="https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-macos-universal.tar.gz"
    curl -sL "${URL}" -o "${TMP}"
    tar -xzf "${TMP}"
    mv "cmake-${CMAKE_VERSION}-macos-universal" "${CMAKE_DEST_DIR}"
    rm -f "${TMP}"
    final_instructions "${CMAKE_DEST_DIR}/CMake.app/Contents/bin"
    ;;

*)
    echo "Unsupported OS: ${OS} — please install CMake manually." >&2
    exit 1
    ;;
esac
