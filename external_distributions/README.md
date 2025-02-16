<div align="center">
  <a href="http://oomph-lib.maths.man.ac.uk">
    <img alt="oomph-lib logo" src="./assets/oomph_logo.png">
  </a>
</div>

<div align="center">
  <a href="./LICENCE">
    <img alt="License: LGPL v2.1" src="https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg">
  </a>
</div>

<br>

<div align="center">
    <table>
    <tr>
        <th>Test platform</th>
        <th><a href="../../../tree/feature-shared-cmake"><code>feature-shared-cmake</code></a></th>
    </tr>
    <tr>
        <td>Ubuntu</td>
        <td>
            <a href="../../../actions/workflows/test-third-party-libs-on-ubuntu.yaml">
                <img alt="Ubuntu" src="../../../actions/workflows/test-third-party-libs-on-ubuntu.yaml/badge.svg?branch=feature-shared-cmake" style="vertical-align: middle">
            </a>
        </td>
    </tr>
    <tr>
        <td>macOS</td>
        <td>
            <a href="../../../actions/workflows/test-third-party-libs-on-macos.yaml">
                <img alt="macOS" src="../../../actions/workflows/test-third-party-libs-on-macos.yaml/badge.svg?branch=feature-shared-cmake" style="vertical-align: middle">
            </a>
        </td>
    </tr>
    </table>
</div>

<!-- Use <h2> tags to omit heading from table of contents -->
<h2>Table of contents</h2>

- [Description](#description)
- [Prerequisites](#prerequisites)
- [Library versions](#library-versions)
- [Example](#example)
- [Build options](#build-options)
  - [Remark on `OOMPH_USE_<LIBRARY>_FROM` variables](#remark-on-oomph_use_library_from-variables)
  - [Extended example](#extended-example)
  - [macOS support](#macos-support)

## Description

This project builds and installs the optional third-party libraries used by [`oomph-lib`](https://github.com/oomph-lib/oomph-lib).

## Prerequisites

Tool    | Version
--------|--------
`CMake` | 3.24

If you cannot obtain a recent enough version of CMake via your favourite package manager, you will need to build it from source. Worry not however, this is a straightforward task. For details on how to do this, see [the instructions here](https://github.com/puneetmatharu/oomph-lib/tree/feature-shared-cmake#building-cmake).

## Library versions

The table below contains the version of each library that you can install with this repository. If you wish to provide your own version of OpenBLAS or Boost, make sure it is the right version.

Library    | Version
-----------|--------
`OpenBLAS` | 0.3.25
`Boost`    | 1.83.0
`CGAL`     | 5.6
`MUMPS`    | 5.6.2
`HYPRE`    | 2.32.0
`Trilinos` | 16.0.0

## Example

To build all libraries without MPI support, simply run

```bash
>>> cmake -G Ninja -B build
>>> cmake --build build
```

Note that we do not need an install step; this is because we only build/install the third-party libraries in this project, which are installed at build-time.

## Build options

The table below contains the flags that you can pass to the `cmake` command to control the build of the third-party libraries. Arguments to the `cmake` command must adhere to the format `-D<FLAG_1>=<VALUE_1>`. For examples on how to do this, see [Extended example](#extended-example).

Option                                      | Description                                                | Default
--------------------------------------------|------------------------------------------------------------|----------------------------------
`OOMPH_ENABLE_MPI`                          | *Enable MPI support?*                                      | `OFF`
`OOMPH_BUILD_OPENBLAS`                      | *Build OpenBLAS?*                                          | `ON`
`OOMPH_BUILD_SUPERLU`                       | *Build SuperLU?*                                           | `ON`
`OOMPH_BUILD_SUPERLU_DIST`                  | *Build SuperLU DIST?*                                      | `ON` if MPI is enabled else `OFF`
`OOMPH_BUILD_CGAL`                          | *Build CGAL (with dep.  Boost)?*                           | `ON`
`OOMPH_BUILD_MUMPS`                         | *Build MUMPS?*                                             | `ON` if MPI is enabled else `OFF`
`OOMPH_BUILD_HYPRE`                         | *Build Hypre?*                                             | `ON`
`OOMPH_BUILD_TRILINOS`                      | *Build Trilinos?*                                          | `ON`
`OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTING` | *Disable testing when building the third-party libraries?* | `OFF`
`OOMPH_THIRD_PARTY_INSTALL_DIR`             | *Base installation directory for third-party libraries.*   | `<project_root>/install/`
`OOMPH_USE_OPENBLAS_FROM`                   | *The path to a preinstalled version of OpenBLAS.*          | `""`
`OOMPH_USE_BOOST_FROM`                      | *The path to a preinstalled version of Boost.*             | `""`

### Remark on `OOMPH_USE_<LIBRARY>_FROM` variables

The arguments to the `OOMPH_USE_<LIBRARY>_FROM` flags must be a folder containing a `lib/` and `include/` folder (and possibly even a `bin/` folder) for that library. As an example, consider the Boost library. Below we show how we locate the Boost installation of Homebrew

```bash
>>> brew --prefix boost
/opt/homebrew/opt/boost

>>> ls -1 /opt/homebrew/opt/boost@1.87/
include/
lib/
share/
INSTALL_RECEIPT.json
README.md
sbom.spdx.json
```

To use this installation of Boost in the build of the third-party libraries, you would call `cmake` like so

```bash
>>> cmake -G Ninja -DOOMPH_USE_BOOST_FROM=/opt/homebrew/opt/boost@1.87/ -B build
```

or, more simply

```bash
>>> cmake -G Ninja -DOOMPH_USE_BOOST_FROM=$(brew --prefix boost) -B build
```

### Extended example

**Example 1:** To build all of the third-party libraries without any testing (at your own peril!), run

```bash
>>> cmake -G Ninja -DOOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS=ON -B build
>>> cmake --build build
```

**Example 2:** To enable MPI support but only build CGAL and OpenBLAS without any testing whatsoever, run

```bash
>>> cmake -G Ninja -DOOMPH_ENABLE_MPI=ON -DOOMPH_BUILD_CGAL=ON -DOOMPH_BUILD_OPENBLAS=ON -DOOMPH_BUILD_MUMPS=OFF -DOOMPH_BUILD_HYPRE=OFF -DOOMPH_BUILD_TRILINOS=OFF -DOOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS=ON -B build
>>> cmake --build build
```

Since, by default, we want to build all of the third-party libraries, you can actually reduce the above commands to

```bash
>>> cmake -G Ninja -DOOMPH_ENABLE_MPI=ON -DOOMPH_BUILD_MUMPS=OFF -DOOMPH_BUILD_HYPRE=OFF -DOOMPH_BUILD_TRILINOS=OFF -DOOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS=ON -B build
>>> cmake --build build
```

**Example 3:** To build everything with your own copy of OpenBLAS and Boost, run

```bash
>>> cmake -G Ninja -DOOMPH_USE_OPENBLAS_FROM=<path-to-installed-openblas> -DOOMPH_USE_BOOST_FROM=<path-to-installed-boost> -B build
>>> cmake --build build
```

See [Remark on `OOMPH_USE_<LIBRARY>_FROM` variables](#remark-on-oomph_use_library_from-variables) for instructions on how to determine the correct value for `<path-to-installed-XXXX>`.

### macOS support

Currently the build of OpenBLAS does not work on macOS with the compilers detected by CMake, but it does work with `gcc`/`g++`. However, we don't recommend that you specify your own compiler; instead you should let CMake pick the compiler for you. For this reason, we recommend that you install OpenBLAS via your desired package manager (e.g. Homebrew or MacPorts). If, for example, you install these libraries with Homebrew, you can specify their paths to CMake like so

```bash
>>> cmake -G Ninja -DOOMPH_USE_OPENBLAS_FROM=$(brew --prefix openblas) -B build
>>> cmake --build build
```
