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
        <th><a href="../../../tree/main"><code>main</code></a></th>
        <th><a href="../../../tree/development"><code>development</code></a></th>
    </tr>
    <tr>
        <td>Ubuntu</td>
        <td>
            <a href="../../../actions/workflows/test-third-party-libs-on-ubuntu.yaml">
                <img alt="Ubuntu" src="../../../actions/workflows/test-third-party-libs-on-ubuntu.yaml/badge.svg?branch=main" style="vertical-align: middle">
            </a>
        </td>
        <td>
            <a href="../../../actions/workflows/test-third-party-libs-on-ubuntu.yaml">
                <img alt="Ubuntu" src="../../../actions/workflows/test-third-party-libs-on-ubuntu.yaml/badge.svg?branch=development" style="vertical-align: middle">
            </a>
        </td>
    </tr>
    <tr>
        <td>macOS</td>
        <td>
            <a href="../../../actions/workflows/test-third-party-libs-on-macos.yaml">
                <img alt="macOS" src="../../../actions/workflows/test-third-party-libs-on-macos.yaml/badge.svg?branch=main" style="vertical-align: middle">
            </a>
        </td>
        <td>
            <a href="../../../actions/workflows/test-third-party-libs-on-macos.yaml">
                <img alt="macOS" src="../../../actions/workflows/test-third-party-libs-on-macos.yaml/badge.svg?branch=development" style="vertical-align: middle">
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

If you cannot obtain a recent enough version of CMake via your favourite package manager, you will need to build it from source. Worry not however, this is a straightforward task. For details on how to do this, see [the instructions here](https://github.com/oomph-lib/oomph-lib/tree/main#building-cmake).

## Library versions

`oomph-lib` depends on and works with a number of third-party libraries. In the table below we list these libraries and the versions that we can install for you. By default, we build "everything we can"; see further below for how to customise the build process and how to select a serial or an MPI build. If you wish to provide your own version of these libraries, make sure you provide the right version.


| Library        | Required/optional | Built by default (serial build)? | Built by default (MPI build)?  | Version |
| ----           | ---               | -----                    | ---                    | ---     | 
| `OpenBLAS`     | required by `oomph-lib`              | Yes (though not on macOS; see below) | Yes (though not on macOS; see below) |  [0.3.25](https://github.com/OpenMathLib/OpenBLAS/tree/v0.3.29)      |
| `SuperLU`       | required by `oomph-lib`               | Yes | Yes | [v6.0.1](https://github.com/xiaoyeli/superlu/tree/v6.0.1) | 
| `METIS`        | required by `oomph-lib` (via `SuperLU`) | Yes | Yes | [commit `a6e6a2cfa92f93a3ee2971ebc9ddfc3b0b581ab2`](https://github.com/KarypisLab/METIS/tree/a6e6a2cfa92f93a3ee2971ebc9ddfc3b0b581ab2)  |              
`GKlib`          | required by `oomph-lib` (via `METIS`)  | Yes | Yes | [commit `6e7951358fd896e2abed7887196b6871aac9f2f8`](https://github.com/KarypisLab/GKlib/tree/6e7951358fd896e2abed7887196b6871aac9f2f8)    |
| `SuperLU_DIST` | required for `oomph-lib` MPI build                   | No | Yes | [v9.1.0](https://github.com/xiaoyeli/superlu_dist/tree/v9.1.0)  
| `ParMETIS`     | required for `oomph-lib` MPI build (via `SuperLU_DIST`)                  | No | Yes | [commit `83bb3d4f5b2af826d0683329cad1accc8d829de2`](https://github.com/puneetmatharu/ParMETIS/tree/83bb3d4f5b2af826d0683329cad1accc8d829de2) | 
| `CGAL`         | optional, highly recommended                        | Yes | Yes | [6.0.1](https://github.com/CGAL/cgal/tree/v6.0.1)                                                                 |
| `Boost`        | required by `CGAL`                          | Yes | Yes | [1.83.0](https://github.com/boostorg/boost/tree/boost-1.83.0)                                                                             |
| `MUMPS`        | optional                                  | Yes | Yes | [5.6.2](https://github.com/puneetmatharu/mumps/tree/v5.6.2.5)                                                                             |
| `HYPRE`        | optional                                  | Yes | Yes | [2.32.0](https://github.com/hypre-space/hypre/tree/v2.32.0)                                                                               |
| `Trilinos`     | optional                                  | Yes | Yes | [16.0.0](https://github.com/trilinos/Trilinos/tree/trilinos-release-16-0-0)                                                               |


## Example

To build all libraries without MPI support, simply run

```bash
>>> cmake -G Ninja -B build
>>> cmake --build build
```
To build all libraries with MPI support (assuming your computer has MPI installed, of course), do
```bash
>>> cmake -G Ninja -DOOMPH_ENABLE_MPI=ON -B build
>>> cmake --build build
```
Note that we do not need an install step; this is because we only build/install the third-party libraries in this project, which are installed at build-time.

## Build options

The table below contains the flags that you can pass to the `cmake` command to control the build of the third-party libraries. Arguments to the `cmake` command must adhere to the format `-D<FLAG>=<VALUE>`. For examples on how to do this, see [Extended example](#extended-example).

Option                                   | Description                                               | Default
-----------------------------------------|-----------------------------------------------------------|----------------------------------
`CMAKE_INSTALL_PREFIX`                   | *Base installation directory for third-party libraries.*  | `<project_root>/install/`
`OOMPH_ENABLE_MPI`                       | *Enable MPI support?*                                     | `OFF`
`OOMPH_BUILD_OPENBLAS`                   | *Build OpenBLAS?*                                         | `ON`
`OOMPH_BUILD_SUPERLU`                    | *Build SuperLU?*                                          | `ON`
`OOMPH_BUILD_SUPERLU_DIST`               | *Build SuperLU DIST?*                                     | `ON` if MPI is enabled else `OFF`
`OOMPH_BUILD_CGAL`                       | *Build CGAL (with dep.  Boost)?*                          | `ON`
`OOMPH_BUILD_MUMPS`                      | *Build MUMPS?*                                            | `ON` if MPI is enabled else `OFF`
`OOMPH_BUILD_HYPRE`                      | *Build Hypre?*                                            | `ON`
`OOMPH_BUILD_TRILINOS`                   | *Build Trilinos?*                                         | `ON`
`OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS` | *Enable testing when building the third-party libraries?* | `OFF`
`OOMPH_USE_OPENBLAS_FROM`                | *The path to a preinstalled version of OpenBLAS.*         | `""`
`OOMPH_USE_BOOST_FROM`                   | *The path to a preinstalled version of Boost.*            | `""`

> **Note:** `CMAKE_INSTALL_PREFIX` is a built-in variable that is typically used to control the installation directory, hence why it does not have an `OOMPH_` prefix.

### Remark on `OOMPH_USE_<LIBRARY>_FROM` variables

The arguments to the `OOMPH_USE_<LIBRARY>_FROM` flags must be a folder containing a `lib/` and `include/` folder (and possibly even a `bin/` folder) for that library. As an example, consider the Boost library. Below we show how we locate the Boost installation of Homebrew

```bash
>>> brew --prefix boost
/opt/homebrew/opt/boost

>>> ls -1 /opt/homebrew/opt/boost
include/
lib/
share/
INSTALL_RECEIPT.json
README.md
sbom.spdx.json
```

To use this installation of Boost in the build of the third-party libraries, you would call `cmake` like so

```bash
>>> cmake -G Ninja -DOOMPH_USE_BOOST_FROM=/opt/homebrew/opt/boost -B build
```

or, more simply

```bash
>>> cmake -G Ninja -DOOMPH_USE_BOOST_FROM=$(brew --prefix boost) -B build
```

### Extended example

**Example 1:** To build all of the third-party libraries without any testing, run

```bash
>>> cmake -G Ninja -B build
>>> cmake --build build
```

**Example 2:** To enable MPI support and only build CGAL and OpenBLAS with all available testing, run

```bash
>>> cmake -G Ninja -DOOMPH_ENABLE_MPI=ON -DOOMPH_BUILD_CGAL=ON -DOOMPH_BUILD_OPENBLAS=ON -DOOMPH_BUILD_MUMPS=OFF -DOOMPH_BUILD_HYPRE=OFF -DOOMPH_BUILD_TRILINOS=OFF -DOOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS=ON -B build
>>> cmake --build build
```

Since, by default, we build all of the third-party libraries, you can actually reduce the above commands to

```bash
>>> cmake -G Ninja -DOOMPH_ENABLE_MPI=ON -DOOMPH_BUILD_MUMPS=OFF -DOOMPH_BUILD_HYPRE=OFF -DOOMPH_BUILD_TRILINOS=OFF -DOOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS=OFF -B build
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
