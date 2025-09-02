<div align="center">
  <a href="http://oomph-lib.maths.man.ac.uk">
    <img alt="reviewdog" src="./doc/figures/oomph_logo.png">
  </a>
</div>

<div align="center">
  <a href="./LICENCE">
    <img alt="License: LGPL v2.1" src="https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg">
  </a>
  <a href="../../actions/workflows/build-and-publish-docs.yaml">
    <img alt="Documentation" src="../../actions/workflows/build-and-publish-docs.yaml/badge.svg?branch=main">
  </a>
</div>

<br>

<div align="center">
    <table style="width:55%">
    <tr>
        <th>Test platform</th>
        <th><a href="../../tree/main"><code>main</code></a></th>
        <th><a href="../../tree/development"><code>development</code></a></th>
        <th><a href="../../tree/feature-shared-cmake"><code>feature-shared-cmake</code></a></th>
    </tr>
    <tr>
        <td>Ubuntu</td>
        <td>
            <a href="../../actions/workflows/self-tests-ubuntu.yaml">
                <img alt="Ubuntu" src="../../actions/workflows/self-tests-ubuntu.yaml/badge.svg?branch=main" style="vertical-align: middle">
            </a>
        </td>
        <td>
            <a href="../../actions/workflows/self-tests-ubuntu.yaml">
                <img alt="Ubuntu" src="../../actions/workflows/self-tests-ubuntu.yaml/badge.svg?branch=development" style="vertical-align: middle">
            </a>
        </td>
        <td>
          <div align="center">
              <a href="../../actions/workflows/self-tests-ubuntu.yaml">
                  <img alt="Ubuntu" src="../../actions/workflows/self-tests-ubuntu.yaml/badge.svg?branch=feature-shared-cmake" style="vertical-align: middle">
              </a>
          </div>
        </td>
    </tr>
    <tr>
        <td>macOS</td>
        <td>
            <a href="../../actions/workflows/self-tests-macos.yaml">
                <img alt="macOS" src="../../actions/workflows/self-tests-macos.yaml/badge.svg?branch=main" style="vertical-align: middle">
            </a>
        </td>
        <td>
            <a href="../../actions/workflows/self-tests-macos.yaml">
                <img alt="macOS" src="../../actions/workflows/self-tests-macos.yaml/badge.svg?branch=development" style="vertical-align: middle">
            </a>
        </td>
        <td>
          <div align="center">
              <a href="../../actions/workflows/self-tests-macos.yaml">
                  <img alt="macOS" src="../../actions/workflows/self-tests-macos.yaml/badge.svg?branch=feature-shared-cmake" style="vertical-align: middle">
              </a>
          </div>
        </td>
    </tr>
    </table>
</div>

<!-- Use <h2> tags to omit heading from table of contents -->
<h2>Table of contents</h2>

- [Overview](#overview)
- [Prerequisites](#prerequisites)
  - [CMake](#cmake)
  - [Ninja](#ninja)
  - [[macOS only: OpenBLAS]](#macos-only-openblas)
  - [Compilers](#compilers)
  - [Other tools (and how to install all the prerequisites)](#other-tools-and-how-to-install-all-the-prerequisites-in-one-go)
- [Building, installing and uninstalling `oomph-lib`](#building-installing-and-uninstalling-oomph-lib)
  - [Step 1: Installing the third-party libraries](#step-1-installing-the-third-party-libraries)
  - [Step 2: Installing oomph-lib](#step-2-installing-oomph-lib)
    - [Option 1: Default installation](#option-1-default-installation)
    - [Option 2: Specifying a custom installation location](#option-2-specifying-a-custom-installation-location)
    - [Option 3: Installation with root privileges](#option-3-installation-with-root-privileges)
    - [Further build options](#further-build-options)
- [Recommended alternative: Building with `oomph_build.py`](#recommended-alternative-building-with-oomph_buildpy)
  - [When to use the build script](#when-to-use-the-build-script)
  - [Using the script](#using-the-script)
  - [Command line options](#command-line-options)
  - [Running the script: example and output](#running-the-script-example-and-output)
- [Testing the installation](#testing-the-installation)
  - [Testing the entire installation](#testing-the-entire-installation)
  - [How to analyse failed tests](#how-to-analyse-failed-tests)
  - [Testing specific demo driver codes and how to analyse failed tests](#testing-specific-demo-driver-codes-and-how-to-analyse-failed-tests)
  - [Selective testing](#selective-testing)
    - [Filtering by `TEST_NAME`](#filtering-by-test_name)
    - [Filtering tests by keywords appearing in demo driver codes](#filtering-tests-by-keywords-appearing-in-demo-driver-codes)
    - [Disabling a test](#disabling-a-test)
- [How to run demo driver codes (and how to modify them when working on coding exercises, say)](#how-to-run-demo-driver-codes-and-how-to-modify-them-when-working-on-coding-exercises-say)
- [Adding your own driver code](#adding-your-own-driver-code)
- [Linking a stand-alone project to `oomph-lib`](#linking-a-stand-alone-project-to-oomph-lib)
  - [Creating the `CMakeLists.txt` file](#creating-the-cmakeliststxt-file)
  - [Customising driver codes](#customising-driver-codes)
    - [Customising driver codes via the `oomph_add_executable()` function](#customising-driver-codes-via-the-oomph_add_executable-function)
    - [Customising targets using native CMake commands; hashed target names](#customising-targets-using-native-cmake-commands-hashed-target-names)
    - [Using native CMake commands throughout](#using-native-cmake-commands-throughout)
- [Creating new demo driver directories](#creating-new-demo-driver-directories)
  - [Step 1: Rename and edit the driver code](#step-1-rename-and-edit-the-driver-code)
  - [Step 2: Update the CMakeLists.txt file](#step-2-update-the-cmakeliststxt-file)
  - [Step 3: Run and develop your own code](#step-3-run-and-develop-your-own-code)
  - [Step 4: Update the validation script and the validata](#step-4-update-the-validation-script-and-the-validata)
  - [Step 5: Register the new directory so it is included in a top-level build](#step-5-register-the-new-directory-so-it-is-included-in-a-top-level-build)
- [Creating a new library (or updating an existing one)](#creating-a-new-library-or-updating-an-existing-one)
  - [Updating an existing library](#updating-an-existing-library)
  - [Adding a new library](#adding-a-new-library)
- [The `doc` directory](#the-doc-directory)
- [Cheat sheet: Autotools vs. CMake](#cheat-sheet-autotools-vs-cmake)
  - [Build and install the entire library using a helper script](#build-and-install-the-entire-library-using-a-helper-script)
    - [Autotools](#autotools)
    - [CMake](#cmake-1)
  - [Alternative: Build and install the entire library step by step](#alternative-build-and-install-the-entire-library-step-by-step)
    - [Autotools](#autotools-1)
    - [CMake](#cmake-2)
  - [Run all self-tests](#run-all-self-tests)
    - [Autotools](#autotools-2)
    - [CMake](#cmake-3)
  - [Run all self-tests in a given demo driver directory](#run-all-self-tests-in-a-given-demo-driver-directory)
    - [Autotools](#autotools-3)
    - [CMake](#cmake-4)
  - [Debug demo driver code](#debug-demo-driver-code)
    - [Autotools](#autotools-4)
    - [CMake](#cmake-5)
  - [Edit a file in the library and rebuild](#edit-a-file-in-the-library-and-rebuild)
    - [Autotools](#autotools-5)
    - [CMake](#cmake-6)
- [Dos and Don'ts](#dos-and-donts)
- [FAQs](#faq)
  - [When configuring my driver code CMake can't find the package configuration file. Now what?](#when-configuring-my-driver-code-cmake-cant-find-the-package-configuration-file)
  - [Are there any complete worked examples of the build process?](#are-there-any-complete-worked-examples-of-the-build-process)
  - [What happened to the `user_src` directory?](#im-used-to-the-autotools-based-version-of-oomph-lib-what-happened-to-the-user_src-directory)
  - [What happened to the `bin` directory?](#im-used-to-the-autotools-based-version-of-oomph-lib-what-happened-to-the-bin-directory)
- [Additional information for developers](#additional-information-for-developers)
  - [Use symbolic links for header files](#use-symbolic-links-for-header-files)
  - [Paranoia and range checking are incompatible with Release mode](#paranoia-and-range-checking-are-deemed-to-be-incompatible-with-release-mode)
  - [How to add additional compiler macros to `oomph-lib` (and to stand-alone driver codes)](#how-to-add-additional-compiler-macros-to-oomph-lib-and-to-stand-alone-driver-codes)
  - [Creating robust `validata` for self tests](#creating-robust-validata-for-self-tests)
    - [A self-test fails even though the output files produced by the code are correct](#a-self-test-fails-even-though-the-output-files-produced-by-the-code-are-correct)
    - [Handling non-deterministic output](#handling-non-deterministic-output)
    - [Careful with driver codes that use `triangle` to generate meshes](#careful-with-driver-codes-that-use-triangle-to-generate-meshes)
  - [How to update third-party libraries to later versions](#how-to-update-third-party-libraries-to-later-versions)
- [Appendix](#appendix)
  - [CMake resources](#cmake-resources)
  - [Building CMake](#building-cmake)
    - [Ubuntu](#ubuntu)
    - [macOS](#macos)



## Overview

This document provides instructions how to install `oomph-lib` using our CMake-based build tools. We provide support for the following operating systems:

Operating system | Support provided?
-----------------|------------------
Ubuntu           | Yes
macOS            | Yes
Windows          | No

The [`oomph-lib` homepage](http://www.oomph-lib.org) provides much more detail on the library's capabilities, tutorials, and licencing information.

To learn more about contributing to `oomph-lib`, please see the document [`CONTRIBUTING.md`](CONTRIBUTING.md) which contains a detailed description of our GitHub-based workflow for adding and maintaining code.

## Prerequisites

We assume you have checked out `oomph-lib` from its [GitHub repository](https://github.com/oomph-lib/oomph-lib). The document [`CONTRIBUTING.md`](CONTRIBUTING.md) provides detailed instructions for how to interact with GitHub, but to get started you can simply clone the repository onto your computer

```bash
git clone https://github.com/oomph-lib/oomph-lib.git
```

`oomph-lib` uses CMake and Ninja to build and install the project. You will also need a sufficiently up-to-date compiler. The minimum requirements are as follows:

### CMake

You will need CMake 3.24+. You can install `cmake` via your native package manager, e.g.

```bash
# On Ubuntu:
sudo apt-get install cmake
```

or

```bash
# On macOS:
brew install cmake
```

If your package manager does not provide a sufficiently recent version of `cmake`, you will need to build it from source. You can find instructions on how to do this [below](#building-cmake) and on the [CMake website](https://cmake.org/install/).

> [!NOTE]
> We have heavily customised the CMake-based build process to facilitate the installation and usage of the library. Newcomers to CMake who wish to google details of the CMake syntax should note that any functions starting with the `oomph_` prefix are our own and are therefore not described in the native CMake documentation.

### Ninja

We strongly advocate the use of [Ninja](https://github.com/ninja-build/ninja)
for its automatic parallelisation of the build process. Ninja creates clear,
human-readable build files and allows for fast incremental builds.
Throughout this document we will assume that Ninja is installed.
Again, the installation is straightforward if you don't already have it on on your machine. For example, on Ubuntu

```bash
sudo apt install ninja-build
```

will do the trick.

### [macOS only: OpenBLAS]

> [!IMPORTANT]
> Building OpenBLAS as part of the third-party libraries build is not currently supported on macOS. Instead, you need to install it with a package manager, e.g.
>
> ```bash
> brew install openblas
> ```
>
> then set
>
> ```bash
> -DOOMPH_USE_OPENBLAS_FROM=$(brew --prefix openblas)
> ```
>
> during the project configuration step (see below). Adding OpenBLAS installation support on macOS is a work in progress.

### Compilers

`oomph-lib` is written entirely in C++ and we require a compiler that can handle the C++17 standard. We also use third-party libraries written in C and Fortran so your compiler suite must be able to handle these too. On Ubuntu you may have to install `gfortran`, using

```bash
sudo apt install gfortran
```
If you wish to use `oomph-lib`'s parallel capabilities you need a working MPI installation on your machine. We tend to use [OpenMPI](https://www.open-mpi.org/) which on Ubuntu can be installed using
```bash
sudo apt install openmpi-bin libopenmpi-dev
```

### Other tools (and how to install all the prerequisites in one go)
If you also want to build the documentation (which will give you a local copy of the `oomph-lib` webpages) you need a few other tools. On Ubuntu the following command should give you all you need:
```bash
sudo apt-get install git cmake ninja python3 doxygen gfortran g++ texlive texlive-latex-extra texlive-font-utils openmpi-bin libopenmpi-dev
```


## Building, installing and uninstalling `oomph-lib`

`oomph-lib` relies on the following third-party libraries:

| Library                          | Version                                                                                                                                   |
|----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `OpenBLAS` (**required**)        | [0.3.25](https://github.com/OpenMathLib/OpenBLAS/tree/v0.3.29)                                                                            |
| `Boost` (**highly recommended**) | [1.83.0](https://github.com/boostorg/boost/tree/boost-1.83.0)                                                                             |
| `CGAL` (**highly recommended**)  | [6.0.1](https://github.com/CGAL/cgal/tree/v6.0.1)                                                                                         |
| `GKlib`                          | [commit `8bd6bad750b2b0d90800c632cf18e8ee93ad72d7`](https://github.com/KarypisLab/GKlib/tree/8bd6bad750b2b0d90800c632cf18e8ee93ad72d7)    |
| `METIS`                          | [commit `e0f1b88b8efcb24ffa0ec55eabb78fbe61e58ae7`](https://github.com/KarypisLab/METIS/tree/e0f1b88b8efcb24ffa0ec55eabb78fbe61e58ae7)    |
| `ParMETIS`                       | [commit `8ee6a372ca703836f593e3c450ca903f04be14df`](https://github.com/KarypisLab/ParMETIS/tree/8ee6a372ca703836f593e3c450ca903f04be14df) |
| `SuperLU`                        | [v6.0.1](https://github.com/xiaoyeli/superlu/tree/v6.0.1)                                                                                 |
| `SuperLU_DIST`                   | [v9.1.0](https://github.com/xiaoyeli/superlu_dist/tree/v9.1.0)                                                                            |
| `MUMPS`                          | [5.6.2](https://github.com/puneetmatharu/mumps/tree/v5.6.2.5)                                                                             |
| `HYPRE`                          | [2.32.0](https://github.com/hypre-space/hypre/tree/v2.32.0)                                                                               |
| `Trilinos`                       | [16.0.0](https://github.com/trilinos/Trilinos/tree/trilinos-release-16-0-0)                                                               |

> [!IMPORTANT]
> If you are an Apple user, make sure you read our instructions for installing [OpenBLAS](#macos-only-openblas).

To facilitate the installation of these (as well as various optional third-party libraries), we provide the option to build them as part of our overall build process. We provide a detailed description of the two-stage build process below, but strongly encourage users to use our `oomph_build.py` script that performs all these actions in one operation. Details are described in the section [Building with `oomph_build.py`](#recommended-alternative-building-with-oomph_buildpy) and we suggest that new users jump straight there.

For everybody still reading, the two stages are:

1. [Installing the third-party libraries](#Installing-the-third-party-libraries)
1. [Installing oomph-lib](#Installing-oomph-lib)

### Step 1: Installing the third-party libraries

When cloning the code from GitHub, the third-party libraries (required and optional) are contained in the directory `external_distributions`. We recommend installing all third-party libraries which is as easy as

```bash
# Go to the third-party distributions directory
cd external_distributions

# Configure the third party library build
cmake -G Ninja -B build

# Build and install them (this takes a while). Note that no separate
# install step is required here. By default the libraries are
# installed locally in the external_distributions/install directory
cmake --build build
```

> [!TIP]
> Make sure you delete the `build` directory if it already exists from a previous installation. Having it there can confuse `cmake`.

At the end of the configure step you get a summary of the build options used, e.g.

```bash
********************************************************************************
OOMPH-LIB THIRD-PARTY LIBRARIES OPTIONS:
********************************************************************************
  ⦿  CMAKE_INSTALL_PREFIX: '/home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install'
  ⦿  CMAKE_BUILD_TYPE: 'Release'
  ⦿  OOMPH_ENABLE_MPI: 'OFF'
  ⦿  OOMPH_BUILD_OPENBLAS: 'ON'
  ⦿  OOMPH_BUILD_SUPERLU: 'ON'
  ⦿  OOMPH_BUILD_SUPERLU_DIST: 'OFF'
  ⦿  OOMPH_BUILD_CGAL: 'ON'
  ⦿  OOMPH_BUILD_MUMPS: 'ON'
  ⦿  OOMPH_BUILD_HYPRE: 'ON'
  ⦿  OOMPH_BUILD_TRILINOS: 'ON'
  ⦿  OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS: 'OFF'
  ⦿  OOMPH_USE_OPENBLAS_FROM: ''
  ⦿  OOMPH_USE_GKLIB_FROM: ''
  ⦿  OOMPH_USE_METIS_FROM: ''
  ⦿  OOMPH_USE_PARMETIS_FROM: ''
  ⦿  OOMPH_USE_BOOST_FROM: ''
********************************************************************************
```

You'll see that by default we build the libraries in Release mode (i.e. with full optimisation) and in serial (i.e. without any MPI support), and then install them in `external_distributions/install`. This is done deliberately to make sure that even users who do not have root privileges on their computers can install the libraries.

Given that these are third-party libraries it's unlikely that you'll want to build them in (slow) Debug mode but you could do so by specifying the flag `-DCMAKE_BUILD_TYPE=Debug` at the configure stage. If you have MPI installed on your machine and wish to use the parallel capabilities of some of the third-party libraries, specify `-DOOMPH_ENABLE_MPI=ON`. You can also specify a different location for the install directory using `-DCMAKE_INSTALL_PREFIX=/home/joe_cool/oomph-lib_third_party_libraries`, say. Finally, if your computer already has some (possibly optimised) version of the third-party libraries installed, you can use them instead of rebuilding those libraries yourself. So, if you have an optimised installation of OpenBLAS installed at `/home/joe_cool/super_fast_openblas`, say, you can use that by configuring with `-DOOMPH_USE_OPENBLAS_FROM=/home/joe_cool/super_fast_openblas`.

> [!TIP]
> It is always a good idea (and sometimes required) to use absolute paths when specifying such directories!

Here's an example:

```bash
# Configure in Debug mode and with MPI support. Use an existing OpenBLAS
# installation and install the others in a specific directory:
cmake -G Ninja -B build \
   -DCMAKE_BUILD_TYPE=Debug \
   -DOOMPH_ENABLE_MPI=ON \
   -DCMAKE_INSTALL_PREFIX=/home/joe_cool/oomph-lib_third_party_libraries \
   -DOOMPH_USE_OPENBLAS_FROM=/home/joe_cool/super_fast_openblas

# Now build/install it
cmake --build build
```

> [!IMPORTANT]
> When specifying the directory that contains an existing third-party library, make sure that it contains the `lib` and `include` directories of the relevant library, so for OpenBLAS you should see this:

```bash
ls -l /home/joe_cool/super_fast_openblas
total 12
drwxrwxr-x 2 joe_cool  joe_cool  4096 May 27 11:53 include
drwxrwxr-x 4 joe_cool  joe_cool  4096 May 27 11:53 lib
```

Following the successful installation of the third-party libraries, we provide a summary of how to use them in the subsequent `oomph-lib` build.

### Step 2: Installing oomph-lib

Now that the required third-party libraries are available, we proceed to the actual `oomph-lib` build process. Again we stress that inexperienced (or lazy) users may prefer to use the `oomph_build.py` script, described [below](#recommended-alternative-building-with-oomph_buildpy). For those who want to do the installation manually, we discuss three different options. The process is obviously very similar to the generic CMake build process already used for the third-party libraries. The key here is to specify the location of the third-party libraries -- either those just built or those already pre-installed elsewhere. For simplicity we initially assume that the third-party libraries were built as described above. The install directory (irrespective of its location) then contains a file `cmake_flags_for_oomph_lib.txt` which contains all the relevant information about the location of the third-party libraries. It can be passed directly to the configure stage of the `oomph-lib` build. In fact, when the build process for the third-party libraries finishes it lists a variety of equivalent CMake commands that can be used to configure the subsequent `oomph-lib` build.

#### Option 1: Default installation

To configure, build and install the library, `cd` into the root directory of the cloned `oomph-lib` project (checked out from GitHub) and run the following commands:

```bash
# Configure and generate the build system. -G specifies the build
# system generator; -B specifies the build directory (here "build").
# The final argument specifies the location of the third-party
# libraries (built in the  previous step; the command shown here
# assumes that the third party libraries were installed in the
# default location external_distributions/install. Update this if
# you installed them somewhere else).
cmake -G Ninja -B build $(cat external_distributions/install/cmake_flags_for_oomph_lib.txt)

# Build the oomph-lib libraries (i.e. compile the sources
# and turn them into libraries); specify the directory
# that was created at the configure stage (here "build").
# This takes a while...
cmake --build build

# Install, again specifying the build directory created above
# (here "build")
cmake --install build
```

Note that, by default, the last step installs the headers and the generated library files into a directory named `install/` in the `oomph-lib` root directory. This is done deliberately to make sure that even users who do not have root privileges on their computers can install the library.

To uninstall the library simply delete the `build/` and `install/` directories in the `oomph-lib` root directory.

#### Option 2: Specifying a custom installation location

By default, `oomph-lib` will be installed to the `install/` subdirectory of the root `oomph-lib` directory. To specify a different installation location, pass `-DCMAKE_INSTALL_PREFIX=<install-location>` to `cmake` during the configure step. Note that `<install-location>` must be an **absolute** path. Here's an example where we specify different build and install directories:

```bash
# Configure and build (write build information into the
# directory "joe_build") and use the install prefix to
# declare where the library should ultimately be installed.
cmake -G Ninja -B joe_build -DCMAKE_INSTALL_PREFIX=/home/joe_cool/oomph_lib_install_dir

# Build the library; specify the directory created during the
# configure stage.
cmake --build joe_build

# ...and install it; again specify the build directory created
# above (here "joe_build"). The library will now be installed
# in /home/joe_cool/oomph_lib_install_dir
cmake --install joe_build
```

To make sure that the library can be found when building driver codes (see below) it is easiest to set the `oomphlib_ROOT` environment variable:

```bash
# Export the custom installation directory
export oomphlib_ROOT=/home/joe_cool/oomph_lib_install_dir
```

This is, of course, only a temporary assigment for the current shell; add the command to your `.bashrc` file (or `.zshrc` file, if using Zsh) to reassign it automatically for every session. Alternatively, you can specify the location of the install directory with the `-Doomphlib_ROOT` flag when configuring the driver codes; see [below](#testing-the-entire-installation).

To uninstall the library, delete the `joe_build` directory in the `oomph-lib` root directory and the directory `/home/joe_cool/oomph_lib_install_dir`. Also remove the
addition of the installation directory to the `PATH` from your `.bashrc` or `.zshrc` file (if you added it there).

#### Option 3: Installation with root privileges

If you have root privileges on your computer, you can install `oomph-lib` in the standard location `/usr/local/`:

```bash
# Configure; using a local build directory but prepare
# for installation in /usr/local
cmake -G Ninja -B build -DOOMPH_ALLOW_INSTALL_AS_SUPERUSER=ON

# Build (i.e. compile the sources and create the libraries)
cmake --build build

# Install -- here you need to use your superuser powers (by
# using sudo) because this will install the library
# into /usr/local
sudo cmake --install build
```

Note that `/usr/local/` is a standard location for libraries, it is therefore not necessary to add the location of the install directory to the `PATH` variable.

Deleting things from `/usr/local` is scary, therefore it's best to uninstall `oomph-lib` using

```bash
cd build/
sudo ninja oomph_uninstall
```

where you again have to use your root privileges (`sudo`).

#### Further build options

You can customise your build by passing flags of the form `-D<FLAG>` to `cmake` during the configuration stage. For reference, the table below contains a list of options that the user can control. The build and installation steps will remain the same.

| FLAG                                          | Description                                                                                         | Default                          |
|-----------------------------------------------|-----------------------------------------------------------------------------------------------------|----------------------------------|
| `CMAKE_BUILD_TYPE`                            | The build type (`Debug`, `Release`, `RelWithDebInfo` or `MinSizeRel`)                               | `Release`                        |
| `OOMPH_PRINT_SETTINGS_AFTER_INSTALL`          | Display the oomph-lib settings at the end of the install step                                       | `ON`                             |
| `OOMPH_ALLOW_INSTALL_AS_SUPERUSER`            | Allow the user to install to the default system install path (if `CMAKE_INSTALL_PREFIX` is not set) | `OFF`                            |
| `OOMPH_INSTALL_HEADERS_AS_SYMLINKS`           | Install symlinks to the oomph-lib headers instead of copying them                                   | `OFF`                            |
| `OOMPH_DONT_SILENCE_USELESS_WARNINGS`         | Display (harmless) warnings from external_src/ and src/ that are silenced                           | `OFF`                            |
| `OOMPH_ENABLE_MPI`                            | Enable the use of MPI for parallel processing                                                       | `OFF`                            |
| `OOMPH_MPI_NUM_PROC`                          | Number of processes to use with MPI-enabled tests                                                   | 2                                |
| `OOMPH_ENABLE_MPI_OVERSUBSCRIPTION`           | Enable oversubscription with MPI                                                                    | `OFF`                            |
| `OOMPH_ENABLE_PARANOID`                       | Enable the `PARANOID` flag in Debug                                                                 | `ON` if `Debug` build else `OFF` |
| `OOMPH_ENABLE_RANGE_CHECKING`                 | Enable `RANGE_CHECKING` flag in Debug                                                               | `ON` if `Debug` build else `OFF` |
| `OOMPH_ENABLE_SANITIZER_SUPPORT`              | Enable sanitizer support in Debug builds                                                            | `OFF`                            |
| `OOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER` | Use MumpsSolver as the default linear solver (when available)                                       | `OFF`                            |
| `OOMPH_SUPPRESS_TRIANGLE_LIB`                 | Suppress build of oomph-lib's copy of the triangle library                                          | `OFF`                            |
| `OOMPH_SUPPRESS_TETGEN_LIB`                   | Suppress build of oomph-lib's copy of the tetgen library                                            | `OFF`                            |

For instance, to build `oomph-lib` with MPI support, use

```bash
cmake -G Ninja -B build -DOOMPH_ENABLE_MPI=ON
```

Note that by default we build the library in Release mode, so the above command will build a version of the library that is compiled with full optimisation and no Debug support.

Other options are used mainly when developing new code. In this case, we recommend building with range checking (which makes the code run very slowly but detects array overruns) and "paranoia" which activates many internal self-tests throughout the library. Both will issue explicit error messages if problems are detected. However, these only tend to be useful if you can then reconstruct how the offending line of code was reached. For this purpose, we recommend building the library in Debug mode to ensure that a debugger like `ddd` can be used to backtrace from the error.

We therefore recommend doing code development with the following configure options:

```bash
cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=Debug -DOOMPH_ENABLE_PARANOID=ON -DOOMPH_ENABLE_RANGE_CHECKING=ON
```

The last two flags define the C++ macros `PARANOID` and `RANGE_CHECKING` which should be used as follows:

```c++
[...]

double some_function(const double& x)
 {
#ifdef PARANOID
   if (x<0.0)
   {
    throw OomphLibError("Muppet! I can't take the square root of a negative number!",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
#endif
   return sqrt(x);
 }

 [...]
```

The `RANGE_CHECKING` macro is mainly used within `oomph-lib`'s own containers but can/should also be used for any other newly-created objects where an array-like index has to be in a certain range.

Finally, developers should install the library headers as symbolic links to the actual header files. If you don't do this you'll go insane while fixing bugs in copies of header files that immediately get over-written when the (supposed to be fixed) library is reinstalled. So install the library as follows if you intend to work on it.

```bash
cmake -G Ninja -B build -DOOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON
```

Obviously don't do this if you want to install the library in a permanent location and then delete the sources!

## Recommended alternative: Building with `oomph_build.py`

This section explains how to use the `oomph_build.py` script -- a Python helper script located in the root of the repository – to configure, build, and install `oomph-lib` in an automated way.

### When to use the build script

The `oomph_build.py` script is provided as a convenient one-stop solution for building and installing `oomph-lib`. It essentially wraps the manual CMake steps (configuration, build, and installation) into a single command, with sensible defaults and additional helper options. New users are encouraged to use this script for a hassle-free installation, as it reduces the chance of missing a step or incorrectly specifying an option. Advanced users or those integrating `oomph-lib` into other build systems may still prefer the direct CMake approach (as in Options 1–3 above) for finer control, but the script covers the common use cases.

> [!IMPORTANT]
> In addition to the standard requirements (CMake and Ninja, as described above), you will need a working Python 3 installation to run `oomph_build.py`. Ensure that `python3` is available in your `PATH`. If you plan to use the script's documentation-building feature (described below), make sure `doxygen` is installed as well.

### Using the script

To use the build script, `cd` to the top-level `oomph-lib` directory (the same directory that contains `oomph_build.py`) and run it with Python. For example:

```python
python3 oomph_build.py
```

By default, this will build the third-party libraries, configure the `oomph-lib` project, compile its libraries, and install them to the local `install/` directory (within the `oomph-lib` root directory). The script will create a build directory (named `build/` by default) if one does not exist. Note that CMake can get confused by existing `build` or `install` directories from previous installations. Therefore we generally stop the process with an error if these directories already exist. (It is possible to reuse them by specifying suitable command line options; see below.) The script thus combines the equivalent of the `cmake -B build`, `cmake --build build`, and `cmake --install build` commands into one automated sequence.

After running `oomph_build.py`, if it completes without errors, `oomph-lib` should be built and installed in the specified location. You can then proceed to test the installation as described in the next section.

> [!NOTE]
> On Unix-like systems, you can make the script executable by running `chmod +x oomph_build.py` and then invoking it with `./oomph_build.py`. Using `python3 ...` is equivalent to `./...` (after making the script executable).

### Command line options

The `oomph_build.py` script supports several command-line options to customise its behavior. You can see a summary by running

```bash
python3 oomph_build.py --help
```

Below is a list of the options and their purpose:

- **`--build-doc`**: By default, the script does not build the documentation (to save time and space). If you want to generate the full HTML (and PDF) documentation for oomph-lib, include the `--build-doc` flag. This will build the documentation (using Doxygen and LaTeX) as part of the build process. *Note:* Building documentation can be time-consuming and requires Doxygen (and a LaTeX distribution for PDFs) to be installed. If these tools are missing, the doc build will fail – in that case, either install the necessary tools or omit this option.

- **`--oomph-CMAKE_INSTALL_PREFIX`**: Use this option when you intend to provide a custom installation location.
- **`--oomph-OOMPH_ALLOW_INSTALL_AS_SUPERUSER`**: Use this option when you intend to install `oomph-lib` system-wide (by installing it in `/usr/local/`) using root privileges. This flag tells the script to configure the installation prefix to the system's default location (`/usr/local`) and to set any required internal CMake switches such as`OOMPH_ALLOW_INSTALL_AS_SUPERUSER=ON`. If you also specify `--oomph-CMAKE_INSTALL_PREFIX`, this flag will have no effect.)

> [!IMPORTANT]
> When using this option, you will need to run the installation step with administrative privileges. The script will attempt to perform the installation step with sudo if possible, or it will remind you to re-run the script as root for the install phase. It is generally recommended to run oomph_build.py `--oomph-CMAKE_INSTALL_PREFIX` as a normal user for the build, and let it prompt for a password or instruct you for the install, rather than running the entire build as root. (Building as a non-root user helps avoid permission issues in the build directory.)

- **`--wipe-tpl`**, **`--wipe-root`**, **`--wipe-doc`**: These options tell the script to remove the specified build/installation directories that would be written to when building the third-party libraries, root project, and documentation, respectively. Use the `--wipe-*` flags if you want a completely clean rebuild. For example

  ```bash
  # Wipe the default build/installation directories used when building
  # the third-party libraries, root project, and documentation
  python3 oomph_build.py --wipe-tpl --wipe-root --wipe-doc
  ```

  will delete the current `build/` directory and the `install/` directory (if they exist) before configuring a fresh build. `--wipe-doc` will remove the `doc/build/` directory. These options are useful if you suspect a previous build is causing issues or if you want to reclaim space and rebuild from scratch.

> [!CAUTION]
> Wiping will permanently delete those directories
(and all compiled files or installed files therein), so use these flags with care.

- **`--reuse-tpl`**, **`--reuse-root`**, **`--reuse-doc`**: These flags are the counterpart to the `--wipe-*` options, explicitly instructing the script to reuse existing directories. Using `--reuse-*` flag means the script will allow the existing build/installation directories and will update or rebuild files as needed. These options can save time on iterative builds, but be mindful that reusing directories might carry over old artifacts.

> [!TIP]
> If you ever run into strange problems during the configuration/build process, try wiping and doing a clean new build.

In summary, by default the script uses a safe approach (an error is thrown if the build/installation directories already exist), requiring you to re-run `oomph_build.py` with `--wipe-*` or `--reuse-*` flags. If a completely fresh build environment is desired, use the `--wipe-*` options.

### Running the script: example and output

When you run `oomph_build.py`, it will print a few high-level messages indicating its progress through the build stages. You can obtain (much) more detailed output by running `oomph_build.py` with the `--verbose` flag.

If the configuration is successful, you should see the message

```bash
-- Build files have been written to: [..]/build
```

If there are errors (e.g., a required tool is not found or a configuration check failed), CMake will report them; the script will stop and propagate that error.

Next, the script invokes the build stage (equivalent to running `cmake --build build`). If you specified the `--verbose` flag you will see the compiler output as the `oomph-lib` libraries compile. This stage may take some time, especially on a first build. If it succeeds, all necessary library files will be produced in the build directory.

Finally, the script runs the installation stage. You'll see messages about installing libraries and header files to the target directory.

## Testing the installation

### Testing the entire installation

`oomph-lib` comes with an extensive list of well-documented example driver codes situated in the `demo_drivers/` directory. The driver codes in these directories are also used to validate the library.

Once the library is installed you can (and should!) enter the `demo_drivers/` directory and run the self-tests as follows:

```bash
# Go to the demo_drivers directory
cd demo_drivers

# Configure. Note that the directory specified here
# with the -B flag is a new directory that's
# created in the demo_drivers directory. It's generally
# a good idea to call it "build" but here we use
# a different name just to demonstrate that it has
# nothing to do with the build directory we used to
# build the library itself
cmake -G Ninja -B build_for_demo_drivers

# Go into the build directory specified at the configure
# stage
cd build_for_demo_drivers

# Now run the self-tests; here in parallel using
# 4 cores (processes)
ctest -j4
```

While `ctest` is running, we provide an overview of the progress and document the run-times of individual tests. Finally, a summary of the passed/failed tests is provided.

Note that running the tests creates a large amount of data, so once all the tests have completed (successfully!) you can delete the build directory:

```bash
# Exit the build_for_demo_drivers directory
cd ..

# Delete the build directory
rm -rf build_for_demo_drivers
```

Again, remember that if you installed `oomph-lib` to a non-standard location (i.e. neither under `install/` in the `oomph-lib` root directory, nor in `/usr/local/`) you must have added the installation directory to your `oomphlib_ROOT` environment variable (see [Option 2: Specifying a custom installation location](#option-2-specifying-a-custom-installation-location)). Alternatively, you can specify it during the configure stage, by modifying the above procedure to:

```bash
# Go to demo_drivers directory
cd demo_drivers

# Configure. Specify the (non-standard) location of the
# install directory
cmake -G Ninja -B build -Doomphlib_ROOT=/home/joe_cool/oomph_lib_install_dir

# Go into the build directory specified at the configure
# stage
cd build

# Now run the self-tests; here in parallel using
# 4 cores (processes)
ctest -j4
```

### How to analyse failed tests

Following the completion of the self-tests we provide a summary of run/passed/failed tests. Further information can be found in the log file
`demo_drivers/build/Testing/Temporary/LastTest.log` which provides a verbose log of all the tests built/run. Search for the string 'FAILED` within that file to find the failed tests.

To re-run only the failed tests with verbose output use the `--rerun-failed` and `--output-on-failure` flags:

```bash
cd demo_drivers/build
ctest --rerun-failed --output-on-failure
```

To investigate/debug/fix the failed test, `cd` into the relevant `demo_drivers` directory and follow the instructions in the next section.

### Testing specific demo driver codes and how to analyse failed tests

It is also possible to test a single driver code (or rather all the demo drivers below a specific subdirectory of `demo_drivers`):

```bash
# Go to a specific sub-directory in demo_drivers
cd demo_drivers/poisson/one_d_poisson/

# Configure (-B specifies the directory in which tests will be
# built/run). Here we assume that oomph-lib was installed in its
# default location (the install directory in the oomph-lib
# root directory; if you specified a different directory
# for the installation you have to specify it here, using the
# -Doomphlib_ROOT flag.
cmake -G Ninja -B build

# Go into the new build directory
cd build

# Run test(s)
ctest
```

The test procedure performed when running `ctest` in the `build` directory involves the following steps:

- The executables are created by Ninja.
- The script `validate.sh` (accessed via a symbolic link) is executed. The script typically runs the executables, sometimes repeatedly, with different command line options. Selected output from the test runs is then usually concatenated into a file.
- The file with the representative output is compared (using the script `scripts/fpdiff.py`, allowing for some specifiable floating point differences) against the compressed reference data in the `validata` directory (also accessed via a symbolic link). The test is deemed to have passed if the differences between computed and reference data are sufficiently small. A log of the run is contained in the file `validation.log`.

The tests are run in the `Validation` directory where the full output of the runs can be inspected. If a test fails, first establish which code (and/or which command line options) caused the failure by inspecting the `validate.sh` script. Inspect the associated output in the `Validation` directory and then start your debugging by re-running the failing code with the appropriate command line arguments (if any). Good luck! If you find an error please report it by submitting a GitHub Issue in the `oomph-lunch` [repository](https://github.com/oomph-lib/oomph-lib). If you find a bug fix, even better: submit a pull request and help us make the library even better!

### Selective testing

Developers often wish to test if their latest changes to the library have broken specific driver codes. There are various mechanisms to identify potentially affected codes and to run the self-tests only on these:

#### Filtering by `TEST_NAME`

If you inspect the `CMakeLists.txt` files in the `demo_drivers` directories you'll notice that each executable is associated with a `TEST_NAME`, specified in the `oomph_add_test(...)` function. The test name usually mirrors the name of the directory containing the driver code (relative to the `demo_drivers` directory).  For example, in `demo_drivers/poisson/one_d_poisson/CMakeLists.txt` you will see the following:

```cmake
oomph_add_test(
  TEST_NAME poisson.one_d_poisson
  DEPENDS_ON one_d_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)
```

The `TEST_NAME` can be specified as a regular expression to `ctest`, using the `-R`/`--tests-regex` flag. Only tests for which the `TEST_NAME` key contains the specified regular expression will be run. CMake has its own conventions for forming regular expressions and we suggest you consult the relevant CMake documentation. However, for simple cases this is extremely straightforward; for example

```bash
cd demo_drivers/build
ctest -R 'poisson.one_d_poisson'
```

will run all the self-tests whose `TEST_NAME` includes the string `poisson.one_d_poisson`. This approach does, of course, assume that you know the relevant `TEST_NAME`s. You could explore these for all the existing demo driver codes by using the command

```bash
find . -name 'CMakeLists.txt' -exec grep -H TEST_NAME {} \;
```

#### Filtering tests by keywords appearing in demo driver codes

Searching for relevant demo driver codes by their associated `TEST_NAME`, as explained in the previous section, is useful but the test name rarely includes all the features of the relevant code. A common scenario for developers is that, following changes to a specific class, they wish to test if all demo driver codes that explicitly use this class still work. To facilitate this we provide a helper script

```bash
scripts/filter_cmake_tests.py
```

that generates the `ctest` command required to run the self-tests in all `demo_drivers` directories containing driver codes that contain a specified string. So, to generate the `ctest` command that tests all demo drivers that contain the string 'QAdvectionDiffusionElement', say, run the following command

```bash
cd demo_drivers
../scripts/filter_cmake_tests.py --root . --keywords 'QAdvectionDiffusionElement'
```

Note the `--root` argument which specifies the location of the `demo_drivers` directory relative to the current directory.

#### Disabling a test

To (temporarily) disable a test in a given demo driver directory, edit the `CMakeLists.txt` file and
set the `DISABLED` property to `TRUE`, identifying the test using its `TEST_NAME`, as shown here:

```cmake
# Test definition
oomph_add_test(TEST_NAME poisson.one_d_poisson ...)

# Disable
set_tests_properties(poisson.one_d_poisson PROPERTIES DISABLED YES)
```

Disabled tests are listed explicitly in the final summary provided when `ctest` is run on the entire `demo_drivers` directory.

## How to run demo driver codes (and how to modify them when working on coding exercises, say)

The best way to get started with `oomph-lib` is to explore some of the demo driver codes. Many of these codes are explained in great detail in the associated [tutorials](https://oomph-lib.github.io/oomph-lib/doc/example_code_list/html/index.html) which typically end with a few exercises that encourage you to modify the code.

We will illustrate how to do this for the code in `demo_drivers/poisson/one_d_poisson_generic_only` which is the subject of [a specific tutorial](https://oomph-lib.github.io/oomph-lib/doc/quick_guide/html/index.html). When checked out of the GitHub repository, the directory contains the following files:

```bash
>>> ls -1
CMakeLists.txt
one_d_poisson_generic_only.cc
validata
validate.sh
```

The script `validate.sh` and the directory `validata` are only used during the self-tests (for obvious purposes). There is a driver code called `one_d_poisson_generic_only.cc`. Instructions for how to turn this into an executable are encoded in the local `CMakeLists.txt` file which contains

```bash
[...]

oomph_add_executable(
  NAME one_d_poisson_generic_only
  SOURCES one_d_poisson_generic_only.cc
  LIBRARIES oomph::meshes oomph::generic)

[...]
```

where we have omitted some boilerplate CMake code that is not relevant in the present context. The `oomph_add_executable(...)` statement defines:

- the name of the executable, here `one_d_poisson_generic_only`
- the local source code it depends on, here `one_d_poisson_generic_only.cc`
- the `oomph-lib` libraries used, here the `generic` library and `meshes` library (prefixed with the `oomph::` namespace identifier).

The only information required for the current task is the name of the executable.

To build and run the executable, we start with the usual configure step which specifies a (local) build directory:

```bash
# Configure a build directory
cmake -G Ninja -B build
```

We then `cd` into that directory and use `ninja` to create the executable, using `ninja <name_of_executable>`,  so here

```bash
# Enter the newly-created build directory
cd build

# Build the executable (the name is specified by the statement
# in the CMakeLists.txt file). NOTE: If you don't specify an
# executable, ninja will build all the executables listed in the
# CMakeLists.txt file.
ninja one_d_poisson_generic_only

# Now create any directories that may be needed
# and run the code
mkdir RESLT
./one_d_poisson_generic_only
```

Given that this code also doubles as a self-test, we can run the test in the same directory, simply using

```bash
# Run self-test
ctest
```

You can now edit the source, while working your way through the exercises in the tutorial. Ninja will only recompile the code (when you run `ninja`) if any of the sources have changed.

```bash
# Edit the source code, e.g increasing the number of
# elements from 10 to 20
emacs ../one_d_poisson_generic_only.cc

# Rebuild and run
ninja one_d_poisson_generic_only
./one_d_poisson_generic_only
```

Once the driver code is modified, the self-tests will obviously not pass any more (try running `ctest` again!). It may therefore be a better idea to work with a copy of the code, so revert back to the original version and make a copy

```bash
# Back to demo_drivers/poisson/one_d_poisson_generic_only
cd ..

# Revert code to its pristine status
git checkout one_d_poisson_generic_only.cc

# Make a copy
cp one_d_poisson_generic_only.cc my_code.cc
```

The self-tests should now pass again, because they operate with the original version of the code. The newly-created driver code now needs to be identified in the `CMakeLists.txt` file by adding

```bash
[...]

oomph_add_executable(
  NAME my_code
  SOURCES my_code.cc
  LIBRARIES oomph::meshes oomph::generic)

[...]
```

It is not necessary to reconfigure the build directory, so we can build and execute the new driver code straightaway (To be more precise: Ninja realises that the `CMakeLists.txt` has changed and therefore automatically triggers a reconfiguration.)

```bash
# Jump back into the existing build directory
cd build

# Build the new executable
ninja my_code

# Run it
./my_code
```

Now edit `my_code.cc` as you work your way through the exercises and rebuild the executable every time you change something using `ninja my_code`.

*Note for emacs users:*

If you're used to compiling from within emacs you've probably added something like this to your `.emacs` file:

```bash
(global-set-key `[F4]` 'compile)
(define-key global-map `[F12]` 'next-error)
```

With these keybindings, the `[F4]` function key will do a `make -k` in the current directory (you can edit the command line that appears at the bottom of the emacs window to add the name of a specific executable). If any errors occur during the compilation, `[F12]` will go from error to error in the source file(s). Very useful!

The fact that CMake separates the source and build directories means this won't quite work any more. However, assuming you adopt the common convention of calling your build directory `build/`, adding this

```bash
(setq compile-command "cd build; ninja")
```

to your `.emacs` file will produce equivalent behaviour. You can now edit the source code in its directory; `[F4]` will then compile in the build directory (again, you can specify the name of a specific executable by editing the command line that appears at the bottom of the emacs window); `[F12]` will then go through errors in the source file(s).

## Adding your own driver code

Developing your own code in an existing demo driver directory is a quick-and-dirty way to get started, especially since you are most likely to begin by modifying an existing driver code anyway. However, long-term this is not a sensible solution. One slightly more attractive alternative is to create a new directory, just for your code, in the `user_drivers` directory. This approach has the advantage of not interfering with existing `oomph-lib` driver codes and the associated test machinery. We already provide a sample directory `user_drivers/joe_cool` that shows you what to do. Simply add your own directory there and in the first instance copy across the files from the existing `joe_cool` directory. So assuming you're called Jack Cool you could do:

```bash
cd user_drivers
mkdir jack_cool
cp joe_cool/one_d_poisson.cc joe_cool/CMakeLists.txt jack_cool
```

You can then configure/build the code as usual:

```bash
cd jack_cool
cmake -G Ninja -B build
cd build
ninja
./one_d_poisson
```

It is not necessary to add the new sub-directory to `user_drivers/CMakeLists.txt`, though if you do, the modified version must not be checked into the library's GitHub repository! Now add additional driver codes to your directory `user_drivers/jack_cool` and add them to `user_drivers/jack_cool/CMakeLists.txt` file, following the embedded instructions.

## Linking a stand-alone project to `oomph-lib`

Adding your code into a new subdirectory in `user_drivers` directory has the advantage of your code being embedded into the `oomph-lib` build machinery. This is particularly useful for developers who tend to modify their driver codes and the library at the same time.

If you are not a developer and you simply wish to use `oomph-lib` as a library, you may prefer to have your driver code in a completely separate location. Given that you are then outside the `oomph-lib` framework you'll have to provide all kinds of additional information, particularly the location of the `oomph-lib` libraries and associated third-party libraries. CMake helps with that.

To illustrate the advantage of this approach, assume that you have a stand-alone driver code `non_oomph_lib_one_d_poisson.cc` that does not require the `oomph-lib` library. Using `g++` as the compiler you may want to compile this code using

```bash
g++ -o non_oomph_lib_one_d_poisson non_oomph_lib_one_d_poisson.cc
```

If the code has to be linked against an external library, e.g. `blas` which is installed in `/usr/local/blas` you'd specify this by adding

```bash
g++ -o non_oomph_lib_one_d_poisson non_oomph_lib_one_d_poisson.cc -I/usr/local/blas/include -L/usr/local/blas/lib -lblas
```

`oomph-lib` contains a lot of libraries and third-party libraries, so the equivalent command for a driver code `one_d_poisson.cc` that uses `oomph-lib` becomes

```bash

g++ -O3 -DNDEBUG -DKOKKOS_DEPENDENCE CMakeFiles/one_d_poisson_1f4eea6.dir/one_d_poisson.cc.o -o one_d_poisson  -Wl,-rpath,/home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/openblas/lib  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/libpoisson.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/libgeneric.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/liboomph_zlib.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/liboomph_hsl.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/liboomph_crbond_bessel.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/liboomph_triangle.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/install/lib/oomphlib/liboomph_tetgen.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/superlu/lib/libsuperlu.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/metis/lib/libmetis.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/gklib/lib/libGKlib.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/hypre/lib/libHYPRE.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/openblas/lib/libopenblas.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libModeLaplace.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libanasaziepetra.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libanasazi.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libml.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libifpack.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libamesos.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libaztecoo.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libtrilinosss.a  -lm  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libepetraext.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libtriutils.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libepetra.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchoskokkoscomm.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchoskokkoscompat.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchosremainder.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchosnumerics.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/openblas/lib/libopenblas.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/openblas/lib/libopenblas.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchoscomm.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchosparameterlist.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchosparser.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libteuchoscore.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libkokkossimd.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libkokkosalgorithms.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libkokkoscontainers.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/trilinos/lib/libkokkoscore.a  /usr/lib/x86_64-linux-gnu/libdl.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/mumps/lib/libsmumps.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/mumps/lib/libdmumps.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/mumps/lib/libmumps_common.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/mumps/lib/libpord.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/mumps/lib/libmpiseq.a  /home/joe_cool/oomph-lib_cmake/oomph-lib/external_distributions/install/openblas/lib/libopenblas.so  -lgfortran  -lquadmath
```

This is not something you want to enter (or maintain) by hand, so the trick is to use CMake to automatically obtain the relevant information from the `oomph-lib` installation.

### Creating the `CMakeLists.txt` file

To illustrate how this is done, assume you have a stand-alone directory that contains the driver codes and any associated header or include files needed to build the executable. Here's an example of such a project (checked out directly from its own GitHub repository)

```bash
>>> ls -1
CMakeLists.txt
glued_mesh_stuff.h
mesh_gluing2.cc
mesh_gluing.cc
```

It contains two driver codes, `mesh_gluing.cc` and `mesh_gluing2.cc`, both of which use the same include file, `glued_mesh_stuff.h`.

The `CMakeLists.txt` file contains the following:

```bash
#------------------------------------------------------------------------------

# Specify minimum cmake version; die if you can't find it
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

# Name of the project, followed by languages used;
# oomph-lib is written entirely in C++ but includes
# some third-party Fortran and C code
project(mesh_gluing CXX C Fortran)

# The find_package statement declares that:
# -- the project needs oomphlib (obviously!)
# -- we get the configuration (e.g. definition of
#    oomph_add_executable(...), used below) from the
#    oomph-lib project
# -- oomph-lib is required, so if you can't find
#    it die!
#
# Furthermore, we assume here that we've installed oomph-lib
# either in a standard location (/usr/local/) or that the
# install directory is available via the PATH
# environment variable oomphlib_ROOT, or that it will be
# be provided via the command line using the
# -Doomphlib_ROOT flag when configuring the project with
# cmake.
find_package(oomphlib CONFIG REQUIRED)

# Now add the first executable: Specify the name of the
# executable, the sources and the oomph-lib libraries
# required. (This is a solid mechanics problem so we need
# oomph-lib's generic library, the meshes library, the
# constitutive equations library, and the solid mechanics library)
oomph_add_executable(
  NAME mesh_gluing
  SOURCES mesh_gluing.cc glued_mesh_stuff.h
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# Now add the second executable: Specify the name of the
# executable, the sources and the oomph-lib libraries
# required. (This is a solid mechanics problem so we need
# oomph-lib's generic library, the meshes library, the
# constitutive equations library, and the solid mechanics library)
oomph_add_executable(
  NAME mesh_gluing2
  SOURCES mesh_gluing2.cc glued_mesh_stuff.h
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ------------------------------------------------------------------------------
```

> [!NOTE]
> In the above example we have included the header file `glued_mesh_stuff.h` into the `SOURCES` 
> variable. This is not strictly necessary in the sense that the code will compile 
> even if the included header files are not listed. However, they are needed for CMake to detect
> the executable's dependency on these files. It is therefore good practice to include them. 
> (This behaviour is different from the Autotools which scan code for included files to detect such dependencies automatically.)

Note that this directory is completely unconnected to the `oomph-lib` directory. To be able to find the `oomph-lib` installation directory, if it is not in one of the system-wide standard locations such as `/usr/local`, it must be declared somehow.

This is most easily done by exporting the `oomphlib_ROOT` environment variable, as discussed [above](#option-2-specifying-a-custom-installation-location). Note that this is necessary even if `oomph-lib` is installed to its default installation directory, `install/`, in the `oomph-lib` root directory -- remember that the current directory is completely unconnected to the `oomph-lib` installation. Alternatively, you can specify the directory where `oomph-lib` is installed when configuring the current project.

So, assuming `oomph-lib` was installed to `/home/joe_user/oomph_lib_install` the driver code can be built using

```bash
# Go to stand-alone driver directory
cd ~/mesh_gluing

# Option 1: Specify the oomph-lib installation directory during
# the configure step. Best because it's fully transparant
cmake -G Ninja -B -Doomphlib_ROOT=/home/joe_user/oomph_lib_install

# Option 2: Define/export the oomphlib_ROOT variable and configure.
# Not quite so good. It saves you having to type the install directory
# every time, but you may forget what you've set the variable
# to and thus acccidentally use an old installation. 
export oomphlib_ROOT=/home/joe_user/oomph_lib_install
cmake -G Ninja -B build

# Option 3: Open the CMakeLists.txt file and add the path to the
# oomph-lib installation to the find_package(...) call by replacing
#   find_package(oomphlib CONFIG REQUIRED)
# with
#   find_package(oomphlib CONFIG REQUIRED PATHS "/home/joe_cool/oomph_lib_install")
# then configure the project. This is considered bad practice since
# the build process will only work correctly on the present machine (or some
# other machine where oomph-lib happens to have been installed in this 
# specific directory)
cmake -G Ninja -B build

# Build; ninja without arguments builds both driver codes
cd build/
ninja

# Run 'em!
mkdir RESLT
./mesh_gluing
```

Note the three options for specifying the install directory (but stick to Option 1!).

Given that the install directory is quite deep, users sometimes worry if they 
have specified the right level of that directory tree. If you've only just installed `oomph-lib` 
it's easy to reconstruct:
- If you've built `oomph-lib` by running `oomph_build.py` (or the corresponding `cmake` 
  commands) in `/home/joe_cool/oomph-lib`, say, without explicitly specifying an install 
  directory, the relevant directory is `/home/joe_cool/oomph-lib/install`.
- If you specified a different install directory, either by running
  ```bash      
  ./oomph_build.py [...] --oomph-CMAKE_INSTALL_PREFIX=/home/joe_cool/local/oomph-lib
  ```
  or the raw CMake equivalent

  ```bash  
  cmake -G Ninja -B build [...] -DCMAKE_INSTALL_PREFIX=/home/joe_cool/local/oomph-lib
  ``` 
  then this is the one you use.

If you're unsure (because the installation was too long ago, say) you should check 
that the directory you're about to specify as the install directory contains the 
file `oomphlibConfig.cmake`. This is the file that contains the key information
that allows CMake to work with the installed library. So doing this
```bash
# Check that oomphlibConfig.cmake is there
find /home/joe_cool/local/oomph-lib -name 'oomphlibConfig.cmake'
```
should find the file
```bash
/home/joe_cool/local/oomph-lib/lib/cmake/oomphlib/oomphlibConfig.cmake
```
if the file can't be found (or it appears at a different level in the directory 
tree) the configuration will fail.

> [!NOTE]
> We provide a separate GitHub repository
>
> <div align="center"><a   href="https://github.com/oomph-lib/stand_alone_oomph-lib_user_code">https://github.com/oomph-lib/stand_alone_oomph-lib_user_code</a></div>
> <br>
>
> which contains a different driver code and a heavily annotated `CMakeLists.txt` file to explain the process in much more detail.

### Customising driver codes

#### Customising driver codes via the `oomph_add_executable()` function

The `CMakeLists.txt` shows how to specify the sources and the name of the executable to be built via the arguments `SOURCES` and `NAME`, respectively. They provide the minimum information required by the compiler. You may, of course, wish to provide additional information to the compiler. For instance, if your code was a stand-alone code that does not not involve `oomph-lib` or any other third-party libraries, you could compile it using

```bash
g++ -o non_oomph_lib_one_d_poisson -Wall -Werror -DREFINEABLE non_oomph_lib_one_d_poisson.cc
```

The additional flags passed to the compiler are of two types:

- `-Wall -Werror` are compiler flags that specify that the compiler should show all warnings and treat warnings as errors, respectively.
- `-DREFINEABLE` passes the compiler macro `REFINEABLE` to the compiler (or, technically, to the preprocessor). Such macros are typically used to activate/de-activate parts of the C++ code, as in this code fragment

  ```c++

  [...]

  #ifdef REFINEABLE
      std::cout << "I'm doing something because the REFINEABLE macro has been defined\n";
  #else
      std::cout << "The REFINEABLE macro hasn't been defined. Doing nothing!\n";
  #endif

  [...]

  ```

  Specifying `-DREFINEABLE` during the compilation is equivalent to having the statement

  ```c++

  [...]

  #define REFINEABLE

  [...]

  ```

  in the code.

To accommodate such additional flags, the `oomph_add_executable()` function also accepts the optional flags

- `CXX_OPTIONS`: this should be used for compiler flags (e.g. `-Wall`, `-g`). Note that these tend to be compiler specific (i.e. not every compiler will use `-g` to compile with debug information).
- `CXX_DEFINITIONS`: this should be used to define compiler macros such as `REFINEABLE`. Note that this keyword does not include the `-D` prefix; CMake will automatically prepend that for you (if this is what the compiler requires).

For example

```cmake
oomph_add_executable(
  NAME one_d_poisson
  SOURCES one_d_poisson.cc
  LIBRARIES oomph::poisson
  CXX_OPTIONS -Wall -Werror
  CXX_DEFINITIONS REFINEABLE)
```

will replicate the `g++` command shown above (but also link the code against `oomph-lib` and any associated third-party libraries).

> [!NOTE]
> Given that most compilers pass macros via the `-D` option, you could also add the flag `-DREFINEABLE` to the `CXX_OPTIONS` and omit the `CXX_DEFINITIONS`. It will achieve the same thing but may upset the CMake purists.

In the above example we've hard-coded the use of the `REFINEABLE` macro into the `CMakeLists.txt` file, so it will be applied for every build of that executable. What if you want to control its use from the command line when configuring your project? The temptation is to use the `-DCMAKE_CXX_FLAGS` command line option for CMake, as in
```bash
# Note: don't do this!
cmake -G Ninja -B build -Doomphlib_ROOT=/home/joe_cool/oomph_lib_install_dir -DCMAKE_CXX_FLAGS="-DREFINEABLE"
```
However, this will overwrite all the C++ compiler flags from the `oomph-lib` build (e.g. paranoia, range checking, debug vs. release mode, etc.). This is extremely dangerous and should be avoided. The proper way to handle this is to put the logic into the `CMakeLists.txt` file:
```bash

```


#### Customising targets using native CMake commands; hashed target names

If you are comfortable with CMake, you may wish to control the target properties of executables in a `CMakeLists.txt` file using native CMake commands. When doing this it is important to realise that the `NAME` specified in the call to `oomph_add_executable(...)` is not the name used by CMake itself. (The reason is technical: we create a modified, hashed version of the name which includes the path relative to the `demo_drivers` directory to avoid clashes between demo driver codes).

To use native CMake commands when customising the properties of an executable, you first have to obtain the hashed name that is used by CMake. For this purpose
we provide the `oomph_get_target_name(...)` function provided by `oomph-lib`:

```cmake
# Test definition
oomph_add_executable(NAME one_d_poisson [...])

# Get the hashed target name and store it in HASHED_TARGET_NAME
oomph_get_target_name(one_d_poisson HASHED_TARGET_NAME)

# Do something to the target, referring to the executable via
# the hashed name actually used by CMake. (Here we specify that the code
# contains C++20 syntax.)
set_target_properties(${HASHED_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
```

#### Using native CMake commands throughout

If you are comfortable with CMake and feel the `oomph_add_executable()` command does not provide the flexibility you require, you may wish to define the build of your own executable using the native CMake functions. If so, you will need to make sure that you add the compile definitions contained in the variable `OOMPH_COMPILE_DEFINITIONS` to the target -- this is what passes the required information about the `oomph-lib` installation to the compiler.

To illustrate this, the equivalent of

```cmake
oomph_add_executable(
  NAME one_d_poisson
  SOURCES one_d_poisson.cc
  LIBRARIES oomph::poisson
  CXX_OPTIONS -Wall -Werror
  CXX_DEFINITIONS REFINEABLE)
```

is the given by the following sequence of native CMake commands:

```cmake
# Specify the target name and source(s)
add_executable(one_d_poisson one_d_poisson.cc)

# Link against oomph-lib's poisson library
target_link_libraries(one_d_poisson PRIVATE oomph::poisson)

# Compiler options
target_compile_options(one_d_poisson PRIVATE -Wall -Werror)

# Preprocessor macros
target_compile_definitions(one_d_poisson PRIVATE REFINEABLE)

# Pull in the information about oomph-lib build
target_compile_definitions(one_d_poisson PRIVATE ${OOMPH_COMPILE_DEFINITIONS})
```

> [!NOTE]
> The `PRIVATE` keyword here tells CMake that the compile definitions only apply to the target one_d_poisson itself. Since this target is an executable, nothing else can depend on it, so using `PUBLIC` or `INTERFACE` would have no practical effect. However, you must specify one of `PRIVATE`, `PUBLIC`, or `INTERFACE`, or CMake will generate an error.

## Creating new demo driver directories

Next we consider how to create a new demo driver directory. This is mainly of interest for developers but, as already mentioned, you can also use this approach to create a directory in the `oomph-lib` directory tree to develop your own code.

As always, it is easiest to copy an existing directory and then modify it. Here we demonstrate how to create a new directory and how to develop the associated self-tests in the `demo_drivers/poisson` directory:

```bash
# Go the relevant directory
cd demo_drivers/poisson

# Make a new demo driver directory
mkdir my_one_d_poisson

# Copy across the source for the driver code:
cp one_d_poisson/one_d_poisson.cc my_one_d_poisson/

# ...and the CMakeLists.txt file
cp one_d_poisson/CMakeLists.txt my_one_d_poisson/

# ...and the validation script
cp one_d_poisson/validate.sh my_one_d_poisson/

# ...and the validata
cp -r one_d_poisson/validata my_one_d_poisson/
```

The `my_one_d_poisson` directory is now ready to be used but will, of course, just replicate what the original one did. Let's modify its content:

### Step 1: Rename and edit the driver code

Change the code:

```bash
# Go there
cd my_one_d_poisson

# Rename
mv one_d_poisson.cc my_fancy_new_poisson_code.cc

# Edit the code: turn it to something even more exciting!
emacs my_fancy_new_poisson_code.cc
```

### Step 2: Update the CMakeLists.txt file

Update the messages printed when you start and finish configuring the sub-directory:

```bash
# OLD
message(VERBOSE "Entered one_d_poisson subdirectory")

# NEW
message(VERBOSE "Entered my_one_d_poisson subdirectory")
```

```bash
# OLD
message(VERBOSE "Leaving one_d_poisson subdirectory")

# NEW
message(VERBOSE "Leaving my_one_d_poisson subdirectory")
```

Note that these messages will only be printed if you configure with `-DCMAKE_MESSAGE_LOG_LEVEL=VERBOSE`, i.e. you run

```bash
cmake -G Ninja -B build -DCMAKE_MESSAGE_LOG_LEVEL=VERBOSE
```

This can be helpful when building several demo drivers at the same time as you will be able to see which subprojects are entered.

Next, update the project name to the name of the enclosing directory. (You can, of course, call the project anything you want but this is a sensible default that must be adopted for any demo driver directories to be added to `oomph-lib`.)

```bash
# OLD
project(one_d_poisson C CXX Fortran)

# NEW
project(my_one_d_poisson C CXX Fortran)
```

The project name is followed by the languages the sources are written in. `oomph-lib` itself is written entirely in C++ but it also includes some third-party C and Fortran sources.

Update the name of the executable and the underlying sources, as well as the `oomph-lib` libraries required:

```bash
# OLD
oomph_add_executable(
  NAME one_d_poisson
  SOURCES one_d_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# NEW
oomph_add_executable(
  NAME my_fancy_new_poisson_code
  SOURCES my_fancy_new_poisson_code.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

[...]
```

### Step 3: Run and develop your own code

You can now build and run the code, using the steps explained earlier:

```bash
# Configure
cmake -G Ninja -B build

# Build in build directory
cd build/
ninja my_fancy_new_poisson_code

# Run it
./my_fancy_new_poisson_code
```

Of course, you'll be working on your code for a while until it does what you want it to do. Then it's time for the final two steps:

### Step 4: Update the validation script and the validata

At the moment, we haven't touched the `validate.sh` script yet, so running `ctest` will fail -- the script is likely to look for the name of the old execuable, etc. So update the `validate.sh` script and the associated validation data in the `validata` directory. The idea behind the script should (hopefully) be self-explanatory: it runs the code, collates some representative output data into a file which is then compared (allowing for small floating-point errors) against reference data in `validata`. You'll have to decide what to do!

To make sure the self-test is run properly when issuing `ctest` you'll also have to update the test specification in the `CMakeLists.txt` file:

```bash
# OLD
oomph_add_test(
  TEST_NAME poisson.one_d_poisson
  DEPENDS_ON one_d_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# NEW
oomph_add_test(
  TEST_NAME poisson.my_one_d_poisson # <-- this mimics the directory structure
  DEPENDS_ON my_one_d_poisson # <-- name of the executable required for the test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} # <-- the arguments to pass to validate.sh ('OOMPH_ROOT_DIR' is set when 'find_package(oomphlib)' is called; it is required by the validate.sh scripts to find fpdiff.py and validate_ok_count)
  TEST_FILES validate.sh validata  # <-- any files that should be available to the build directory for the test
)
```

### Step 5: Register the new directory so it is included in a top-level build

At the moment, the newly-created demo driver directory is not visited when the demo drivers are built (or self-tested) from the top level. So if you build the demo drivers in the `demo_drivers/poisson` directory like this:

```bash
cd demo_drivers/poisson

# Configure in verbose mode to see which directories are entered
cmake -G Ninja -B build -DCMAKE_MESSAGE_LOG_LEVEL=VERBOSE
```

you will see that the `my_one_d_poisson` directory was not entered, i.e. you will not see the lines:

```bash
--   Entered my_one_d_poisson subdirectory
--   Leaving my_one_d_poisson subdirectory
```

To ensure we enter the `my_one_d_poisson` directory, open `demo_drivers/poisson/CMakeLists.txt` and add `my_one_d_poisson` to the `SUBDIRS` variable:

```cmake
[...]

set(SUBDIRS
  poisson_with_singularity
  one_d_poisson
  my_one_d_poisson                    # <--- new line
  one_d_poisson_generic_only
  one_d_poisson_adapt
  one_d_poisson_hp_adapt
  two_d_poisson
  two_d_poisson_flux_bc
  two_d_poisson_flux_bc2
  two_d_poisson_flux_bc_adapt
  fish_poisson
  two_d_poisson_adapt
  two_d_poisson_hp_adapt
  fish_poisson2
  elastic_poisson
  eighth_sphere_poisson
  eighth_sphere_poisson_hp_adapt
  spectral)

[...]
```

If you now configure again

```bash
cd demo_drivers/poisson
cmake -G Ninja -B build -DCMAKE_MESSAGE_LOG_LEVEL=VERBOSE
```

you will see the following lines

```bash
--   Entered my_one_d_poisson subdirectory
--   Leaving my_one_d_poisson subdirectory
```

and doing

```bash
cd build
ctest
```

will now also run the self-tests in the newly-created demo driver. Yay!

## Creating a new library (or updating an existing one)

The `oomph-lib` libraries are defined in the `src/` directory and are separated by functionality: `generic/` contains the generic machinery; most other directories contain equation-specific elements, meshes, solvers, etc. If you have developed some new functionality that would benefit from being added as a new library, you simply create a new sub-directory in `src/`. It will typically contain header and source files. The latter need to be compiled and linked into a library. The headers and the library then need to be installed. As usual, the instructions that allow CMake to do this are encoded in the `CMakeLists.txt` file.

We illustrate this for `oomph-lib`'s `unsteady_heat` library. In a clean checkout from GitHub it contains the following files:

```bash
>>> ls -1
CMakeLists.txt
refineable_unsteady_heat_elements.cc
refineable_unsteady_heat_elements.h
Tunsteady_heat_elements.cc
Tunsteady_heat_elements.h
unsteady_heat_elements.cc
unsteady_heat_elements.h
unsteady_heat_flux_elements.h
```

Here's the `CMakeLists.txt` file:

```bash
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_heat subdirectory")

# Define the headers
set(HEADERS
    unsteady_heat_elements.h
    unsteady_heat_flux_elements.h
    refineable_unsteady_heat_elements.h
    Tunsteady_heat_elements.h)

# Define the sources
set(SOURCES
    unsteady_heat_elements.cc
    refineable_unsteady_heat_elements.cc
    Tunsteady_heat_elements.cc)

# Set the name of the library
set(LIBNAME unsteady_heat)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation.
# No need to touch this!
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
```

### Updating an existing library

To work on an existing library, simply edit the source/header files and/or add new ones. If you have added any new files, add them in the appropriate places in the `CMakeLists.txt` file.

To rebuild/install the library you return to the `oomph-lib` root directory and go through the standard configure/build/install sequence:

```bash
# Run this in the oomph-lib root directory:
cmake -G Ninja -B build
cmake --build build
cmake --install build
```

In fact, if you have already built and installed `oomph-lib` (using the directory `build/` in the `oomph-lib` root directory) you can simply do the following

```bash
# Run this in the oomph-lib root directory:

# Go to the build directory
cd build

# Reconfigure/build/install -- ninja is clever enough to pick up any
# new/modified files, including changed headers, sources and CMakeLists.txt
# files. The latter also makes it possible for ninja to detect newly added
# directories!
ninja install
```

### Adding a new library

To add a new library simply create a new directory in `oomph-lib`'s `src/` directory
and populate it with the required header and source files.
The listing of the `CMakeLists.txt` file for the `unsteady_heat` library shown above should suffice to understand what information is required. However, here's an empty template with a slightly more detailed annotation which illustrates what a `CMakeLists.txt` file has to specify to build a new library whose sources would live in the directory `src/funky_new_equation`. Things that need to be populated/changed are identified by `CHANGE_ME`.

```bash
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")

# CHANGE_ME: Message for verbose output
message(VERBOSE "Entered funky_new_equation subdirectory")

# CHANGE_ME: List your header files here; they get copied across into the
#            include directory when the library is installed.
# Define the headers
set(HEADERS
    # <header-file-1>
    # ...
    # <header-file-N>
)

# CHANGE_ME: List source files here; they get compiled and linked into a library
#            which is then copied across into the lib directory when the
#            library is installed
# Define the sources
set(SOURCES
    # <source-file-1>
    # ...
    # <source-file-N>
)

# CHANGE_ME: library name is usually the same as the name of the
#            sub-directory within src.
# Set the name of the library
set(LIBNAME funky_new_equation)

# CHANGE_ME: Add any other libraries that are required by this library here.
#            Most libraries will require the generic library to get access to
#            oomph-lib's basic objects; libraries that implement multi-physics
#            interactions may require additional libraries.
# Link against the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
# Boilerplate; no need to change this
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

# CHANGE_ME: Message for verbose output
message(VERBOSE "Leaving funky_new_equation subdirectory")
# ------------------------------------------------------------------------------
```

To make sure the library gets built automatically with the rest of `oomph-lib`,
simply add its sub-directory name to the `SUBDIRS` variable in `src/CMakeLists.txt`

```bash
[...]

set(SUBDIRS
    generic
    meshes
    poisson
    funky_new_equation       # <--- new line
    unsteady_heat

[...]
```

It is important to note that order matters. If the `funky_new_equation` library depends on the `generic` library, it should be listed after the `generic` library, so that the `generic` library is already defined when we try to define the `funky_new_equation` library.

When you're done, return to the `oomph-lib` root directory and go through the
usual configure, build and install sequence

```bash
# Run this in the oomph-lib root directory
cmake -G Ninja -B build
cmake --build build
cmake --install build
```

which should be quick because CMake will realise that only the new library needs to be built.

The alternative, already described above,

```bash
# Run this is in the oomph-lib root directory

# Go to build directory
cd build

# Configure/build/install anything that's new
ninja install
```

should also work if `oomph-lib` has already been configured, using a build directory called `build` in the `oomph-lib` root directory.

## The `doc` directory

`oomph-lib`'s documentation is contained in the directory `doc`. Its directory structure is roughly equivalent to that of the `demo_drivers/` directory which contains the driver codes that are explained in the various tutorials.

To build the documentation either specify the flag `--build-doc` to `oomph_build.py` or build the documentation manually, following the usual pattern:

```bash
cd doc
cmake -G Ninja -B build
cmake --build build
```

There's no need for an installation.

Note that building the documentation takes a long time and simply duplicates the documentation already available on the `oomph-lib` webpage. It is therefore mainly of interest to developers or other volunteers who want to provide new (or fix broken) tutorials. Once the build process is complete you can navigate your local copy of the `oomph-lib` webpages starting from `doc/html/index.html`.

The way the documentation is built/installed is slightly different from how it's done when building the library itself, so simply deleting the `build` directory won't get rid of everything that's been created. To clean up do the following

```bash
cd doc/build
ninja clean
cd ..
rm -rf build
```

To add new tutorials to the `doc` directory, follow the same steps used for adding new demo drivers. Detailed instructions on what's required and how the documentation is generated are available on the `oomph-lib` webpage.

## Cheat sheet: Autotools vs. CMake

Earlier releases of `oomph-lib` used Autotools to build/install the library. To ease the conversion to the new way of doing things, the list below provides a comparison of key operations. Note that we also provide a script, `scripts/oomph_generate_cmake_script.py`, that (attempts to) converts an `oomph-lib` Autotools `Makefile.am` into a corresponding `CMakeLists.txt`. It may need a bit work/checking but it tends to do most of the work.

In all the examples shown below we assume that `$oomph_home_dir` is the `oomph-lib` home directory:

### Build and install the entire library using a helper script

#### Autotools

```bash
# Go to oomph-lib home directory
cd $oomph_home_dir

# Run the autogen script which prompts the user for configure options
./autogen.sh
```

#### CMake

```bash
# Go to oomph-lib home directory
cd $oomph_home_dir

# Run the oomph_build.py script, configured via command line
# options; no interactivity. Use ./oomph_build.py --help to
# find out what's available. To do a default build (No MPI,
# build all third-party libraries, build in Release mode, do a
# local installation) do
python3 ./oomph_build.py
```

### Alternative: Build and install the entire library step by step

#### Autotools

```bash
# Go to home directory
cd $oomph_home_dir

# Run manual config/make/install sequence, installing in specfied
# directory. Note that with our autotools based build process, third-party
# libraries  were an integral part of the distribution and were thus built
# automatically.
./configure --prefix /home/joe_cool/oomph_lib_install_dir
make
make install
```

#### CMake

```bash
# Go to home directory
cd $oomph_home_dir

# Build third-party libraries and install them locally; no separate
# "install" step required.
cd external_distributions
cmake -G Ninja -B build
cmake --build build

# Run manual config/build/install sequence for oomph-lib itself
# specify the location of the third party libraries by parsing
# the file generated automatically during the previous step
cd $oomph_home_dir
cmake -G Ninja -B build $(cat external_distributions/install/cmake_flags_for_oomph_lib.txt)
cmake --build build
cmake --install build
```

### Run all self-tests

#### Autotools

```bash
# Go to home directory
cd $oomph_home_dir

# Run the tests
make check
```

#### CMake

```bash
# Go to demo_drivers directory
cd $oomph_home_dir/demo_drivers

# Configure
cmake -G Ninja -B build

# Run the tests (using 8 cores)
cd build
ctest -j 8
```

### Run all self-tests in a given demo driver directory

#### Autotools

```bash
# Go to demo driver directory
cd $oomph_home_dir/demo_drivers/one_d_poisson

# Run the test
make check
```

#### CMake

```bash
# Go to demo driver directory
cd $oomph_home_dir/demo_drivers/one_d_poisson

# Configure
cmake -G Ninja -B build

# Run the test
cd build
ctest
```

### Debug demo driver code

#### Autotools

```bash
# Go to demo driver directory
cd $oomph_home_dir/demo_drivers/one_d_poisson

# Run the test
make check

# Edit
emacs one_d_poisson.cc

# Run again (make realises that code needs to be recompiled)
make check

# Keep going until it works... (of course you may do the
# compilation/rerun in your editor/IDE; see below)
```

#### CMake

```bash
# Go to demo driver directory
cd $oomph_home_dir/demo_drivers/one_d_poisson

# Configure (note that switching on the Debug build only helps if the
# library itself was compiled with that option too; by default it is
# built in Release mode with full optimisation)
cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=Debug

# Run the test (note that the executable
# got built in the build directory!)
cd build
ctest

# Edit source code
emacs ../one_d_poisson.cc

# Rebuild executable(s)
ninja

# ...or just rebuild a specific one
ninja one_d_poisson

# If you want to see compiler errors
# in detail do
ninja --verbose one_d_poisson

# Now rerun the test
ctest

# Keep going until it works... (of course you may do the
# compilation/rerun in your editor/IDE; see below)
```

### Edit a file in the library and rebuild

#### Autotools

```bash
# Go to relevant src directory
cd $oomph_home_dir/src/poisson

# Edit the relevant file
emacs poisson_elements.h

# Recompile the library
make

# You may go through this repeatedly
# until the code compiles...

# EMACS/IDE ASIDE: If you do this in emacs where you can trigger
# the compilation with `[F4]` (or some other user-defined key) the
# command to be spawned is "make" which emacs enters by default
# (you can edit this); if the compilation fails, emacs allows
# you to jump to the source line where the compilation error occurs
# using `[F12]` (on my setup; you can customise this), so the typical
# sequence is:
# - edit code
# - save
# - `[F4]` (spawns "make")
# - [RETURN]
# - `[F12]` (jump to error; edit)
# - `[F12]` (next error; edit)
# - ...
# - save
# - `[F4]` (spawns "make")
# - [RETURN]
# etc.

# Install the library
make install

# If the change also affects other
# dependent libraries (e.g. if you
# edited something in src/generic
# or a library that is used by others,
# e.g. for multi-physics elements,
# rebuild those libraries too. Just to be
# on the safe side, do all of them (make
# knows about the dependencies and will
# only rebuild those that are affected)
cd $oomph_home_dir/src
make
make install
```

#### CMake

```bash
# Go to relevant src directory
cd $oomph_home_dir/src/poisson

# Edit the relevant file
emacs poisson_elements.h

# Go back to the home directory and
# rebuild (CMake knows what needs to be
# recompiled)
cd $oomph_home_dir
cmake --build build

# You may go through this repeatedly
# until the code compiles...

# EMACS/IDE ASIDE: If you do this in emacs where you can trigger
# the compilation with `[F4]` (or some other user-defined key),
# emacs provides the default command "make". Edit this to
#
#   cd $oomph_home_dir; cmake --build build
#
# then hit return. If the compilation fails, emacs allows you to
# jump to the source line where the compilation error occurs using
# `[F12]` (on my setup; you can customise this), so the typical
# sequence is:
# - edit code
# - save
# - `[F4]` (spawns "make")
# - change make to
#
#     cd $oomph_home_dir; cmake --build build
#
# - [RETURN]
# - `[F12]` (jump to error; edit)
# - `[F12]` (next error; edit)
# - ...
# - save
# - `[F4]` (emacs now remembers the previous
#         command)
# - [RETURN]
# etc.

# Install the library
cd $oomph_home_dir
cmake --install build
```

## Dos and Don'ts

- Do not (re-)define the `PARANOID` or `RANGE_CHECKING` macros in your driver code because this could lead to inconsistencies between the header files and the installed libraries, potentially leading to nasty and hard-to-diagnose segmentation faults. The macros used when building the library are automatically imported into your driver code, so their status is available. Just don't assign them yourself!
- Do (re-)build the library with `PARANOID` or `RANGE_CHECKING` enabled if you develop any new machinery or work on a new driver code. The code will run more slowly but the warnings will save you years of your life!
- Do not use the `-DCMAKE_CXX_FLAGS` flag during the CMake configure stage to pass C++ compiler macros to the build process. This applies to the `oomph-lib` build and to the configuration of your own stand-alone projects. If you use this flag in the latter, you will overwrite flags that were defined when configuring `oomph-lib`, with potentially disastrous consequences. Please refer to [this section](#how-to-add-additional-compiler-macros-to-oomph-lib-and-to-stand-alone-driver-codes) for instructions on how to pass compiler macros to the `oomph-lib` build and [this repository](https://github.com/oomph-lib/stand_alone_oomph-lib_user_code) for an illustration how to do this for stand-alone driver codes.

## FAQ

### When configuring my driver code CMake can't find the package configuration file

If you get an error message like this
```bash
CMake Error at CMakeLists.txt:86 (find_package):
  Could not find a package configuration file provided by "oomphlib" with any
  of the following names:

    oomphlibConfig.cmake
    oomphlib-config.cmake

  Add the installation prefix of "oomphlib" to CMAKE_PREFIX_PATH or set
  "oomphlib_DIR" to a directory containing one of the above files.  If
  "oomphlib" provides a separate development package or SDK, be sure it has
  been installed.
```
when configuring your driver code then you've either specified the wrong 
installation directory or you haven't specified one at all when you should 
have. In the latter case you may have used the command
```bash
cmake -G Ninja -B build
```
when you should have done 
```bash
cmake -G Ninja -B build -Doomphlib_ROOT=/home/joe_cool/oomph-lib_playground/oomph_lib_installation
```
say. (This assumes that you've installed `oomph-lib` in 
`/home/joe_cool/oomph-lib_playground/oomph_lib_installation`). Please consult 
the section [Linking a stand-alone project to `oomph-lib`](#linking-a-stand-alone-project-to-oomph-lib) for detailed instructions.

### Are there any complete worked examples of the build process?
Yes, there are! The script
```bash
scripts/how_to_build_and_test_example.bash
``` 
was written in parallel to this documentation. It demonstrates many variants of the entire build process, starting with downloading the sources from GitHub, installing the library (using various methods and settings), running driver codes within `oomph-lib` and linking stand-alone projects to it. 

### I'm used to the Autotools-based version of `oomph-lib`; what happened to the `user_src` directory?
When rewriting the build system to CMake, we retained the `user_drivers` directory to allow users to build their own code within the overall `oomph-lib` framework (but outside the `demo_drivers` directory; it doesn't belong there). The `user_src` directory was initially provided as test-bed within which users could develop their own libraries. The change to GitHub and CMake makes this unnecessary: if you want to provide your own library within the overall `oomph-lib` framework, simply create a fork of `oomph-lib` on GitHub, clone this onto your computer and create a new branch, within which you can create a new library, placing the code in a new directory, `src/joe_cools_new_library`, say. See the section [Adding your own library](#adding-a-new-library)
for instructions. If it all works and is likely to be useful for others, send us a pull request and we'll consider including it into `oomph-lib`'s development branch.

### I'm used to the Autotools-based version of `oomph-lib`; what happened to the `bin` directory?
We renamed it to `scripts` because... Ahem; can't remember.


## Additional information for developers

### Use symbolic links for header files

Use the `--oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON` flag when building the library with `oomph_lib.py` 

```bash
./oomph_build.py [...] --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON
```

or use `-E create_symlink` when configuring the `oomph-lib` build directly

```bash
cmake -G Ninja -B build -E create_symlink
```

With these flags, the header files in the `install` directories (usually copied from the `src` directory) are replaced by symbolic links to the original files in the `src` directory. This means that, if you build a library and are alerted to an error in a header file you don't accidentally edit a copy (that is then promptly overwritten when you reinstall the library).

### Paranoia and range checking are deemed to be incompatible with Release mode

As mentioned before, we strongly recommend building `oomph-lib` with the `PARANOID` and `RANGE_CHECKING` flags when developing new code. (Recall that these are enabled with `-DOOMPH_ENABLE_PARANOID=ON` and `-DOOMPH_ENABLE_RANGE_CHECKING=ON` when configuring the `oomph-lib` build with CMake). These compiler macros activate a large number of internal sanity checks that help you trap errors. The (hopefully) helpful messages that are issued when an error is encountered should facilitate debugging the new code. However, since the sanity checks add to the runtime, our build machinery does not allow the flags to be enabled when `oomph-lib` is built in `Release` mode (which is the default setting, if the mode is not specified explicitly). Note that using range checking and paranoia without having compiled the code in `Debug` mode (which adds the `-g` flag to the compiler) is unlikely to be helpful anyway. The executable may tell you that a function was called with illegal arguments, say, but without being able to backtrack where this function was called from this information is unlikely to be helpful! 

### How to add additional compiler macros to oomph-lib (and to stand-alone driver codes)
`oomph-lib` uses the compiler macros `PARANOID` and `RANGE_CHECKING` to isolate optional code that facilitate debugging; the relevant code is only compiled if the macros are defined. The same technique can be used to isolate "work in progress" code that shouldn't (yet) be used by the general public. (Of course, there's a separate question if "work in progress" code should even be committed to the repository. The answer is generally no, and code reviews should at least question the wisdom of this!). Assuming the relevant macro is called `USE_NEW_CODE`, use the`-DOOMPH_EXTRA_COMPILE_DEFINES` flag when configuring the `oomph-lib` build, i.e. 
```bash
cmake -G Ninja -B build [...] -DOOMPH_EXTRA_COMPILE_DEFINES="USE_NEW_CODE" 
```
or use the `--oomph-OOMPH_EXTRA_COMPILE_DEFINES` flag when building `oomph-lib` using the build script:
```bash
./oomph_build.py [...] --oomph-OOMPH_EXTRA_COMPILE_DEFINES="USE_NEW_CODE" 
```
Note that the `-D` is omitted. You can pass multiple flags in the same quoted string, separating the definitions by spaces, e.g. 
```bash
cmake -G Ninja -B build [...] -DOOMPH_EXTRA_COMPILE_DEFINES="USE_NEW_CODE USE_OTHER_NEW_CODE" 
```
or 
```bash
./oomph_build.py [...] --oomph-OOMPH_EXTRA_COMPILE_DEFINES="USE_NEW_CODE USE_OTHER_NEW_CODE" 
```
when using the build script.


Please refer to [this repository](https://github.com/oomph-lib/stand_alone_oomph-lib_user_code) for an illustration of how to do this for stand-alone driver codes that use `oomph-lib`. The `CMakeLists.txt` file in that repository illustrates various ways of doing this, but whatever you do, do not use the `-DCMAKE_CXX_FLAGS` flag when configuring your project.


### Creating robust `validata` for self tests

It is important to make sure that the `validata` is robust. While the `scripts/fpdiff.py`
allows for floating point tolerances (both absolute and relative) it
cannot cope with data appearing in a different order. This can be very
confusing to debug because the code then produces the correct results
and plotting packages will display them correctly, but yet the tests fail.

Such problems are usually caused by the use of wildcards
in the validations scripts, the comparison of data that is stored
in certain STL containers, or the use of non-deterministic algorithms
(especially the mesh generation with `triangle`).

#### A self-test fails even though the output files produced by the code are correct

Self-tests are performed by the `validate.sh` shell script, which typically runs the executable and
concatenates selected output files to a single file whose contents
are compared against the reference file in the `validata` directory.
While it is tempting to write

```bash
cat RESLT* /soln0.dat > results_file.dat
```

it is important to realise that the order in which the files are concatenated is machine- and/or operating-system dependent. If the above command is run in a directory with the following structure

```bash
|-- RESLT
|   |-- soln0.dat
|   `-- trace.dat
|-- RESLT_elastic
|   |-- soln0.dat
|-- trace.dat
```

some operating systems will expand the command to

```bash
cat RESLT/soln0.dat RESLT_elastic/soln0.dat > results_file.dat
```

while others will execute

```bash
cat RESLT_elastic/soln0.dat RESLT/soln0.dat > results_file.dat
```

In this case the self-test will report a failure, even though the solution files are correct. The `validate.sh` scripts should therefore not contain any wildcards.

#### Handling non-deterministic output

Similar problems can arise if the validation data includes
data that is stored in certain STL containers such as sets.
The order in which items are stored in such containers may
vary from machine to machine and from compiler to compiler.
If such data is to be included in a self-test the data should be
sorted first, based on a user-controllable sorting criterion.

#### Careful with driver codes that use `triangle` to generate meshes

We provide a wrapper to the third-party `triangle` code to generate unstructured 2D meshes. When using these it is important to realise that the meshes generated can change, depending on the machine the code is run on and even the optimisation level of the code. `validata` for such codes must therefore use intergral data such as
norms of the solution, rather than element-by-element output.

### How to update third-party libraries to later versions

The third-party libraries in `external_distributions` are pulled in from their respective GitHub repositories. The version of the libraries is typically encoded in the `GIT_TAG` variable name in the files `external_distributions/cmake/OomphGetExternal<library_name>.cmake` file. An upgrade should (usually) just require a change to that tag, though if a more recent version of the library requires different build steps the relevant `*.cmake` file will, of course, have to be modified further. Make sure that the relevant self-tests still produce the same results after upgrading to a new version of a third-party library!

## Appendix

### CMake resources

For more information on CMake you may wish to consult the following resources:

- The excellently-written "Professional CMake: A Practical Guide" by Scott Craig.
- The [Awesome CMake](https://github.com/onqtam/awesome-cmake) repository.
- [An Introduction to Modern CMake](https://cliutils.gitlab.io/modern-cmake/).
- ...and the list goes on (so add more!).

### Building CMake

If your computer doesn't have a sufficiently up-to-date version of CMake you will have to build it yourself. Here's how to do this:

#### Ubuntu

For Linux, CMake provides an installer script to help you download and install CMake:

```bash
# Pick version and installation location
CMAKE_VERSION=3.24.4
CMAKE_DEST_DIR=~/.cmake-${CMAKE_VERSION}

# Download installer script and run
mkdir ${CMAKE_DEST_DIR}
wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh -O cmake.sh
bash cmake.sh --prefix=${CMAKE_DEST_DIR} --exclude-subdir

# Update environment
export PATH=${CMAKE_DEST_DIR}/bin:$PATH
```

### macOS

To accommodate both Intel-based and Arm-based Macs, CMake provides a "universal binary". To download CMake via this approach, use the following commands:

```bash
# Pick version and download location
CMAKE_VERSION=3.24.4
CMAKE_DEST_DIR=~/.cmake-${CMAKE_VERSION}

# Download package and move to desired location
wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-macos-universal.tar.gz
tar xvfz cmake-${CMAKE_VERSION}-macos-universal.tar.gz
mv cmake-${CMAKE_VERSION}-macos-universal ${CMAKE_DEST_DIR}
rm -f cmake-${CMAKE_VERSION}-macos-universal.tar.gz

# Update environment
export PATH=~/${CMAKE_DEST_DIR}/CMake.app/Contents/bin:${PATH}
```

> [!TIP]
> To make the changes to the `$PATH` variable permanent, add the `export PATH` commands to the end of your shell start-up script, e.g. `.bashrc` or `.zshrc`.
