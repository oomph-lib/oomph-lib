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

- [Description](#description)
  - [Compatibility](#compatibility)
- [Prerequisites](#prerequisites)
  - [CMake](#cmake)
  - [Ninja](#ninja)
- [Building, installing and uninstalling `oomph-lib`](#building-installing-and-uninstalling-oomph-lib)
  - [Option 1: Default installation](#option-1-default-installation)
  - [Option 2: Specifying a custom installation location](#option-2-specifying-a-custom-installation-location)
  - [Option 3: Installation with root privileges](#option-3-installation-with-root-privileges)
- [Testing the installation](#testing-the-installation)
  - [Testing the entire installation](#testing-the-entire-installation)
  - [Testing specific demo driver codes](#testing-specific-demo-driver-codes)
- [How to run demo driver codes (and how to modify them when working on coding exercises)](#how-to-run-demo-driver-codes-and-how-to-modify-them-when-working-on-coding-exercises)
- [Linking a stand-alone project to `oomph-lib`](#linking-a-stand-alone-project-to-oomph-lib)
- [Creating new demo driver directories](#creating-new-demo-driver-directories)
  - [Step 1: Rename and edit the driver code](#step-1-rename-and-edit-the-driver-code)
  - [Step 2: Update the CMakeLists.txt file](#step-2-update-the-cmakeliststxt-file)
  - [Step 3: Run and develop your own code](#step-3-run-and-develop-your-own-code)
  - [Step 4: Update the validation script and the validata](#step-4-update-the-validation-script-and-the-validata)
  - [Step 5: Register the new directory so it is included in a top-level build](#step-5-register-the-new-directory-so-it-is-included-in-a-top-level-build)
- [Creating a new library (or updating an existing one)](#creating-a-new-library-or-updating-an-existing-one)
  - [Updating an existing library](#updating-an-existing-library)
  - [Adding a new library](#adding-a-new-library)
- [Make a new doc directory](#make-a-new-doc-directory)
- [Build options](#build-options)
- [CMake Presets](#cmake-presets)
  - [`CMakePresets.json`](#cmakepresetsjson)
  - [`CMakeUserPresets.json` example](#cmakeuserpresetsjson-example)
- [Advanced testing](#advanced-testing)****
  - [Filtering by regex](#filtering-by-regex)
  - [Disabling a test](#disabling-a-test)
- [Customising driver codes](#customising-driver-codes)
  - [Customising targets](#customising-targets)
- [To be documented](#to-be-documented)
- [A deeper dive into the build system](#a-deeper-dive-into-the-build-system)
- [Helpful CMake resources](#helpful-cmake-resources)
- [Additional instructions](#additional-instructions)
  - [Building CMake](#building-cmake)
    - [Ubuntu](#ubuntu)
    - [macOS](#macos)
- [Community](#community)

## Description

The [`oomph-lib` homepage](http://www.oomph-lib.org) provides much more detail on
installation instructions, tutorials, and licencing information. Provided you
have downloaded a distribution that contains the documentation and you have the
required tools (mainly `doxygen`; get it from [here](http://www.doxygen.org))
available on your machine, the installation procedure will create a local copy
of the `oomph-lib` webpages and the entire online documentation in the `doc`
directory. In particular, `doc/html/index.html` is a local copy of the `oomph-lib`
homepage.

To learn more about contributing to `oomph-lib`, please see
[`CONTRIBUTING.md`](CONTRIBUTING.md) which contains a detailed description of
the workflow.

### Compatibility

Operating system | Support provided?
-----------------|------------------
Ubuntu           | Yes
macOS            | Yes
Windows          | No

## Prerequisites

`oomph-lib` uses CMake and Ninja to build and install the project.
Make sure you have sufficiently recent versions of these programs installed on your computer:

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

If your package manager does not provide a recent enough version of `cmake`, you will need to build it from source. You can find instructions on how to do this [below](#building-cmake) and on the [CMake website](https://cmake.org/install/).

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

## Building, installing and uninstalling `oomph-lib`

### Option 1: Default installation

To configure, build and install the library, `cd` into the root directory of the cloned `oomph-lib` project (checked out from GitHub) and run the following commands:

```bash
# Configure and generate the build system; -B specifies the
# build directory (here "build")
cmake -G Ninja -B build

# Build the oomph-lib libraries (i.e. compile the sources
# and turn them into libraries); specify the directory
# that was created at the configure stage (here "build")
cmake --build build

# Install, again specifying the build directory created above
# (here "build")
cmake --install build
```

Note that, by default, the last step installs the headers and generated library files into a directory named `install/` in the `oomph-lib` root directory. This is done deliberately to make sure that even users who do not have root privileges on their computers can install the library.

To uninstall the library simply delete the `build/` and `install/` directories in the `oomph-lib` root directory.

Alternatively, you can use the provided Python build script (`oomph_build.py`) to perform the above configure, build, and install steps in one go. This script automates the process and uses the same defaults (creating a `build/` directory and installing to an `install/` subdirectory). New users may find it convenient to avoid typing multiple commands. (See the section [Building with `oomph_build.py`](#building-with-oomph_buildpy) below for more details.)

### Option 2: Specifying a custom installation location

By default, `oomph-lib` will be installed to the `install/` subdirectory of the root `oomph-lib` folder. To specify a different installation location, pass `-DCMAKE_INSTALL_PREFIX=<install-location>` to `cmake` during the configure step. Note that `<install-location>` must be an **absolute** path. Here's an example where we specify different build and install directories:

```bash
# Configure and build (write build information into the
# directory "joe_build") and use the install prefix to
# declare where the library should ultimately be installed.
cmake -G Ninja -B joe_build -DCMAKE_INSTALL_PREFIX=/home/joe_cool/oomph_lib_install_dir

# Build the library; specify the directory created during the
# configure stage (here "joe_build")
cmake --build joe_build

# ...and install it; again specify the build directory created
# above (here "joe_build"). The library will now be installed
# in /home/joe_cool/oomph_lib_install_dir
cmake --install joe_build
```

To make sure that the library can be found when building driver codes (see below) it is easiest to set the `oomphlib_ROOT` environment variable:

```bash
# Add the custom installation directory to the path
export oomphlib_ROOT=/home/joe_cool/oomph_lib_install_dir
```

This is, of course, only a temporary assigment for the current shell; add the command to your `.bashrc` file (or `.zshrc` file, if using Zsh) to reassign it automatically for every session. Alternatively, you can specify the location of the install directory with the `-Doomphlib_ROOT` flag when configuring the driver codes; see [below](#testing-the-entire-installation).

To uninstall the library, delete the `joe_build` directory in the `oomph-lib` root directory and the directory `/home/joe_cool/oomph_lib_install_dir`. Also remove the
addition of the installation directory to the `PATH` from your `.bashrc` or `.zshrc` file (if you added it there).

**Note:** The `oomph_build.py` script described in [Building with `oomph_build.py`](#building-with-oomph_buildpy) can also handle custom installation paths. For example, instead of the multi-step process above, you could run the Python build script with an option to specify the install location and achieve the same result.

### Option 3: Installation with root privileges

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

Alternatively, the `oomph_build.py` script can facilitate a system-wide installation. You can either (i) pass `--root-OOMPH_ALLOW_INSTALL_AS_SUPERUSER=ON` to the script or (ii) set the installation path explicitly using the `--root-CMAKE_INSTALL_PREFIX` option as `/usr/local` (and enable the necessary CMake flags for a root install). You would still need to run the install step with superuser privileges (the script will prompt or instruct you to do so). See [Building with `oomph_build.py`](#building-with-oomph_buildpy) for more information.

## Building with `oomph_build.py`

This section explains how to use the `oomph_build.py` script - a Python helper script located in the root of the repository – to configure, build, and install `oomph-lib` in an automated way.

### When to use the build script

The `oomph_build.py` script is provided as a convenient one-stop solution for building and installing `oomph-lib`. It essentially wraps the manual CMake steps (configuration, build, and installation) into a single command, with sensible defaults and additional helper options. New users are encouraged to use this script for a hassle-free installation, as it reduces the chance of missing a step or incorrectly specifying an option. Advanced users or those integrating `oomph-lib` into other build systems may still prefer the direct CMake approach (as in Options 1–3 above) for finer control, but the script covers the common use cases.

**Prerequisites:** In addition to the standard requirements (CMake and Ninja, as described above), you will need a working Python 3 installation to run `oomph_build.py`. Ensure that `python3` is available in your `PATH`. If you plan to use the script's documentation-building feature (described below), make sure `doxygen` is installed as well.

### Using the script

To use the build script, change into the top-level `oomph-lib` directory (the same directory that contains `oomph_build.py`) and run it with Python. For example:

```python
python3 oomph_build.py
```

By default, this will configure the project, compile all the libraries, and install them to the local `install/` directory (within the `oomph-lib` repository folder). The script will create a build directory (named `build/` by default) if one does not exist, or reuse it if it does. It combines the equivalent of the `cmake -B build`, `cmake --build build`, and `cmake --install build` commands into one automated sequence.

After running `oomph_build.py`, if it completes without errors, `oomph-lib` should be built and installed in the specified location. You can then proceed to test the installation as described in the next section.

> **Tip:** On Unix-like systems, you can make the script directly executable by running `chmod +x oomph_build.py`, then invoke it with `./oomph_build.py`. Using `python3 ...` or making it executable are equivalent as long as it's run with Python 3.

### Command-Line options

The `oomph_build.py` script supports several command-line options to customize its behavior. You can see a summary by running

```bash
python3 oomph_build.py --help
```

Below is a comprehensive list of the options and their purposes:

- **`--build-doc`**: By default, the script does not build the documentation (to save time and space). If you want to generate the full HTML (and PDF) documentation for oomph-lib, include the `--build-doc` flag. This will invoke the documentation build (using Doxygen and LaTeX) as part of the process. *Note:* Building documentation can be time-consuming and requires Doxygen (and a LaTeX distribution for PDFs) to be installed. If these tools are missing, the doc build will fail – in that case, either install the necessary tools or omit this option.

- **`--root-CMAKE_INSTALL_PREFIX`**: Use this option when you intend to provide a custom installation location.
- **`--root-OOMPH_ALLOW_INSTALL_AS_SUPERUSER`**: Use this option when you intend to install `oomph-lib` system-wide (e.g. to `/usr/local/`) using root privileges. This flag tells the script to configure the installation prefix to the system's default location (`/usr/local`) and to set any required internal CMake switches (such as`OOMPH_ALLOW_INSTALL_AS_SUPERUSER=ON`). (If you also specify `--root-CMAKE_INSTALL_PREFIX`, this flag will have no effect.)
  > **Important:** When using this option, you will need to run the installation step with administrative privileges. The script will attempt to perform the installation step with sudo if possible, or it will remind you to re-run the script as root for the install phase. It’s generally recommended to run oomph_build.py `--root-CMAKE_INSTALL_PREFIX` under a normal user for the build, and let it prompt for a password or instruct you for the install, rather than running the entire build as root. (Building as a non-root user helps avoid permission issues in the build directory.)
- **`--wipe-tpl`**, **`--wipe-root`**, **`--wipe-doc`**: These options tell the script to remove the specified build/installation directories that would be written to when building the third-party libraries, root project, and documentation, respectively. Use the `--wipe-*` flags if you want a completely clean rebuild. For example

  ```bash
  # Wipe the default build/installation directories used when building
  # the third-party libraries, root project, and documentation
  python3 oomph_build.py --wipe-tpl --wipe-root --wipe-doc
  ```

  will delete the current `build/` directory and the `install/` directory (if they exist) before configuring a fresh build. `--wipe-doc` will remove the `doc/build/` directory.. These options are useful if you suspect a previous build is causing issues or if you want to reclaim space and rebuild from scratch. Warning: Wiping will permanently delete those directories (and all compiled files or installed files therein), so use these flags with caution.
- **`--reuse-tpl`**, **`--reuse-root`**, **`--reuse-doc`**: These flags are the counterpart to the `--wipe-*` options, explicitly instructing the script to reuse existing directories. Using `--reuse-*` flag means the script will allow the existing build/installation directories and will instead update or rebuild files as needed. These options can save time on iterative builds, but be mindful that reusing directories might carry over old artifacts - if you run into strange problems, a wipe might be needed.

In summary, by default the script uses a safe approach (an error is thrown if the build/installation directories already exist), requiring you to re-run `oomph_build.py` with `--wipe-*` or `--reuse-*` flags. If a completely fresh build environment is desired, use the `--wipe-*` options.

### Running the script: example and output

When you run `oomph_build.py`, it will print messages indicating its progress through the build stages. For example, you will see output from the CMake configure stage first. This includes messages about detecting compilers, configuring third-party libraries, and generating build files. If the configuration is successful, you should see a message like

```bash
-- Build files have been written to: .../build
```

If there are errors (e.g., a required tool is not found or a configuration check failed), CMake will report them; the script will stop and propagate that error.

Next, the script invokes the build stage (equivalent to running Ninja). If you specified the `--verbose` flag you will see compiler output as the `oomph-lib` libraries compile. This stage may take some time, especially on first build. If it succeeds, all necessary library files will be produced in the build directory.

Finally, the script runs the installation stage. You'll see messages about installing libraries and header files to the target directory.

## Testing the installation

### Testing the entire installation

`oomph-lib` comes with an extensive list of well-documented example driver codes situated in the `demo_drivers/` directory. The driver codes in these folders are also used to validate the library.

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

### Testing specific demo driver codes

It is also possible to test a single driver code (or rather all the demo drivers below a specific subdirectory of `demo_drivers`):

```bash
# Go to a specific sub-directory in demo_drivers
cd demo_drivers/poisson/one_d_poisson/

# Configure (specify directory in which tests will be run)
cmake -G Ninja -B build

# Go into the new build directory
cd build

# Run test(s)
ctest
```

## How to run demo driver codes (and how to modify them when working on coding exercises)

The best way to get started with `oomph-lib` is to explore some of the demo-driver codes. Many of these codes are explained in great detail in the associated [tutorials](https://oomph-lib.github.io/oomph-lib/doc/example_code_list/html/index.html) which typically end with a few exercises that encourage you to modify the code.

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

The only thing to look out for the current task is the name of the executable.

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

It is not necessary to reconfigure the build directory, so we can build and execute the new driver code straightaway:

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
(global-set-key [f4] 'compile)
(define-key global-map [f12] 'next-error)
```

With these keybindings, the [F4] function key will do a `make -k` in the current directory (you can edit the command line that appears at the bottom of the emacs window to add the name of a specific executable). If any errors occur during the compilation, [F12] will go from error to error in the source file(s). Very useful!

The fact that CMake separates the source and build directories means this won't quite work any more. However, assuming you adopt the common convention of calling your build directory "build", adding this

```bash
(setq compile-command "cd build; ninja")
```

to your `.emacs` file will produce equivalent behaviour. You can now edit the source code in its directory; [F4] will then compile in the build directory (again, you can specify the name of a specific executable by editing the command line that appears at the bottom of the emacs window); [F12] will then go through errors in the source file(s).

## Linking a stand-alone project to `oomph-lib`

Developing your own code in an existing demo driver directory is a quick-and-dirty way to get started, especially since you are most likely to start your work by modifying an existing driver code. However, long-term this is not a sensible solution. One slightly more attractive alternative is to create a new directory, just for your code, in the `demo_drivers` directory; described further below. This approach has the advantage of not interfering with existing `oomph-lib` driver codes and the associated test machinery. However, your code isn't really a demo driver so it should really live somewhere else.

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
# Good practice, tell us what you're doing (if asked to)
message(VERBOSE "Entered mesh_gluing subdirectory")

# Specify minimum cmake version; die if you can't find it
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

# Name of the project, followed by languages used
project(mesh_gluing C CXX Fortran)

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
# environment variable.
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

# Say bye bye (if asked to)
message(VERBOSE "Leaving mesh_gluing subdirectory")
# ------------------------------------------------------------------------------
```

Note that this directory is completely unconnected to the `oomph-lib` directory. To be able to find the `oomph-lib` installation directory, if it is not in one of the system-wide standard locations such as `/usr/local`, it must be declared.

This is most easily done by using the `PATH` environment variable, as discussed [above](#option-2-specifying-a-custom-installation-location). Note that this is necessary even if `oomph-lib` is installed to its default installation directory, `install/`, in the `oomph-lib` root directory.

It is also possible to specify the location of the install directory in the `CMakeLists.txt` file, or to provide hints where to search for it. Please read the [CMake documentation of the `find_package()` function](https://cmake.org/cmake/help/latest/command/find_package.html?highlight=find_package).

Anyway, assuming `oomph-lib` was installed to `/home/joe_user/oomph_lib_install` the driver code can now be built using

```bash
# Go to stand-alone driver directory
cd ~/mesh_gluing

# Option 1: Add oomph-lib install directory to PATH and configure
PATH=$PATH:/home/joe_user/oomph_lib_install
cmake -G Ninja -B build

# Option 2: Specify the oomph-lib installation directory during
# the configure step
cmake -G Ninja -B -Doomphlib_ROOT=/home/joe_user/oomph_lib_install

# Option 3: Open the CMakeLists.txt file and add the path to the
# oomph-lib installation to the find_package(...) call by replacing
#   find_package(oomphlib CONFIG REQUIRED)
# with
#   find_package(oomphlib CONFIG REQUIRED PATHS "/home/joe_cool/oomph_lib_install_dir")
# then configure the project
cmake -G Ninja -B build

# Build; ninja without arguments builds both driver codes
cd build/
ninja

# Run 'em!
mkdir RESLT
./mesh_gluing
```

## Creating new demo driver directories

Next we consider how to create a new demo driver directory. This is mainly of interest for developers but, as already mentioned, you can also use this approach to create a directory in the `oomph-lib` directory tree to develop your own code.

As always, it is easiest to copy an existing directory and then modify it. Here we demonstrate how to create a new directory and how to develop the associated self-tests in the `demo_driver/poisson` directory:

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

Next, update the project name to the name of the enclosing folder. (You can, of course, call the project anything you want but this is a sensible default that must be adopted for any demo driver directories to be added to `oomph-lib`.)

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

At the moment, we haven't touched the `validate.sh` script yet, so running `ctest` will fail -- the script is likely to look for the name of the old exectuable, etc. So update the `validate.sh` script and the associated validation data in the `validata` directory. The idea behind the script should (hopefully) be self-explanatory: it runs the code, collates some representative output data into a file which is then compared (allowing for small floating-point errors) against reference data in `validata`. You'll have to decide what to do!

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

At the moment, the newly-created demo driver directory is not visited when the demo-drivers are built (or self-tested) from the top level. So if you build the demo drivers in the `demo_drivers/poisson` directory like this:

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
  my_one_d_poisson # <--- new line
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
    funky_new_equation # <--- new line
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

## Make a new doc directory

`oomph-lib`'s documentation is contained in the directory `doc/`. Its directory structure is roughly equivalent to that of the `demo_drivers/` directory which contains the driver codes that are explained in the documentation.

**WORK IN PROGRESS**

## Build options

You can customise your build by passing flags of the form `-D<FLAG>` to `cmake` during the configuration/generation step. For reference, the table below contains a list of options that the user can control. (Note that the build and installation steps will remain the same.)

Specifying these flags from the command-line can be cumbersome and you may forget which options you used to previously build the project. For this reason, we recommend that you create your own `CMakeUserPresets.json` file, as described in [CMake Presets](#cmake-presets).

**TODO: Discuss desired/not desired options with MH.**

Option                                  | Description                                                                  | Default
----------------------------------------|------------------------------------------------------------------------------|----------
`CMAKE_BUILD_TYPE`                      | The build type (e.g. `Debug`, `Release`, `RelWithDebInfo` or `MinSizeRel`)   | `Release`
`BUILD_SHARED_LIBS`                     | Build using shared libraries; static otherwise  **["SHARED" DOESN'T WORK!]** | OFF
`BUILD_SHARED_LIBS`                     | Build using shared libraries; static otherwise  **["SHARED" DOESN'T WORK!]** | OFF
`OOMPH_DONT_SILENCE_USELESS_WARNINGS`   | Display (harmless) warnings from external_src/ and src/ that are silenced    | OFF
`OOMPH_ENABLE_MPI`                      | Enable the use of MPI for parallel processing                                | OFF
`OOMPH_MPI_NUM_PROC`                    | Number of processes to use with MPI-enabled tests                            | 2
`OOMPH_ENABLE_PARANOID`                 | Enable the PARANOID flag in Debug                                            | OFF
`OOMPH_ENABLE_RANGE_CHECKING`           | Enable RANGE_CHECKING flag in Debug                                          | OFF
`OOMPH_SUPPRESS_TRIANGLE_LIB`           | Suppress build of oomph-lib's copy of the triangle library                   | OFF
`OOMPH_SUPPRESS_TETGEN_LIB`             | Suppress build of oomph-lib's copy of the tetgen library                     | OFF

## CMake Presets

**Work in progress!**

### `CMakePresets.json`

We provide a generic `CMakePresets.json` file in the root directory of the project. To list the available presets, run

```bash
cmake --list-presets
```

We recommend that you can write your own `CMakeUserPresets.json` file. You can inherit your presets from the presets we provide in `CMakePresets.json`. For details on how to do this refer to the [CMake documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

**Remark:** We recommend that you do not use the Ninja Multi-Config generator yet.

**FIXME:** Sort out the clean-up for the multi-config generator. The install manifest doesn't specify the debug config lib files. Hmm...

### `CMakeUserPresets.json` example

```json
{
  "version": 5,
  "configurePresets": [
    {
      "name": "macos_arm64",
      "inherits": "base",
      "displayName": "macos_arm64",
      "generator": "Ninja",
      "binaryDir": "${sourceDir}/build",
      "cacheVariables": {
        "CMAKE_APPLE_SILICON_PROCESSOR": "arm64",
        "CMAKE_BUILD_TYPE": "Release"
      },
      "warnings": {
        "unusedCli": true
      }
    }
  ],
  "buildPresets": [
    {
      "name": "macos_arm64",
      "configurePreset": "macos_arm64",
    }
  ],
  "testPresets": [
    {
      "name": "macos_arm64",
      "inherits": "test-base",
      "configurePreset": "macos_arm64"
    }
  ]
}
```

## Cheat sheet: autotools vs. cmake

Assume that `$oomph_home_dir` is the `oomph-lib` home directory:
<table>
<tr><th>autotools</th><th>cmake</th></tr>

<tr>
<td colspan="2">
<b>Build and install the entire library using a helper script:</b>
</td>
</tr>
<tr>
<td>
<pre>
# Go to home directory
cd $oomph_home_dir<br>
# Run the autogen script which prompts
# for configure options
./autogen.sh

</pre>
</td>
<td>
<pre>
# Go to home directory
cd $oomph_home_dir<br>
# Run the hierher script which prompts
# for configure options
./something_similar.sh
</pre>
</td>
</tr>

<tr>
<td colspan="2">
<b>Alternative: Build and install the entire library step by step:</b>
</td>
</tr>
<tr>
<td>
<pre>
# Go to home directory
cd $oomph_home_dir<br>
# Run manual config/build/install
# sequence
./configure hierher check this
make
make install
</pre>
</td>
<td>
<pre>
# Go to home directory
cd $oomph_home_dir<br>
# Run manual config/build/install
# sequence
cmake -G Ninja -B build
cmake --build build
cmake --install build
</pre>
</td>
</tr>

<tr>
<td colspan="2">
<b>- Run all self-tests:</b>
</td>
</tr>
<tr>
<td>
<pre>
# Go to home directory
cd $oomph_home_dir<br>
# Run the tests
make check
</pre>
</td>
<td>
<pre>
# Go to demo_drivers directory hierher Puneet: does this also do the tests in self_tests?
cd $oomph_home_dir/demo_drivers<br>
# Configure
cmake -G Ninja -B build<br>
# Run the tests (using 8 cores)
cd build
ctest -j 8
</pre>
</td>
</tr>

<tr>
<td colspan="2">
<b>Run all self-tests in a given demo driver directory:</b>
</td>
</tr>
<tr>
<td>
<pre>
# Go to demo driver directory
cd $oomph_home_dir/demo_driver/one_d_poisson<br>
# Run the test
make check
</pre>
</td>
<td>
<pre>
# Go to demo driver directory
cd $oomph_home_dir/demo_driver/one_d_poisson<br>
# Configure
cmake -G Ninja -B build<br>
# Run the test
cd build
ctest
</pre>
</td>
</tr>

<tr>
<td colspan="2">
<b>Debug demo driver code:</b>
</td>
</tr>
<tr>
<td>
<pre>
# Go to demo driver directory
cd $oomph_home_dir/demo_driver/one_d_poisson<br>
# Run the test
make check<br>
# Edit
emacs one_d_poisson.cc<br>
# Run again (make realises that
# code needs to be recompiled)
make check<br>
# Keep going until it works...
# (of course you may do the
# compilation/rerun in your
# editor/IDE; see below)
</pre>
</td>
<td>
<pre>
# Go to demo driver directory
cd $oomph_home_dir/demo_driver/one_d_poisson<br>
# Configure (note that switching on the
# debug build only helps if the library
# itself was compiled with that option too;
#  by default it is built in Release mode
# with full optimisation)
cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=Debug<br>
# Run the test (note that the executable
# got built in the build directory!)
cd build
ctest<br>
# Edit source code
emacs ../one_d_poisson.cc<br>
# Rebuild executable(s)
ninja <br>
# ...or just rebuild a specific one
ninja one_d_poisson<br>
# If you want to see compiler errors
# in detail do
ninja --verbose one_d_poisson<br>
# Now rerun the test
ctest<br>
# Keep going until it works...
# (of course you may do the
# compilation/rerun in your
# editor/IDE; see below)
</pre>
</td>
</tr>

<tr>
<td colspan="2">
<b>Edit a file in the library and rebuild:</b>
</td>
</tr>
<tr>
<td>
<pre>
# Go to relevant src directory
cd $oomph_home_dir/src/poisson<br>
# Edit the relevant file
emacs poisson_elements.h <br>
# Recompile the library
make <br>
# You may go through this repeatedly
# until the code compiles...
<br>
# EMACS/IDE ASIDE: If you
# do this in emacs where you can trigger
# the compilation with [F4] (or some
# other user-defined key) the command to
# be spawned is "make" which emacs enters
# by default (you can edit this); if the
# compilation fails, emacs allows
# you to jump to the source line where
# the compilation error occurs
# using [F12] (on my setup; you can
# customise this), so the typical
# sequence is
# - edit code
# - save
# - [F4] (spawns "make")
# - [RETURN]
# - [F12] (jump to error; edit)
# - [F12] (next error; edit)
# - ...
# - save
# - [F4] (spawns "make")
# - [RETURN]
# etc.
<br>
# Install the library
make install<br>
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
make; make install</b>
</pre>
</td>
<td>
<pre>
# Go to relevant src directory
cd $oomph_home_dir/src/poisson<br>
# Edit the relevant file
emacs poisson_elements.h <br>
# Go back to the home directory and
# rebuild (CMake knows what needs to be
# recompiled)
cd $oomph_home_dir
cmake --build build <br><br>
# You may go through this repeatedly
# until the code compiles...
<br>
# EMACS/IDE ASIDE: If you If you
# do this in emacs where you can trigger
# the compilation with [F4] (or some
# other user-defined key), emacs
# provides the default command "make".
# Edit this to
#
#   cd $oomph_home_dir; cmake --build build
#
# then hit return. If the
# compilation fails, emacs allows
# you to jump to the source line where
# the compilation error occurs
# using [F12] (on my setup; you can
# customise this), so the typical
# sequence is
# - edit code
# - save
# - [F4] (spawns "make")
# - change make to
#
#     cd $oomph_home_dir; cmake --build build
#
# - [RETURN]
# - [F12] (jump to error; edit)
# - [F12] (next error; edit)
# - ...
# - save
# - [F4] (emacs now remembers the previous
#         command)
# - [RETURN]
# etc.
<br>
# Install the library
cd $oomph_home_dir
cmake --install build<br>
</pre>
</td>
</tr>

</table>

## Advanced testing

**TODO:** Add support for `self_test`.

You can filter tests based on the values of `TEST_NAME` in the `oomph_add_test()` test definition. To extract these values, open the `CMakeLists.txt` file in the directory of the test you wish to run. For example, in `demo_drivers/poisson/one_d_poisson/CMakeLists.txt` you will see the following:

```cmake
oomph_add_test(
  TEST_NAME poisson.one_d_poisson
  DEPENDS_ON one_d_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)
```

### Filtering by regex

An alternative approach for filtering tests is to specify a regular expression to the `-R`/`--tests-regex` flag. Only tests for  which the `TEST_NAME` key matches the regular expression will be run. For example

```bash
ctest -R poisson.one_d_poisson
```

will cause all tests containing `poisson.one_d_poisson` in the `TEST_NAME` to be run. To run, say, only `poisson.one_d_poisson` and `gzip.one_d_poisson`, you could use a regex recipe of the form:

```bash
ctest -R '(poisson|gzip)\.one_d_poisson$'
```

### Disabling a test

To temporarily disable a test, you need to set the `DISABLED` property to `TRUE`
using the argument to `TEST_NAME`:

```cmake
# Test definition
oomph_add_test(TEST_NAME poisson.one_d_poisson ...)

# Disable
set_tests_properties(poisson.one_d_poisson PROPERTIES DISABLED YES)
```

## Customising driver codes

You may wish to provide additional information to the build of your executable. A few notable options provided by this function are

- `CXX_OPTIONS`: Compiler flags (e.g. `-Wall`, `-O3`). However, this is likely to only affect your executable and not the library.
- `CXX_DEFINITIONS`: Preprocessor definition(s). Arguments to this keyword do not require a `-D` prefix; CMake will automatically prepend it for you.

For example

```cmake
oomph_add_executable(
  NAME one_d_poisson
  SOURCES one_d_poisson.cc
  LIBRARIES oomph::poisson
  CXX_OPTIONS -Wall -Werror
  CXX_DEFINITIONS REFINEABLE)
```

If you are comfortable with CMake and feel the `oomph_add_executable()` command does not provide the flexibility that you require, you may wish to specify your own executable using the standard CMake functions. If so, you will need to make sure that you add the compile definitions in `OOMPH_COMPILE_DEFINITIONS` to the target, e.g.

```cmake
add_executable(<target-name> <source-1> ... <source-N>)
target_compile_features(<target-name> INTERFACE cxx_std_17)
target_link_libraries(<target-name> PRIVATE oomph::poisson)
target_compile_definitions(<target-name> ${OOMPH_COMPILE_DEFINITIONS})
```

where `<target-name>` is the name of your executable. This imports the compile
definitions defined by `oomph-lib` (during its build) that are needed to make
sure all of the code required is available to your executable.

### Customising targets

For those of you comfortable with CMake, you may wish to control the target properties of executables in a `CMakeLists.txt` file. You may also notice that you are unable to apply target-based CMake commands because CMake is unable to recognise the name of the target you have provided. The reason for this is that inside `oomph_add_executable(...)` we create a unique target name for each executable/test by appending the SHA1 hash of the path to the target. This allows us to provide a unified self-test build (from the base `demo_drivers` folder) that avoid clashes between target names. We do rely on the user never creating two targets with the same name in the same folder but this should always be the case. To use target-based commands on a particular target, create a (SHA1) hash of the path, shorten it to 7 characters, then append it to the original target name and use that name for your commands:

```cmake
# Test definition
oomph_add_executable(NAME <executable-name> ...)

# Construct target name
string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")      # Create hash
string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)            # Shorten to 7 characters
set(OOMPH_TARGET_NAME <executable-name>_${PATH_HASH})  # Append hash

# Do something to the target...
set_target_properties(${OOMPH_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
```

Alternatively, you can use the `oomph_get_target_name(...)` function provided by `oomph-lib`:

```cmake
# Test definition
oomph_add_executable(NAME one_d_poisson ...)

# Get the hashed target name and store it in OOMPH_TARGET_NAME
oomph_get_target_name(one_d_poisson OOMPH_TARGET_NAME)

# Do something to the target...
set_target_properties(${OOMPH_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
```

**Note:** You do not need to append the path hash to test names as, unlike
targets, CMake allows tests to share the same name.

## To be documented

**In progress**:

- [ ] Document CTest usages:
  - [ ] Reading output of failed tests: `cat build/Testing/Temporary/LastTest.log`
  - [ ] List of failed tests: `cat build/Testing/Temporary/LastTestsFailed.log`
  - [ ] Repeat failed test and log output: `--rerun-failed --output-on-failure`
  - [ ] Disable test: `set_tests_properties(<test-name> PROPERTIES DISABLED YES)`
  - [ ] Demand number of processors required by MPI programs:
    - Standard MPI run command: `set_tests_properties(FooWithBar PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})`
    - Variable NP: `set_tests_properties(FooWithBar PROPERTIES PROCESSORS <N>)`, where `N` is the maximum number of processors required

## A deeper dive into the build system

**Work in progress. TEST**

To describe:

- [ ] `OomphLibraryConfig.cmake`
  - [ ] Types of libraries, e.g. regular library, header-only library, etc.
  - [ ] Include paths for each library (`BUILD_INTERFACE`/`INSTALL_INTERFACE`)
  - [ ] Installation
    - [ ] Symlinking headers vs. copying (`OomphCreateSymlinksForHeaders.cmake`)
    - [ ] Combined header (`OomphCreateCombinedHeader.cmake`)
    - [ ] Additional clean-up for symlinks
- [ ] `oomphlibConfig.cmake.in`/`oomphlibConfig.cmake`
  - [ ] Complicated. Save until last...
- [ ] ...
- [ ] `OomphInstallLibrary.cmake`
- [ ] External libraries...

## Helpful CMake resources

For those of you new to CMake, you may wish to consult the following resources:

- The excellently-written "Professional CMake: A Practical Guide" by Scott Craig.
- The [Awesome CMake](https://github.com/onqtam/awesome-cmake) repository.
- [An Introduction to Modern CMake](https://cliutils.gitlab.io/modern-cmake/).
- ...and the list goes on (so add more!).

## Additional instructions

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

#### macOS

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

**Remark:** To make the changes to the `$PATH` variable permanent, add the `export PATH` commands to the end of your shell start-up script, e.g. `.bashrc` or `.zshrc`.

## Community

The original "architects" of `oomph-lib` (in alphabetical order) are
[**Andrew Hazel**](https://github.com/alhazel) (**co-founder**) &
[**Matthias Heil**](https://github.com/MatthiasHeilManchester) (**co-founder**).
Alongside the library founders, `oomph-lib` is currently maintained with the
help of [**Jonathan Deakin**](https://github.com/jondea) and
[**Puneet Matharu**](https://github.com/PuneetMatharu). However, the library has
received (and is still receiving) significant contributions from former/current
project/MSc/PhD students and collaborators. For an exhaustive list, see
[`CONTRIBUTORS.md`](CONTRIBUTORS.md).

If you're interested in joining the team, get in touch. We're always looking for
more help!
