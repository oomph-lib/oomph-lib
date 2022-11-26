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
  <a href="https://codecov.io/gh/puneetmatharu/oomph-lib">
    <img alt="Code coverage" src="https://codecov.io/gh/puneetmatharu/oomph-lib/branch/development/graph/badge.svg?token=ZYGR7Q26I1"/>
  </a>
  <a href="https://www.codacy.com/gh/puneetmatharu/oomph-lib/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=puneetmatharu/oomph-lib&amp;utm_campaign=Badge_Grade">
    <img src="https://app.codacy.com/project/badge/Grade/f64cf70e1f784fc9838f929af4257b41"/>
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
    - [Ubuntu](#ubuntu)
    - [macOS](#macos)
- [Recommended](#recommended)
  - [Ninja](#ninja)
  - [pre-commit](#pre-commit)
- [Usage](#usage)
  - [Building and installing](#building-and-installing)
    - [Doing a fresh rebuild](#doing-a-fresh-rebuild)
    - [Specifying a custom installation location](#specifying-a-custom-installation-location)
    - [Building dependent projects](#building-dependent-projects)
  - [Build options](#build-options)
  - [CMake Presets](#cmake-presets)
    - [`CMakePresets.json`](#cmakepresetsjson)
    - [`CMakeUserPresets.json` example](#cmakeuserpresetsjson-example)
  - [Examples/testing](#examplestesting)
    - [Filtering by label](#filtering-by-label)
    - [Filtering by regex](#filtering-by-regex)
    - [Building a specific demo driver](#building-a-specific-demo-driver)
    - [Clean-up](#clean-up)
    - [Disabling a test](#disabling-a-test)
  - [Development](#development)
    - [Customising targets](#customising-targets)
- [To be documented](#to-be-documented)
- [A deeper dive into the build system](#a-deeper-dive-into-the-build-system)
- [FAQ](#faq)
- [Helpful CMake resources](#helpful-cmake-resources)
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

## Compatibility

Operating system | Support provided?
-----------------|------------------
Ubuntu           | Yes
macOS            | Yes
Windows          | No

## Prerequisites

### CMake

To build this project, you will need CMake 3.22+. You can install `cmake` via your native package manager, e.g. `sudo apt-get install cmake` or `brew install cmake`. If your package manager does not provide a recent enough version of `cmake`, you will need to build it from source. You can find instructions on how to do this below and on the [CMake website](https://cmake.org/install/).

#### Ubuntu

For Linux, CMake provides an installer script to help you download and install CMake:

```bash
# Pick version and installation location
CMAKE_VERSION=3.22.1
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
CMAKE_VERSION=3.22.1
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

## Recommended

### Ninja

We strongly advocate the use of [Ninja](https://github.com/ninja-build/ninja)
for its automatic parallelisation of the build process. Ninja creates clear,
human-readable build files and allows for fast incremental builds.

### pre-commit

We use the [`cmake-format`](https://github.com/cheshirekow/cmake_format) pre-commit hook to automatically format
`CMakeLists.txt` files. To install the pre-commit hooks, see the [Contributing](./CONTRIBUTING.md#pre-commit-hooks-optional) guide.

For this you will need to install `pre-commit`
([available here](https://pre-commit.com/)) using the following

```bash
pip install pre-commit
pre-commit install
```

The `.pre-commit-config.yaml` will take care of the rest. Do not edit the
`.cmake-format.py` file.

## Usage

`oomph-lib` uses CMake to build and install the project.

### Building and installing

To configure, build and install the project using Ninja (recommended), `cd` into the root directory of the cloned `oomph-lib` project and run the following commands:

```bash
cmake -G Ninja -B build   # Configure and generate build system
cmake --build build       # Build
cmake --install build     # Install
```

After the configure step, a `build/` directory will appear with several files in it. During the build step, the individual libraries of `oomph-lib` will be built. Finally, during the install step the headers and generated library files will be installed to the `install/` subdirectory of this project. **It is important to be aware of this design choice as it affects how you use the `oomph-lib` library. For more details, see [Building `oomphlib`-dependent projects](#building-oomphlib-dependent-projects).**

If you would prefer to build `oomph-lib` using the Unix Makefile generator (not recommended), repeat the above steps but omit the `-G Ninja` argument.

To uninstall the project, run the following:

```bash
cd build
ninja uninstall   # replace "ninja" with "make" if you used a Makefile generator
```

If you no longer require any of the build files, you can also delete the `build/` directory.

#### Doing a fresh rebuild

If you wish to do a clean build of the library, you should first uninstall any files that have been installed (refer to the steps above). Once that is done, either delete the `CMakeCache.txt` file inside the `build/` folder or delete the entire `build/` directory itself. Finally, to build `oomph-lib` again, rerun the commands shown at the start of [Building and installing](#building-and-installing).

#### Specifying a custom installation location

By default, `oomph-lib` will be installed to the `install/` subdirectory of the root `oomph-lib` folder. To specify a different installation location, pass `--install-prefix=<install-location>` to `cmake` during the configure step. For example

```bash
cmake -G Ninja -B build --install-prefix=~/oomph_install  # Configure and generate build system
cmake --build build                                       # Build
cmake --install build                                     # Install
```

**Important:** `<install-location>` **must(!)** be an absolute path.

#### Building dependent projects

Recall from [Building and installing](#building-and-installing) that we do not try to install to the user's system paths by default. This choice affects how other projects can consume `oomph-lib`. When trying to configure a CMake project that depends on `oomph-lib`, CMake will check the default system paths for the `oomph-lib` installation. If the library is installed to a non-standard location, i.e. somewhere not in your `PATH` environment variable, CMake will not be able to find it. To get around this issue, you have three choices:

1. Install `oomph-lib` as a superuser:

    ```bash
    cmake -G Ninja -B build -DENABLE_INSTALL_AS_SUPERUSER=ON  # Configure
    cmake --build build                                       # Build
    sudo cmake --install build                                # Install with superuser rights
    ```

    Here, the `-DENABLE_INSTALL_AS_SUPERUSER` flag indicates to CMake that you do not wish it to override the default installation location (to make it install to `install/`). It is only during the install step that you will actually need to use `sudo`.

    When you later wish to uninstall the library, you may need to use `sudo` again, e.g.

    ```bash
    cd build/
    sudo ninja uninstall
    ```

    Now when you try to build a project dependent on `oomph-lib`, CMake will be able to find it without any additional help.

2. Specify the `CMAKE_PREFIX_PATH` variable every time you try to configure a separate subproject:

    ```bash
    # Configure, build and install the library
    cmake -G Ninja -B build --install-prefix=~/oomph_install
    cmake --build build
    cmake --install build

    # Now try to build a project that depends on oomph-lib
    cd ~/some_project_dependent_on_oomphlib/
    cmake -G Ninja -B build -DCMAKE_PREFIX_PATH=~/oomph_install
    ```

3. Add the location of the installation to your `PATH` environment variable:

    ```bash
    # Configure, build and install the library
    cmake -G Ninja -B build --install-prefix=~/oomph_install
    cmake --build build
    cmake --install build

    # Update the PATH environment variable
    # NOTE: This change is only temporary. To make it permanent, add the line below to your .bashrc
    export PATH="$PATH:~/oomph_install"

    # Now try to build a project that depends on oomph-lib
    cd ~/some_project_dependent_on_oomphlib/
    cmake -G Ninja -B build
    ```

    If you use this approach, you should install `oomph-lib` to some location outside of your `oomph-lib` folder; if you use the default installation location and move your `oomph-lib` folder, you will have to update your `PATH` variable to reflect that change.

**Remark:** When you invoke `cmake`, you can specify important variables with flags of the form `-D<variable>=<value>`. These variables are called "cache" variables and take precedence over regular variables and can be used to enable/disable options and other key features during the project configuration.

### Build options

You can customise your build by passing flags of the form `-D<FLAG>` to `cmake` during the configuration/generation step. For reference, the table below contains a list of options that the user can control. (Note that the build and installation steps will remain the same.)

Specifying these flags from the command-line can be cumbersome and you may forget which options you used to previously build the project. For this reason, we recommend that you create your own `CMakeUserPresets.json` file, as described in [CMake Presets](#cmake-presets).

**TODO: Discuss desired/not desired options with MH.**

Option                                    | Description                                                                    | Default
------------------------------------------|--------------------------------------------------------------------------------|--------
`CMAKE_BUILD_TYPE`                        | The build type (e.g. `Debug`, `Release`, `RelWithDebInfo` or `MinSizeRel`)     | `Debug`
`BUILD_SHARED_LIBS`                       | Build using shared libraries; static otherwise  **["SHARED" DOESN'T WORK!]**   | OFF
`BUILD_SHARED_LIBS`                       | Build using shared libraries; static otherwise  **["SHARED" DOESN'T WORK!]**   | OFF
`OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY`   | Build tests with library build                                                 | OFF
`OOMPH_DONT_SILENCE_USELESS_WARNINGS`     | Display (harmless) warnings from external_src/ and src/ that are silenced      | OFF
`OOMPH_ENABLE_MPI`                        | Enable the use of MPI for parallel processing                                  | OFF
`OOMPH_ENABLE_PARANOID`                   | Enable the PARANOID flag in Debug                                              | OFF
`OOMPH_ENABLE_RANGE_CHECKING`             | Enable RANGE_CHECKING flag in Debug                                            | OFF
`OOMPH_ENABLE_SYMBOLIC_LINKS_FOR_HEADERS` | Install symbolic links to the headers instead of copying them                  | ON
`OOMPH_ENABLE_USE_OPENBLAS`               | Use the OpenBLAS implementation of BLAS/LAPACK                                 | OFF
`OOMPH_ENABLE_DOC_BUILD`                  | Suppress Doxygen creation of API documentation **[DOESN'T WORK YET!]**         | OFF
`OOMPH_TRANSITION_TO_VERSION_3`           | Try to build with up-to-date external sources                                  | OFF
`OOMPH_USE_DEPRECATED_SUPERLU`            | Use oomph-lib's deprecated version of SuperLU (4.3)                            | OFF
`OOMPH_USE_DEPRECATED_TETGEN`             | Use oomph-lib's deprecated version of TetGen (1.4)                             | OFF
`OOMPH_SUPPRESS_TRIANGLE_LIB`             | Suppress build of oomph-lib's copy of the triangle library                     | OFF
`OOMPH_SUPPRESS_TETGEN_LIB`               | Suppress build of oomph-lib's copy of the tetgen library                       | OFF
`OOMPH_WANT_NLOHMANN_JSON`                | Enable the [`nlohmann/json`](https://github.com/nlohmann/json) JSON library    | OFF
`OOMPH_WANT_SPDLOG`                       | Enable the [`gabime/spdlog`](https://github.com/gabime/spdlog) logging library | OFF
`OOMPH_WANT_CGAL`                         | Enable we want to build the CGAL library? **[DOESN'T WORK YET!]**              | OFF
`OOMPH_WANT_HYPRE`                        | Enable Hypre library                                                           | OFF
`OOMPH_WANT_MUMPS`                        | Enable MUMPS library [CURRENTLY ONLY WORKING WITH MPI]                         | OFF
`OOMPH_WANT_TRILINOS`                     | Enable Trilinos library  **[DOESN'T WORK YET!]**                               | OFF
`OOMPH_ENABLE_CODE_COVERAGE`              | Enable collection of code coverage results                                     | OFF

### CMake Presets

**Work in progress!**

#### `CMakePresets.json`

We provide a generic `CMakePresets.json` file in the root directory of the project. To list the available presets, run

```bash
cmake --preset list
```

We recommend that you can write your own `CMakeUserPresets.json` file. You can inherit your presets from the presets we provide in `CMakePresets.json`. For details on how to do this refer to the [CMake documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

**Remark:** We recommend that you do not use the Ninja Multi-Config generator yet.

**FIXME:** Sort out the clean-up for the multi-config generator. The install manifest doesn't specify the debug config lib files. Hmm...

#### `CMakeUserPresets.json` example

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
        "CMAKE_BUILD_TYPE": "Release",
        "OOMPH_ENABLE_USE_OPENBLAS": "ON"
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
      "configurePreset": "macos_arm64",
      "output": {
        "labelSummary": true
      }
    }
  ]
}
```

### Examples/testing

**TODO:** Add support for `self_test`.

`oomph-lib` comes with an extensive list of well-documented examples situated in the `demo_drivers/` directory. The driver codes in these folders are also used to validate the library. Before you can run these tests, you must install the `oomph-lib` library using the steps described in [Building and installing](#building-and-installing). To run all of these tests, enter the
`demo_drivers/` folder and run the following:

```bash
cmake -G Ninja -B build   # Configure and generate build system for demo_drivers project
cd build                  # Enter the build/ directory
ctest -j4                 # Execute all tests using 4 processes
```

If you installed `oomph-lib` to a custom location, you will need to tell CMake where it lives during the configuration step. For details on how to do this, refer to [Building dependent projects](#building-dependent-projects).

For the sake of brevity, you can combine all of the above commands into a single command:

```bash
# ctest --build-and-test <source-dir> <build-dir> --build-generator <generator> --test-command <command>
ctest --build-and-test . build --build-generator Ninja --test-command ctest
```

(**Note:** The executable `ctest` is distributed with CMake -- you do not need to install it separately.)

You can filter tests based on the values of `LABELS` or `TEST_NAME` in the `oomph_add_test()` test definition. To extract these values, open the `CMakeLists.txt` file in the directory of the test you wish to run. For example, in `demo_drivers/poisson/one_d_poisson/CMakeLists.txt` you will see the following:

```cmake
oomph_add_test(
  TEST_NAME poisson.one_d_poisson
  TARGET_DEPENDENCIES one_d_poisson
  EXTRA_REQUIRES
  LABELS poisson one_d_poisson)
```

#### Filtering by label

To run the `poisson.one_d_poisson` test based on the `LABELS` key, you can pass the `-L`/`--label-regex` flag to `ctest`, as follows:

```bash
ctest -L poisson -j4         # run all tests with "poisson" in their LABELS
ctest -L one_d_poisson -j4   # run all tests with "one_d_poisson" in their LABELS
```

It is important to note that both of these commands will cause all other tests with similar `LABELS` to be run. This can, however, be particularly helpful when you wish to run a group of tests, e.g. all Poisson-based tests (assuming they have `poisson` under their `LABELS`).

#### Filtering by regex

An alternative approach for filtering tests is to specify a regular expression to the `-R`/`--tests-regex` flag. Only tests for  which the `TEST_NAME` key matches the regular expression will be run. For example

```bash
ctest -R poisson.one_d_poisson
```

will cause all tests containing `poisson.one_d_poisson` in the `TEST_NAME` to be run. To run, say, only `poisson.one_d_poisson` and `gzip.one_d_poisson`, you could use a regex recipe of the form:

```bash
ctest -R '(poisson|gzip)\.one_d_poisson$'
```

#### Building a specific demo driver

A simpler, but slightly more restrictive approach of testing is to only build the demo driver project you want to test. Here, a "project" refers to a folder with a `CMakeLists.txt` file that invokes the `project(...)` command (see e.g. `demo_drivers/poisson/one_d_poisson/CMakeLists.txt`). If you opt for the latter option you will notice that inside each child folder there is a shell script called `validate.sh`, inherited from the old Autotools-based build system, which runs the executables and compares them against the validation data in the `validata` folder.

To test a specific demo driver, you simply step into the folder of the demo driver you want to test and run the commands mentioned in [Examples/testing](#examplestesting). For example, to just run the `one_d_poisson` example, run the commands below:

```bash
cd demo_drivers/poisson/one_d_poisson/  # enter chosen demo driver directory
cmake -G Ninja -B build                 # configure build system for one_d_poisson project
cd build                                # enter build directory
ctest                                   # invoke all tests defined by the one_d_poisson project
```

#### Clean-up

The examples in `demo_drivers/` produce a large amount of data. It is a good idea to remove this data after you run the tests to conserve disk space. To do so, simply delete the `demo_drivers/build/` folder, i.e.

```bash
# Run demo driver tests
cd demo_drivers/
cmake -G Ninja -B build
ctest -j4

# Clean-up
cd ..           # Exit the build/ folder into the parent demo_drivers folder
rm -rf build    # Wipe the self-tests output
```

#### Disabling a test

To temporarily disable a test, you need to set the `DISABLED` property to `TRUE`
using the argument to `TEST_NAME`:

```cmake
# Test definition
oomph_add_test(TEST_NAME poisson.one_d_poisson ...)

# Disable
set_tests_properties(poisson.one_d_poisson PROPERTIES DISABLED YES)
```

### Development

To define your own executable that uses the `oomph-lib` library, you will first
need to import the `oomphlib` package after it has been installed. Once this
has been done, you can define your own executable using the helper function
`oomph_add_executable(...)` ([defined here](cmake/OomphAddExecutable.cmake>)).
For example, to create an executable called `one_d_poisson` from the source
`one_d_poisson.cc` using the Poisson library (`oomph::poisson`), use

```cmake
find_package(oomphlib REQUIRED)
oomph_add_executable(NAME one_d_poisson
                     SOURCES one_d_poisson.cc
                     LIBRARIES oomph::poisson)
```

You may wish to provide additional information to the build of your executable. A few notable options provided by this function are

- `CXX_OPTIONS`: Compiler flags (e.g. `-Wall`, `-O3`). However, this is likely to only affect your executable and not the library.
- `CXX_DEFINITIONS`: Preprocessor definition(s). Arguments to this keyword do not require a `-D` prefix; CMake will automatically prepend it for you.

For example

```cmake
oomph_add_executable(NAME one_d_poisson
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

#### Customising targets

For those of you comfortable with CMake, you may wish to control the target properties of executables in a `CMakeLists.txt` file. You may also notice that you are unable to apply target-based CMake commands because CMake is unable to recognise the name of the target you have provided. The reason for this is that inside `oomph_add_executable(...)` we create a unique target name for each executable/test by appending the SHA1 hash of the path to the target. This allows us to provide a unified self-test build (from the base `demo_drivers` folder) that avoid clashes between target names. We do rely on the user never creating two targets with the same name in the same folder but this should always be the case. To use target-based commands on a particular target, create a (SHA1) hash of the path, shorten it to 7 characters, then append it to the original target name and use that name for your commands:

```cmake
# Test definition
oomph_add_executable(NAME <executable-name> ...)

# Construct target name
string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")      # Create hash
string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)            # Shorten to 7 characters
set(HASHED_TARGET_NAME <executable-name>_${PATH_HASH})  # Append hash

# Do something to the target...
set_target_properties(${HASHED_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
```

Alternatively, you can use the `oomph_get_hashed_target_name(...)` function provided by `oomph-lib`:

```cmake
# Test definition
oomph_add_executable(NAME <executable-name> ...)

# Get the hashed target name
set(HASHED_TARGET_NAME)
oomph_get_hashed_target_name(one_d_poisson HASHED_TARGET_NAME)

# Do something to the target...
set_target_properties(${HASHED_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
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
  - [ ] Demand MPI programs to run in serial: `set_tests_properties(FooWithBar PROPERTIES RUN_SERIAL)`
  - [ ] Note that in all cases here, if you specify the target name, you must rememeber to append the SHA1 path hash.
- [ ] The general workflow:
  - [ ] Creating your own drivers (in e.g. `private/`)

## A deeper dive into the build system

**Work in progress.**

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

## FAQ

<!-- Fill this in... -->

## Helpful CMake resources

For those of you new to CMake, you may wish to consult the following resources:

- The excellently-written "Professional CMake: A Practical Guide" by Scott Craig.
- The [Awesome CMake](https://github.com/onqtam/awesome-cmake) repository.
- [An Introduction to Modern CMake](https://cliutils.gitlab.io/modern-cmake/).
- ...and the list goes on (so add more!).

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
