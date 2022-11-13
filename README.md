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
        <th><a href="../../tree/feature-convert-to-cmake-build-system"><code>convert-to-cmake-build-system</code></a></th>
    </tr>
    <tr>
        <td>Ubuntu</td>
        <td>
            <a href="../../actions/workflows/self-tests-ubuntu.yaml">
                <img alt="Ubuntu" src="../../actions/workflows/self-tests-ubuntu.yaml/badge.svg?branch=main">
            </a>
        </td>
        <td>
            <a href="../../actions/workflows/self-tests-ubuntu.yaml">
                <img alt="Ubuntu" src="../../actions/workflows/self-tests-ubuntu.yaml/badge.svg?branch=development">
            </a>
        </td>
        <td>
          <div align="center">
              <a href="../../actions/workflows/test-ubuntu.yaml">
                  <img alt="Ubuntu" src="../../actions/workflows/test-ubuntu.yaml/badge.svg?branch=feature-convert-to-cmake-build-system">
              </a>
          </div>
        </td>
    </tr>
    <tr>
        <td>macOS</td>
        <td>
            <a href="../../actions/workflows/self-tests-macos.yaml">
                <img alt="macOS" src="../../actions/workflows/self-tests-macos.yaml/badge.svg?branch=main">
            </a>
        </td>
        <td>
            <a href="../../actions/workflows/self-tests-macos.yaml">
                <img alt="macOS" src="../../actions/workflows/self-tests-macos.yaml/badge.svg?branch=development">
            </a>
        </td>
        <td>
          <div align="center">
              <a href="../../actions/workflows/test-macos.yaml">
                  <img alt="macOS" src="../../actions/workflows/test-macos.yaml/badge.svg?branch=feature-convert-to-cmake-build-system">
              </a>
          </div>
        </td>
    </tr>
    </table>
</div>

## Table of contents

- [Table of contents](#table-of-contents)
- [Description](#description)
- [Compatibility](#compatibility)
- [Documentation to-do list](#documentation-to-do-list)
- [Prerequisites](#prerequisites)
  - [CMake](#cmake)
- [Recommended](#recommended)
  - [Ninja](#ninja)
  - [pre-commit](#pre-commit)
- [Usage](#usage)
  - [Building and installing](#building-and-installing)
    - [Specifying a custom installation location](#specifying-a-custom-installation-location)
  - [Uninstalling](#uninstalling)
  - [Build options](#build-options)
  - [CMake Presets](#cmake-presets)
    - [`CMakeUserPresets.json` example](#cmakeuserpresetsjson-example)
  - [Examples/testing](#examplestesting)
    - [Filtering by label](#filtering-by-label)
    - [Filtering by regex](#filtering-by-regex)
    - [Clean-up](#clean-up)
  - [Disabling a test](#disabling-a-test)
  - [Development](#development)
- [The build system in detail](#the-build-system-in-detail)
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

## Documentation to-do list

Finish documenting the following:

- [x] Basic build instructions.
- [x] How to build with Ninja.
- [ ] Adding a new library.
- [ ] Adding a new demo-driver.
- [ ] Packaging with CPack.

## Prerequisites

### CMake

To build this project, you will need CMake 3.24+. You can install `cmake` via your native package manager, e.g. `sudo apt-get install cmake` or `brew install cmake`. If your package manager does not provide a recent enough version of `cmake`, you will need to build it from source. You can find instructions on how to do this [here](https://cmake.org/install/).

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
>>> cmake -G Ninja -B build   # Configure and generate build system
>>> cmake --build build       # Build
>>> cmake --install build     # Install
```

After the configure step, a `build/` directory will appear with several files in it. During the build step, the individual libraries of `oomph-lib` will be built. Finally, the install step will cause the headers and generated library files to be installed to the user's system paths.

If you would prefer to use the Unix Makefile generator (not recommended), repeat the above steps but omit the `-G Ninja` argument.

**Tip:** If you wish to do a clean build of the library, you can either delete the `build/` directory or run the initial `cmake` command with the `--fresh` flag.

**TODO:** *Put a note here on out-of-source builds*.

#### Specifying a custom installation location

By default, `oomph-lib` will be installed to `/usr/local/` on Unix systems. To specify a custom installation location, pass `--install-prefix=<install-location>` to `cmake` during the configure step. For example

```bash
>>> cmake -G Ninja -B build --install-prefix=~/oomph_install   # Configure and generate build system
>>> cmake --build build                                       # Build
>>> cmake --install build                                     # Install
```

**Note:** `<install-location>` **must(!)** be an absolute path.

**Recommendation:** Do not specify a non-standard installation location if you have superuser rights on your machine. When you try to build a project that calls `find_package(oomphlib ...)` (e.g. `demo_drivers`), CMake will check the default system paths for the `oomph-lib` installation. If the library is installed to a non-standard location, CMake will not be able to find it. As a result, you will either need to pass the location of the installation to `cmake` using one of the following options:

1. Specify the `CMAKE_PREFIX_PATH` variable (every time you try to configure a separate subproject!), or
2. You will need to add the location of the installation to your environment `PATH` variable:

```bash
# Build and install the main library
>>> cmake --install-prefix=~/oomph_install -G Ninja -B build
>>> cmake --build build
>>> cmake --install build

# Building, e.g., the demo drivers
cd demo_drivers/

# Option 1: Specify the installation location
>>> cmake -G Ninja -B build -DCMAKE_PREFIX_PATH=~/oomph_install

# Option 2: Update the PATH environment variable
>>> export PATH="$PATH:~/oomph_install"
```

**Remark:** When you invoke `cmake`, you can specify important variables with flags of the form `-D<variable-name>` or `-D <variable-name>`. These variables are called "cache" variables and take precedence over regular variables and can be use to enable/disable options during the project configuration.

### Uninstalling

To uninstall the project, enter the `build` folder, uninstall the installed project files and delete the build folder using the following

```bash
>>> cd build
>>> ninja uninstall   # replace ninja with make if using Makefile generator
>>> cd ..
>>> rm -rf build
```

### Build options

To customise your build, provide arguments of the form `-D<YOUR-FLAG-HERE>` at the CMake configuration/generation step. The table below contains a list of options that the user can control.

Specifying these flags from the command-line can be cumbersome and you may forget what options you used to previously build the project. For this reason, we recommend that you create your own `CMakeUserPresets.json` file, as described in [CMake Presets](#cmake-presets).

**TODO: Discuss with MH which options to make available to the user.**

Option                                    | Description                                                                    | Default
------------------------------------------|--------------------------------------------------------------------------------|--------
`CMAKE_BUILD_TYPE`                        | The build type (e.g. `Debug`, `Release`, `RelWithDebInfo` or `MinSizeRel`)     | `Debug`
`BUILD_SHARED_LIBS`                       | Build using shared libraries; static otherwise                                 | OFF
`OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY`   | Build tests with library build                                                 | OFF
`OOMPH_DONT_SILENCE_USELESS_WARNINGS`     | Display (harmless) warnings from external_src/ and src/ that are silenced      | OFF
`OOMPH_ENABLE_MPI`                        | Enable the use of MPI for parallel processing                                  | OFF
`OOMPH_ENABLE_PARANOID`                   | Enable the PARANOID flag in Debug                                              | OFF
`OOMPH_ENABLE_RANGE_CHECKING`             | Enable RANGE_CHECKING flag in Debug                                            | OFF
`OOMPH_ENABLE_SYMBOLIC_LINKS_FOR_HEADERS` | Replace headers by symbolic links                                              | ON
`OOMPH_ENABLE_USE_OPENBLAS`               | Use the OpenBLAS implementation of BLAS/LAPACK                                 | OFF
`OOMPH_ENABLE_DOC_BUILD`                  | Suppress Doxygen creation of API documentation **[DOESN'T WORK YET!]**         | OFF
`OOMPH_TRANSITION_TO_VERSION_3`           | Try to build with up-to-date external sources                                  | OFF
`OOMPH_USE_DEPRECATED_SUPERLU`            | Use oomph-lib's deprecated version of SuperLU (4.3)                            | OFF
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

We provide a generic `CMakePresets.json` file in the root directory of the project. We recommend that you can write your own `CMakeUserPresets.json` file.

For details on CMake presets refer to the [CMake documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

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

`oomph-lib` comes with an extensive list of well-documented examples situated in the `demo_drivers/` directory. The driver codes in these folders are also used to validate the library. Before you can run these tests, you must install the `oomph-lib` library using the steps described in [Building and installing](#building-and-installing). To run all of these tests, enter the
`demo_drivers/` folder and run the following:

```bash
>>> cd demo_drivers/
>>> cmake -G Ninja -B build   # Configure and generate build system for demo_drivers project
>>> ctest -j4                 # Enter the build folder and execute all tests with 4 jobs
```

(The executable `ctest` is distributed with CMake -- you do not need to install it separately.)

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
>>> ctest -L poisson -j4         # run all tests with "poisson" in their LABELS
>>> ctest -L one_d_poisson -j4   # run all tests with "one_d_poisson" in their LABELS
```

It is important to note that both of these commands will cause all other tests with similar `LABELS` to be run. This can, however, be particularly helpful when you wish to run a group of tests, e.g. all Poisson-based tests (assuming they have `poisson` under their `LABELS`).

#### Filtering by regex

An alternative approach for filtering tests is to specify a regular expression to the `-R`/`--tests-regex` flag. Only tests for  which the `TEST_NAME` key matches the regular expression will be run. For example

```bash
>>> ctest -R poisson.one_d_poisson
```

will cause all tests containing `poisson.one_d_poisson` in the `TEST_NAME` to be run. To run, say, only `poisson.one_d_poisson` and `gzip.one_d_poisson`, you could use a regex recipe of the form:

```bash
ctest -R '(poisson|gzip)\.one_d_poisson$'
```

**TODO:** Patch support for `self_test`.

#### Clean-up

The examples in `demo_drivers/` produce a large amount of data. It is a good idea to remove this data after you run the tests to conserve disk space. To
do so, simply delete the `demo_drivers/build/` folder, i.e.

```bash
# Run demo driver tests
>>> cd demo_drivers/
>>> cmake -G Ninja -B build
>>> ctest -j4

# Clean-up
>>> cd ..           # Exit the build/ folder into the parent demo_drivers folder
>>> rm -rf build    # Wipe the self-tests output
```

The approach described above allows you to test the entire `oomph-lib` build,
all at once. However you may wish to test a smaller subset of these problems or
just one. To do this, you may either:

- (i) provide `ctest` with a filter to select the tests that you wish to run (described further below), or
- (ii) enter any child project and rerun the same commands as above.

Here, a child project refers to any subfolder containing a `CMakeLists.txt`
file that invokes the `project(...)` command (e.g.
`demo_drivers/poisson/one_d_poisson`). If you opt for the latter option you
will notice that inside each child folder there is a shell script called
`validate.sh`, inherited from the old Autotools-based build system, which runs
the executables and compares them against the validation data in the
`validata` folder. **You should not edit the `validate.sh` scipt or the data
in `validata`.**

For those of you comfortable with CMake, you may wish to control the target
properties of executables in a `CMakeLists.txt` file. You may also notice that
you are unable to apply target-based CMake commands because CMake is unable to
recognise the name of the target you have provided. The reason for this is that
inside `oomph_add_executable(...)` we create a unique target name for each
executable/test by appending the SHA1 hash of the path to the target. This allows
us to provide a unified self-test build (from the base `demo_drivers` folder)
that avoid clashes between target names. We do rely on the user never creating
two targets with the same name in the same folder but this should always be the
case. To use target-based commands on a particular target, create a (SHA1) hash
of the path, shorten it to 7 characters, then append it to the original target
name and use that name for your commands:

```cmake
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")           # Create hash
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)                 # Shorten to 7 characters
  set(HASHED_TARGET_NAME <YOUR-EXECUTABLE-NAME>_${PATH_HASH})  # Append hash
```

**Note:** You do not need to append the path hash to test names as, unlike
targets, CMake allows tests to share the same name.

### Disabling a test

To temporarily disable a test, you need to set the `DISABLED` property to `TRUE`
using the argument to `TEST_NAME`:

```cmake
# Test definition
oomph_add_test(TEST_NAME poisson.one_d_poisson ...)

# Disable
set_tests_properties(poisson.one_d_poisson PROPERTIES DISABLED YES)
```

**In progress**:

- [ ] Add a "make self-test" command for the root `oomph-lib` directory which executes all of the self-tests.
- [ ] Document CTest usages:
  - [ ] Parallel execution; append a `-j <N>` flag.
  - [ ] Test filtering:
    - [ ] Filter by labels `-L <label>` or
    - [ ] Regular expression matching; run all tests beginning with poisson: `-R '^poisson'`
    - [ ] Run all but `<test-name>`: `-E <test-name>`
    - [ ] Run ony `<test-1>` and `<test-2>`: `-R '<test-1>|<test-2>'`
  - [ ] Reading output of failed tests: `cat build/Testing/Temporary/LastTest.log`
  - [ ] List of failed tests: `cat build/Testing/Temporary/LastTestsFailed.log`
  - [ ] Repeating failed tests: `--rerun-failed`
  - [ ] Repeat failed test and log output: `--rerun-failed --output-on-failure`
  - [ ] Disable test: `set_tests_properties(<test-name> PROPERTIES DISABLED YES)`
  - [ ] Demand parallel codes to run in serial: `set_tests_properties(FooWithBar PROPERTIES RUN_SERIAL)`;
  - [ ] Providing a RESOURCE_LOCK for parallel codes.
  - [ ] Note that in all cases here, if you specify the target name, you must rememeber to append the SHA1 path hash.

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

You may wish to provide additional information to the build of your executable.
A few notable options provided by this function are

- `CXX_STANDARD`: The C++ standard. The only arguments we currently allow are 11, 14, or 17 (corresponding to C++11, C++14, and C++17, respectively). We currently assume C++11 for all files in the library. Specifying a more modern standard may result in unexpected consequences. Don't say we didn't warn you!
- `CXX_OPTIONS`: Compiler flags (e.g. `-Wall`, `-O3`). However, this is likely to only affect your executable and not the library. (`TODO: Find out about this!`)
- `CXX_DEFINITIONS`: Preprocessor definition(s). Arguments to this keyword do not require a `-D` prefix; CMake will automatically prepend it for you.

For example

```cmake
  oomph_add_executable(NAME one_d_poisson
                       SOURCES one_d_poisson.cc
                       LIBRARIES oomph::poisson
                       CXX_STANDARD 11
                       CXX_OPTIONS -Wall -Werror
                       CXX_DEFINITIONS REFINEABLE)
```

If you are comfortable with CMake and you wish to specify your own executable
using the standard CMake functions then make sure to add the following line
after calling `add_executable`

```cmake
  target_compile_definitions(<your-target> ${OOMPH_COMPILE_DEFINITIONS})
```

where `<your-target>` is the name of your executable. This imports the compile
definitions defined by `oomph-lib` (during its build) that are needed to make
sure all of the code required is available to your executable.

## The build system in detail

**Work in progress.**

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
