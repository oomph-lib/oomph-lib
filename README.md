[![oomph-lib](https://github.com/PuneetMatharu/oomph-lib-cmake/blob/dev/doc/figures/oomph_logo.png)](https://github.com/PuneetMatharu/oomph-lib-cmake/releases)

- [Description](#description)
- [Compatibility](#compatibility)
- [Documentation to-do list](#documentation-to-do-list)
- [Recommended](#recommended)
- [Usage](#usage)
  - [Building and installing](#building-and-installing)
  - [Useful CMake flags](#useful-cmake-flags)
  - [Examples/testing](#examplestesting)
  - [Uninstall](#uninstall)
  - [Development](#development)
  - [Packaging](#packaging)
- [Helpful CMake resources](#helpful-cmake-resources)
- [Authors](#authors)

## Description

The [`oomph-lib` homepage](http://www.oomph-lib.org) provides much more detail
on tutorials, coding conventions, licencing information, etc. Provided you have
downloaded a distribution that contains the documentation and you have the
required tools (mainly ``doxygen``; get it [here](http://www.doxygen.org>)
available on your machine, the installation procedure will create a local copy
of the oomph-lib webpages and the entire online documentation in the ``doc``
directory. In particular, ``doc/html/index.html`` is a local copy of the
oomph-lib homepage.

## Compatibility

The CMake-based build of the oomph-lib library has been tested on MacOS and
Linux thus far; it has not been tested on Windows yet.

## Documentation to-do list

Finish documenting the following:

- [x] Basic build instructions.
- [x] How to build with Ninja.
- [ ] Usage via FetchContent.
- [ ] Adding a new library.
- [ ] Adding a new demo-driver.
- [ ] Packaging with CPack.

## Recommended

We strongly advocate the use of [Ninja](https://github.com/ninja-build/ninja)
for its automatic parallelisation of the build process. Ninja creates clear,
human-readable build files and allows for fast incremental builds.

We use the ``cmake-format`` pre-commit hook to automatically format
``CMakeLists.txt`` files. For this you will need to install ``pre-commit``
([available here](https://pre-commit.com/)) using the following
```bash
  pip install pre-commit
```
The ``.pre-commit-config.yaml`` will take care of the rest. Do not edit the
``.cmake-format.json`` file.

## Usage

### Building and installing

To configure, build and install the project using Ninja (recommended), ``cd``
into the root directory of the oomph-lib project and run the following:
```bash
  cmake -G Ninja -B build              # Configure and generate build system
  cd build && ninja && ninja install   # Build and install
```
If you'd prefer to use Makefile Generators for your build system instead, run
```bash
  cmake -B build
  cd build && make && make install
```
### Useful CMake flags

To customise your build, provide arguments of the form ``-D<YOUR-FLAG-HERE>``
at the CMake configuration/generation step.

To specify a ``Release`` build (i.e. optimised; the default is ``Debug``) use
```bash
  -DCMAKE_BUILD_TYPE=Release
```
To specify the installation location append the flag
```bash
  -DCMAKE_INSTALL_PREFIX=~/my-custom-install
```
To enable the use of MPI (if available on your system) use
```bash
  -DOOMPH_ENABLE_MPI
```
### Examples/testing

``oomph-lib`` comes with an extensive list of well-documented examples situated
in the ``demo_drivers`` directory. The driver codes in these folders are also
used to validate the library. To run all of these tests, enter the
``demo_drivers`` folder and run the following:
```bash
  cmake -G Ninja -B build   # Configure and generate build system for demo_drivers project
  cd build && ctest         # Enter the build folder and execute all tests
```
If you intend to use Makefile Generators, remove "``-G Ninja``" from the first
command. Note that unlike the lightweight unit-tests in the ``tests/`` folder,
these "integration tests" are more intensive and take much longer to complete.
After running all of the self-tests, you may wish to get rid of the output. To
do so, simply delete the ``build`` folder, i.e.
```bash
  cd ..           # Exit the build folder into the parent demo_drivers folder
  rm -rf build    # Wipe the self-tests output
```
The approach described above allows you to test the entire ``oomph-lib`` build,
all at once. However you may wish to test a smaller subset of these problems or
just one. To do this, you may either:

- (i) provide ``ctest`` with a filter to select the tests that you wish to run (described further below), or
- (ii) enter any child project and rerun the same commands as above.

Here, a child project refers to any subfolder containing a ``CMakeLists.txt``
file that invokes the ``project(...)`` command (e.g.
``demo_drivers/poisson/one_d_poisson``). If you opt for the latter option you
will notice that inside each child folder there is a shell script called
``validate.sh``, inherited from the old Autotools-based build system, which runs
the executables and compares them against the validation data in the
``validata`` folder. **You should not edit the ``validate.sh`` scipt or the data
in ``validata``.**

For those of you comfortable with CMake, you may wish to control the target
properties of executables in a ``CMakeLists.txt`` file. You may also notice that
you are unable to apply target-based CMake commands because CMake is unable to
recognise the name of the target you have provided. The reason for this is that
inside ``oomph_add_executable(...)`` and ``oomph_add_test(...)`` we create a
unique target name for each executable/test by appending the SHA1 hash of the
path to the target. This allows us to provide a unified self-test build (from
the base ``demo_drivers`` folder) that avoid clashes between target names. We do
rely on the user never creating two targets with the same name in the same
folder but this should always be the case. To use target-based commands on a
particular target, create a hash of the path, append it to the original target
name and use that name for your commands:
```cmake
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")           # Create hash
  set(HASHED_TARGET_NAME <YOUR-EXECUTABLE-NAME>_${PATH_HASH})  # Append hash
```
**In progress**:

- [ ] Add a "make self-test" command for the root oomph-lib directory which executes all of the self-tests.
- [ ] Document CTest usages:
   - [ ] Parallel execution; append a ``-j <N>`` flag.
   - [ ] Test filtering:
      - [ ] Filter by labels ``-L <label>`` or
      - [ ] Regular expression matching; run all tests beginning with poisson: ``-R '^poisson'``
      - [ ] Run all but ``<test-name>``: ``-E <test-name>``
      - [ ] Run ony ``<test-1>`` and ``<test-2>``: ``-R '<test-1>|<test-2>'``
   - [ ] Reading output of failed tests: ``cat build/Testing/Temporary/LastTest.log``
   - [ ] List of failed tests: ``cat build/Testing/Temporary/LastTestsFailed.log``
   - [ ] Repeating failed tests: ``--rerun-failed``
   - [ ] Repeat failed test and log output: ``--rerun-failed --output-on-failure``
   - [ ] Disable test: ``set_tests_properties(<test-name> PROPERTIES DISABLED YES)``
   - [ ] Demand parallel codes to run in serial: ``set_tests_properties(FooWithBar PROPERTIES RUN_SERIAL)``;
   - [ ] Providing a RESOURCE_LOCK for parallel codes.
   - [ ] Note that in all cases here, if you specify the target name, you must rememeber to append the SHA1 path hash.

### Uninstall

To uninstall the project, enter the ``build`` folder, remove the installed
project files and delete the build folder using the following
```bash
  cd build
  ninja uninstall   # replace ninja with make if using Autotools
  cd ..
  rm -rf build
```
### Development

To define your own executable that uses the oomph-lib library, you will first
need to import the ``oomphlib`` package after it has been installed. Once this
has been done, you can define your own executable using the helper function
``oomph_add_executable(...)`` ([defined here](cmake/OomphAddExecutable.cmake>)).
For example, to create an executable called ``one_d_poisson`` from the source
``one_d_poisson.cc`` using the Poisson library (``oomph::poisson``), use
```cmake
  find_package(oomphlib REQUIRED)
  oomph_add_executable(NAME one_d_poisson
                       SOURCES one_d_poisson.cc
                       LIBRARIES oomph::poisson)
```
You may wish to provide additional information to the build of your executable.
A few notable options provided by this function are

- ``CXX_STANDARD``: The C++ standard. The only arguments we currently allow are 11, 14, or 17 (corresponding to C++11, C++14, and C++17, respectively). We currently assume C++11 for all files in the library. Specifying a more modern standard may result in unexpected consequences. Don't say we didn't warn you!
- ``CXX_OPTIONS``: Compiler flags (e.g. ``-Wall``, ``-O3``). However, this is likely to only affect your executable and not the library. (``TODO: Find out about this!``)
- ``CXX_DEFINITIONS``: Preprocessor definition(s). Arguments to this keyword do not require a ``-D`` prefix; CMake will automatically prepend it for you.

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
after calling ``add_executable``
```cmake
      target_compile_definitions(<your-target> ${OOMPH_COMPILE_DEFINITIONS})
```
where ``<your-target>`` is the name of your executable. This imports the compile
definitions defined by ``oomph-lib`` (during its build) that are needed to make
sure all of the code required is available to your executable.

### Packaging

**Work in progress.**

## Helpful CMake resources

For those of you new to CMake, you may wish to consult the following resources:

* The excellently-written "Professional CMake: A Practical Guide" by Scott Craig.
* The [Awesome CMake](https://github.com/onqtam/awesome-cmake) repository.
* [An Introduction to Modern CMake](https://cliutils.gitlab.io/modern-cmake/).
* ...and the list goes on (so add more!).

## Authors

`oomph-lib` is developed and maintained by Matthias Heil and Andrew Hazel at the
School of Mathematics, University of Manchester, UK, assisted by many students,
postdocs and collaborators. The authors may be contacted by email

> oomph-lib[at]maths.man.ac.uk

or via snail mail

> School of Mathematics
> University of Manchester
> Oxford Road
> Manchester M13 9PL
> UK

Constructive feedback on any aspect of the library is most welcome!
