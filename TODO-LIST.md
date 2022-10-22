# Task list

## Requirements for minimum viable product

* [ ] Discuss option naming!! E.g. `WANT_` or `ENABLE_` etc.
* [ ] Try using Accelerate.Framework BLAS library on M1 Mac
* [ ] Update GitHub self-tests to only run after pushing .h, .c, .cc, CMakeLists.txt code (see e.g. [here](https://github.com/scivision/mumps/blob/v5.5.1.7/.github/workflows/ci.yml)).
* [ ] Convert external library version numbers into variables (for ODR)
* [ ] Make oomph_add_test() a wrapper around oomph_add_executable()?
  * Need to be very careful about how we set up test names and output names to avoid clashes...
  * Will be helpful to create a separate test per executable and so we know whether a specific executable is actually included in a test... (although we'll know automatically from NUM_TESTS)
* [ ] Handle M1 Macs gracefully and update GitHub Actions...
* [ ] Update instructions for how to run demo drivers (need to set CMAKE_APPLE_SILICON_PROCESSOR...?)
* [ ] Implement `fpdiff.py` in C++ (FUN!) (requires gzip processing)
* [ ] Settle on a versioning scheme (e.g. define a `version.h`)
* Documentation for:
  * [ ] Installing CMake (3.24)
  * [ ] Downloading and installing oomph-lib
  * [ ] Incrementing version number (STRONGLY RECOMMEND USING BUMPVERSION.CFG TO KEEP GIT VERSION AND CMAKE VERSION IN LINE!)
  * [ ] Adding demo drivers
    * [ ] Document changes to `validate.sh` (first arg will be `oomph-lib` root; absolute path!!)
    * [ ] ONLY ONE oomph_add_test() per demo_driver directory
  * [ ] Adding *MPI* demo drivers (special flags)
    * [ ] Document additional `validate.sh` arg (mpi run command)
* [ ] Change 'oomphlib' -> 'oomph_lib' (carefully(!) avoid changing headers/sources)
* Add:
  * [x] `spdlog::spdlog`
  * [x] `nlohmann::json`
  * [x] `Hypre`
  * [ ] `MUMPS`
  * [ ] `Trilinos`
    * For Trilinos, link against `Trilinos::all_libs`
  * [ ] Sequential build of `MUMPS`
* [ ] Remove non-working `external_distributions/` crap
* [ ] Patch OomphMPI.cmake to avoid doing a global link_libraries()
* [ ] Write a breakdown of what has been completed and what features are missing for:
  * [ ] `src/`
  * [ ] `external_src/`
  * [ ] `demo_drivers/`
  * [ ] `private/`
* [ ] Talk to MH about updating to latest MUMPS
* [ ] Fix building of documentation
* [ ] Download CMake oomph-lib in a folder that has a space in the name.
* [ ] Find and doc. how users can build a Debug and Release version and easily switch between the two.
* [ ] Sort out a subproject build of private/
* [ ] Add check_...() calls to make sure the C/C++/Fortran/MPI compilers work
* [ ] Add notes on how to use `CMakePresets.json` and how they can define their own `CMakeUserPresets.json` file.
* [ ] Switch to mpic++ and mpicxx compilers for MPI-enabled code(?)
* [ ] Fix broken MPI-enabled demo drivers

## Less urgent

* [ ] Add (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) in conditionals where appropriate. See: <https://cliutils.gitlab.io/modern-cmake/chapters/testing.html>

* [ ] Consider replacing bash scripts with Python scripts for demo_driver tests to make them platform-independent.
* [ ] Address MUMPS 64-bit int (i.e. long long) compatibility issue


## Finished

* [x] Update minimum CMake version
* [x] Fix issue with not being able to rerun `ctest` in a directory twice (at least for `mpi/distribution/adaptive_driven_cavity`).
* [x] Remove Docker build (can be brought back after beta release)
* [x] Remove CUDA crap
* [x] Remove all `if(DEFINED OOMPH_WITH_...)` if-blocks
* [x] Remove OomphPrintCompilerDefinitions.cmake. [EDIT: Changed mind -* it's helpful, so leave it!]
* [x] Test and fix the build of the demo_drivers with the main build (i.e. add_subdirectory(demo_drivers)).
* [x] Add code coverage.
* [x] Remove warnings for external sources on a target basis; see e.g. OomphWarnings.cmake.
* [x] Grab latest version of oomph-lib and add in new code.
* [x] Add Ninja and MPI support to Workflows.
* [x] Test build and self-tests with Linux build system using GitHub Actions and output test results.
* [x] Change all file headers, e.g. in demo drivers ("This file forms part of oomph-lib...") so that they use the correct version and "LastChangedDate". Create a script to automate this! Also update the Copyright notice at the same time.
* [x] Remove set of CMAKE_OSX_DEPLOYMENT_TARGET in demo_drivers/CMakeLists.txt. Find out how else to fix the issue with deployment targets causing a lot of warnings on OSX. [Consensus: I don't care -* it works for me with the SDK on the new MacOS.]
* [x] Ask Matthias about the two_d_unstructured_adaptive_poisson/ validate.sh script.
* [x] Add a change-log. [Already have one, but maybe it's time we actually use it...]
* [x] Set up build of docs. [Handled by Jon's excellent Workflow script.]
* [x] Make sure all demo driver CMake scripts have the "#---* TESTING ---" string
* [x] Sort out .bumpversion.cfg.
* [x] Add CI support.
* [x] Add .clang-format to the build. Need to settle on the formatting choices. Probably worth creating clear examples of each choice.
* [x] Fix issue with running CMake generation step twice when FindMPI is included.
* [x] Create a version.h.in file. [Version info in oomph-lib-config.h]
* [x] Fix support for the oomph-lib-config.h.
* [x] Add PARANOIA and RANGE_CHECKING to possible compile options/decide on how to add it from commandline.
* [x] Retest new version of Scalapack and MUMPS. Patch issue with overridden install paths.
* [x] Sort out construction of oomph-lib-config.h.
* [x] Change conditional statements involving WANT_MPI to OOMPH_HAS_MPI
* [x] Make all "options" in main CMakeLists.txt file OOMPH_ specific so that if it is imported in, there won't be any issues with cache variables.
* [x] Test out demo drivers build with MPI
* [x] Change oomph_parmetis_* libtype to static.
* [x] Replace SuperLU files with an up-to-date version.
* [x] Address simple warning in generic/matrices.cc: "2561:2: warning: 'delete' applied to a pointer that was allocated with 'new[]'; did you mean 'delete[]'? [-Wmismatched-new-delete] delete dist_nrow_local;" (may be more than one case)
* [x] Sort out the mpi build! Add MPI::MPI_CXX to libraries everywhere and set MPI_RUN_COMMAND.
* [x] Add RESOURCE_LOCK properties to MPI-enabled drivers.
* [x] Update SuperLU to v5.2.2 and run demo drivers self-tests (1 fail in shell/plate; need to speak to Matthias).
* [x] Fix demo_drivers/biharmonic alone for now; come back and fix the issue with topologically_rectangular_domain.cc (duplicate symbols). [**A:** Needed to declare functions defined outside a class but inside the header as "inline".]
* [x] Think about how users can specify their own libraries to be built (but not in src/).
* [x] Disable SYMBOLIC_LINKS_FOR_HEADERS if not a Unix-based system.
* [x] For ideas on how to handle ``make check``, see [this link](https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/EmulateMakeCheck)
* [x] Test out whether adding WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES as a compile definition for certain executables actually produce warnings.
* [x] Decide on what to do with the non-mesh src/meshes sources
* [x] Need to make drivers default to the project C++ standard; set this in the config file.
* [x] Need a way to define demo driver requirements, i.e. files that it doesn't link to but it depends on, e.g. a folder of vmtk_files.
* [x] Sort out a subproject build of demo_drivers
* [x] Add a config header and export it with the library!
* [x] Add an "include_guard()" to CMake modules files likely to be included multiple times.
* [x] Specify required dependencies of the package. See: <https://cmake.org/cmake/help/git-master/guide/importing-exporting/index.html#id8>
