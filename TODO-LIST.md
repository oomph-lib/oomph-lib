# Task list

## Urgent and important
- [ ] Switch to mpic++ and mpicxx compilers for MPI-enabled code
- [ ] Fix broken MPI-enabled demo drivers

## Urgent but not important
- [ ] Empty.

## Not urgent but important
- [ ] Add function `oomphlib(<your-target>)` that implements `target_compile_definitions(<your-target> ${OOMPH_COMPILE_DEFINITIONS})` if the user doesn't use the oomph-lib CMake commands to add drivers, see e.g. `cotire()`.
- [ ] Add code coverage.
- [ ] Add Docker support; see https://github.com/FEniCS/dolfinx
- [ ] Test the build of the demo_drivers with the main build (i.e. add_subdirectory(demo_drivers)).
- [ ] Add sanitiser support (e.g. address/memory, etc.).
- [ ] Sort out external distribution build
- [ ] Add check_...() calls to make sure the C/C++/Fortran compiler work
- [ ] Add (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) in conditionals where appropriate. See: https://cliutils.gitlab.io/modern-cmake/chapters/testing.html
- [ ] Sort out building of oomph-lib with shared libraries. Currently breaks when it gets to oomph_superlu
- [ ] Sort out *sequential* build of MUMPS
- [ ] Sort out a subproject build of private/
- [ ] Find and doc. how users can build a Debug and Release version and easily switch between the two.
- [ ] Remove OomphPrintCompilerDefinitions.cmake.
- [ ] Download CMake oomph-lib in a folder that has a space in the name.

## Not urgent and not important
- [ ] Sort out external distribution build with up-to-date sources.
- [ ] Sort out external sources build with up-to-date sources.
- [ ] Consider replacing bash scripts with Python scripts for demo_driver tests to make them platform-independent.
- [ ] Address MUMPS 64-bit int (i.e. long long) compatibility issue

## Finished
- [x] Remove warnings for external sources on a target basis; see e.g. OomphWarnings.cmake.
- [x] Grab latest version of oomph-lib and add in new code.
- [x] Add Ninja and MPI support to Workflows.
- [x] Test build and self-tests with Linux build system using GitHub Actions and output test results.
- [x] Change all file headers, e.g. in demo drivers ("This file forms part of oomph-lib...") so that they use the correct version and "LastChangedDate". Create a script to automate this! Also update the Copyright notice at the same time.
- [x] Remove set of CMAKE_OSX_DEPLOYMENT_TARGET in demo_drivers/CMakeLists.txt. Find out how else to fix the issue with deployment targets causing a lot of warnings on OSX. [Consensus: I don't care -- it works for me with the SDK on the new MacOS.]
- [x] Ask Matthias about the two_d_unstructured_adaptive_poisson/ validate.sh script.
- [x] Add a change-log. [Already have one, but maybe it's time we actually use it...]
- [x] Set up build of docs. [Handled by Jon's excellent Workflow script.]
- [x] Make sure all demo driver CMake scripts have the "#---- TESTING ---" string
- [x] Sort out .bumpversion.cfg.
- [x] Add CI support.
- [x] Add .clang-format to the build. Need to settle on the formatting choices. Probably worth creating clear examples of each choice.
- [x] Fix issue with running CMake generation step twice when FindMPI is included.
- [x] Create a version.h.in file. [Version info in oomph-lib-config.h]
- [x] Fix support for the oomph-lib-config.h.
- [x] Add PARANOIA and RANGE_CHECKING to possible compile options/decide on how to add it from commandline.
- [x] Retest new version of Scalapack and MUMPS. Patch issue with overridden install paths.
- [x] Sort out construction of oomph-lib-config.h.
- [x] Change conditional statements involving WANT_MPI to OOMPH_HAS_MPI
- [x] Make all "options" in main CMakeLists.txt file OOMPH_ specific so that if it is imported in, there won't be any issues with cache variables.
- [x] Test out demo drivers build with MPI
- [x] Change oomph_parmetis_* libtype to static.
- [x] Replace SuperLU files with an up-to-date version.
- [x] Address simple warning in generic/matrices.cc: "2561:2: warning: 'delete' applied to a pointer that was allocated with 'new[]'; did you mean 'delete[]'? [-Wmismatched-new-delete] delete dist_nrow_local;" (may be more than one case)
- [x] Sort out the mpi build! Add MPI::MPI_CXX to libraries everywhere and set MPI_RUN_COMMAND.
- [x] Add RESOURCE_LOCK properties to MPI-enabled drivers.
- [x] Update SuperLU to v5.2.2 and run demo drivers self-tests (1 fail in shell/plate; need to speak to Matthias).
- [x] Fix demo_drivers/biharmonic alone for now; come back and fix the issue with topologically_rectangular_domain.cc (duplicate symbols). [**A:** Needed to declare functions defined outside a class but inside the header as "inline".]
- [x] Think about how users can specify their own libraries to be built (but not in src/).
- [x] Disable SYMBOLIC_LINKS_FOR_HEADERS if not a Unix-based system.
- [x] For ideas on how to handle ``make check``, see [this link](https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/EmulateMakeCheck)
- [x] Test out whether adding WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES as a compile definition for certain executables actually produce warnings.
- [x] Decide on what to do with the non-mesh src/meshes sources
- [x] Need to make drivers default to the project C++ standard; set this in the config file.
- [x] Need a way to define demo driver requirements, i.e. files that it doesn't link to but it depends on, e.g. a folder of vmtk_files.
- [x] Sort out a subproject build of demo_drivers
- [x] Add a config header and export it with the library!
- [x] Add an "include_guard()" to CMake modules files likely to be included multiple times.
- [x] Specify required dependencies of the package. See: https://cmake.org/cmake/help/git-master/guide/importing-exporting/index.html#id8

