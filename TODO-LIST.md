# Task list

## Requirements for minimum viable product

### Third party library support

* Add support for:
  * [x] `spdlog::spdlog`
  * [x] `nlohmann::json`
  * [x] `Hypre`
  * [x] `MUMPS`
  * [ ] `Trilinos`
    * For Trilinos, link against `Trilinos::all_libs`
  * [ ] Sequential build of `MUMPS`
  * [ ] Specifying the location of already-installed libraries
    * Could possibly specify installation directory with `-D <PACKAGE>_DIR`? E.g. `-D TRILINOS_DIR`, `-D HYPRE_DIR` etc.
    * [ ] `spdlog::spdlog`
    * [ ] `nlohmann::json`
    * [ ] `Hypre`
    * [ ] `MUMPS`
    * [ ] `Trilinos`
      * For Trilinos, link against `Trilinos::all_libs`
    * [ ] Sequential build of `MUMPS`
* [ ] Convert external library version numbers into variables? (ODR)
* [ ] Address MUMPS 64-bit int (i.e. long long) compatibility issue

### Improvements

* [ ] *Optional:* Convert `oomph_gzip` library to single-file?

### Pure-CMake

* [ ] Think about *installing* `fpdiff.py` incase the user wants to install `oomph-lib` but wipe everything but the demo drivers.
  * *If* I decide to do this:
    * [ ] Make a copy of `fpdiff.py` in each folder (seems reasonable as 200 x 15KB = 3MB)
    * [ ] Need to update `OomphAddTest.cmake` CAREFULLY as we need to know where the `fpdiff.py` file will live
      * Maybe create a variable to store the default location (e.g. `PATH_TO_FPDIFF_PY`) that either gets appended to or overwritten when the library and the `fpdiff.py` script is installed
      * **NOTE:** If we append, we may make changes to one script and affect the demo drivers in an unclear way
    * [ ] Update each `validate.sh` script to
      * [ ] Remove the `OOMPH_ROOT_DIR` argument from every `validate.sh`
      * [ ] Remove the `OOMPH_ROOT_DIR` usage from `fpdiff.py` (user can assume that it will be installed in the same folder as the `validate.sh`; i.e. outside the `Validation/` directory)
* [ ] Add two ways to specify arguments to `validate.sh` from `oomph_add_test`
  * One should
* [ ] Add a `oomph_add_pure_cpp_test()` with an `ARGS` command for arguments to pass
  * Should define the executable AND the test target
  * Will allow us to create an individual test per test for finer granularity
* [ ] Add `OOMPH_BLAS_LIB`/`OOMPH_LAPACK_LIB` variables to store the path to the chosen BLAS lib **then** tidy up the Hypre `CMakeLists.txt`
* [ ] *Add presets:*
  * [ ] For MPI configuration
  * [ ] For Intel-based Macs (`--preset macos`) and Arm-based (`--preset macos_arm64`)
    * Should just need to define `CMAKE_APPLE_SILICON_PROCESSOR="arm64"` for the latter
* [ ] Print additional helpful info during configure step, e.g. install destination
* [ ] *Optional:* Update all demo drivers to stop piping output to `validation.log`; use `ctest --o validation.log --output-on-failure`
  * Ideally we'd still output `validation.log` info to this file even if the test doesn't fail... (Just incase there's a bug where a test fails but CTest doesn't catch it)
* [ ] Update `OomphMPI.cmake` to avoid doing a global `link_libraries()`
* [ ] Use `PROJECT_IS_TOP_LEVEL` to enable/disable tests

### Features to add or patch

* [ ] Add Intel MKL BLAS support
* [ ] Replace current BLAS/LAPACK with OpenBLAS
* [ ] Add control for setting number of procs to use in `mpiexec` command
* [ ] Add `self_test/` directory to CMake build
  * [ ] Run tests with C++ fpdiff (requires C++17 support)
  * [ ] Implement with `oomph_add_pure_cpp_test()`
* [ ] Graceful support M1 Macs gracefully and update GitHub Actions...
* [ ] Remove non-working `external_distributions/` crap
* [ ] Fix building of documentation
* [ ] Download CMake oomph-lib in a folder that has a space in the name.
* [ ] Add `check_...()` calls to make sure the C/C++/Fortran/MPI compilers work
* [ ] Fix broken MPI-enabled demo drivers
* [ ] Let the user specify their own BLAS/LAPACK libraries
  * [ ] Create an `OomphPickBLASAndLAPACK.cmake`
  * [ ] Will need to locally set `CMAKE_PREFIX_PATH` and unset `BLA_VENDOR`
* [ ] Test out Ninja multi-config support
  * [ ] If broken, patch any issues, e.g. `ctest --config $<CONFIG>`
* [ ] Sort out a subproject build of `private/`
* [ ] Update GitHub self-tests to only run after pushing .h, .c, .cc, CMakeLists.txt code (see e.g. [here](https://github.com/scivision/mumps/blob/v5.5.1.7/.github/workflows/ci.yml)).
* [ ] Update MPI demo drivers to use 'RESOURCE_GROUPS' to declare the number of CPUs they'll use (will need to store)

### Needs investigation

* [ ] Get timing comparison for Accelerate.Framework BLAS and OpenBLAS on M1 Mac

### Documentation

* [ ] Installing CMake (3.24)
* [ ] Downloading and installing oomph-lib
* [ ] Update instructions for how to run demo drivers
  * [ ] Add info on how to set C++ standard and add compile definitions
  * [ ] macOS info:
    * [ ] Add warning that the user might need to set `CMAKE_APPLE_SILICON_PROCESSOR`
    * [ ] Ask users to `brew uninstall coreutils` if it is not required (over)
    * [ ] Add documentation on how users can also just set the `CMAKE_APPLE_SILICON_PROCESSOR` environment variable (see [CMake GitLab](https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/CMakeDetermineSystem.cmake#L35))
* [ ] Adding demo drivers
  * [ ] Document changes to `validate.sh` (first arg will be `oomph-lib` root; absolute path!!)
  * [ ] ONLY ONE oomph_add_test() per demo_driver directory
* [ ] Adding *MPI* demo drivers (special flags)
  * [ ] Document additional `validate.sh` arg (mpi run command)
* [ ] Document how to use `CMakePresets.json` and `CMakeUserPresets.json` file
* [ ] Find and doc. how users can build both a `Debug` and `Release` version and easily switch between the two
* [ ] Write a breakdown of what has been completed and what features are missing for:
  * [ ] `src/`
    * [ ] `src/meshes/`
      * New build format for `src/meshes`
      * Headers and template headers ONLY; no library artefact! (So no linking to `oomph::meshes` needed)
      * We can automatically find meshes by `#include`-ing `meshes/<mesh-header>`; the `target_include_directories()` command in `OomphLibraryConfig.cmake` has been set up cleverly for that
  * [ ] `external_src/`
    * [ ] New `oomph_gzip` folder
    * [ ] BLAS and LAPACK build can be skipped with `find_package(BLAS)` and `find_package(LAPACK)`
  * [ ] `demo_drivers/`
  * [ ] `private/`
* [ ] Write info on how to debug issue:
  * [ ] Run commands for MPI tests (e.g. `mpiexec -np 2`)
* [ ] Write a breakdown of all new features important changes, e.g.
  * [ ] **Either update Notion or create a CMake Changelog**
  * [ ] C++ implementation of `fpdiff.py`
  * [ ] `validate.sh` scripts now take the path to the root directory
* [ ] Incrementing version number (**strongly recommend using `bumpversion.cfg` to keep git version and cmake version in sync!**)

### Fine-tuning for beta release

Include tasks here that likely need some collaboration with Matthias

* [ ] Talk to MH about updating to C++17
* [ ] Properly review to-do list and see if there's anything missing!
* [ ] Talk to MH about latest MUMPS
  * [ ] Do we need Scotch/METIS?
* [ ] Discuss option naming!! E.g. `WANT_` or `ENABLE_` etc.
* [ ] Possibly change `oomphlib` -> `oomph_lib`? (carefully(!) avoid changing headers/sources)
* [ ] Settle on a versioning scheme (e.g. should we define a `version.h`?)
* [ ] Decide how to handle `OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY` and whether to leave it in
  * Useful for quick debugging of build with demo drivers
  * ~~E.g. if demo drivers are built with the library, do not define the `install()` functions; force the user to work in the build directory!~~
* [ ] Discuss required base presets in `CMakePresets.json`

### Less urgent

## Finished

* [x] ~~Make oomph_add_test() a wrapper around oomph_add_executable()?~~
  * Need to be very careful about how we set up test names and output names to avoid clashes...
  * Will be helpful to create a separate test per executable and so we know whether a specific executable is actually included in a test... (although we'll know automatically from NUM_TESTS)
    * Can't do this at the moment because all tests are run with `validate.sh` so we can only have one test per directory.
    * Can get around this by using the new C++ fpdiff though...
* [x] Implement `fpdiff.py` in C++ (FUN!)
  * [x] Add `zlib` to `external_src` for `.gz` decompression
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
* [x] ~~Switch to mpic++ and mpicxx compilers for MPI-enabled code(?)~~
  * No! We should be letting CMake pick the compilers itself!
