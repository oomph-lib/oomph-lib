# Task list

## Fine-tuning for beta release

Include tasks here that likely need some collaboration with Matthias

* [ ] Update doc/FAQ/...
  * [ ] E.g. "The oomph-lib distribution includes some third-party libraries. How do I get the code to link against optimised local versions of these libraries that are already installed on my machine?"
* [ ] Kill all scripts that will no longer be required, e.g. `move_external_libraries_and_distributions_to_permanent_location.bash`
* [ ] Discuss switching to OpenBLAS only (need to change OOMPH_USE_BLAS/LAPACK_FROM)

* [ ] Go through all new commits to the main repository in the last 1-2 years and make sure any changes in Makefiles have been reflected in their corresponding CMakeLists.txt file
* [ ] Possibly change `oomphlib` -> `oomph_lib`? (carefully(!) avoid changing headers/sources)
  * **CHANGE IT**

* [ ] Test with `OOMPH_TRANSITION_TO_VERSION_3`
* [ ] Add demo driver test target to "all" so `ninja` in the build directory will also cause required test files to be copied over
  * **YES DO THIS!**

### Improvements

* [ ] Discuss the Doxygen search engine; `SEARCHENGINE=NO` at the moment...
* [ ] *Optional:* Convert `oomph_gzip` library to single source file?

### GitHub-y

* [ ] Add a `destruct_test.yaml` to test lots of different configurations (e.g. with and without MPI, w/ and w/o Hypre, etc.)

### Pure-CMake

* [ ] Propagate compiler flags to users executables
* [ ] Patch the use of shared libs
* [ ] Update MPI tests to specify [PROCESSORS](https://cmake.org/cmake/help/latest/prop_test/PROCESSORS.html#prop_test:PROCESSORS) rather than forcing everything to run sequentially
* [x] Update `oomph_add_test()` to take args for `validate.sh` and explicitly pass e.g. "${OOMPH_ROOT_DIR}" and "${OOMPH_MPI_RUN_COMMAND}"
  * Should help make things transparent and easier for users to control
* [ ] Update every test name to be the path to the demo driver; e.g. `mpi.distribution.adaptive_driven_cavity`
* [ ] Patch `ninja oomph_uninstall` command for Ninja Multi-Config
* [ ] Update build system to make sure building a demo driver target without `ctest` also runs the `check_...` target to build the target and copy files over
* [ ] *Add presets:*
  * [ ] For MPI configuration
  * [ ] For Intel-based Macs (`--preset macos`) and Arm-based (`--preset macos_arm64`)
    * Should just need to define `CMAKE_APPLE_SILICON_PROCESSOR="arm64"` for the latter
* [ ] *Optional:* Update all demo drivers to stop piping output to `validation.log`; use `ctest --o validation.log --output-on-failure`
  * Ideally we'd still output `validation.log` info to this file even if the test doesn't fail... (Just incase there's a bug where a test fails but CTest doesn't catch it)

### Features to check, add or patch

* [ ] Fix MPI-enabled demo drivers that are broken on macOS
* [ ] Add `check_...()` calls to make sure the C/C++/Fortran/MPI compilers work
* [ ] Download CMake `oomph-lib` to a folder that has a space in the name

### Documentation

* [ ] Update `doc/the_distribution/`:
  * [ ] Remove "how to link against oomph-lib from outside the automake/autoconf framework"
  * [ ] Update third-party library version numbers or remove their information entirely
* [ ] Update instructions for how to run demo drivers
  * [ ] macOS info:
    * [ ] Add warning that the user might need to set `CMAKE_APPLE_SILICON_PROCESSOR`
    * [ ] Ask users to `brew uninstall coreutils` if it is not required (over)
    * [ ] Add documentation on how users can also just set the `CMAKE_APPLE_SILICON_PROCESSOR` environment variable (see [CMake GitLab](https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/CMakeDetermineSystem.cmake#L35))
* [ ] Adding *MPI* demo drivers (special flags)
  * [ ] Document additional `validate.sh` arg (mpi run command)
* [ ] Document how to use `CMakePresets.json` and `CMakeUserPresets.json` file
* [ ] Doc. how users can build both a `Debug` and `Release` version and easily switch between the two
  * Very simple; just build and install both build types and set `CMAKE_BUILD_TYPE` when configuring a project that runs `find_package(oomphlib)`
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
* [ ] Write a breakdown of all new features and important changes, e.g.
  * [ ] C++ implementation of `fpdiff.py`
  * [ ] `validate.sh` scripts now take the path to the root directory
* [ ] Incrementing version number (**recommend using `bumpversion.cfg` to keep git version and cmake version in sync; can make a makefile command for this 'make bump major' or 'make bump minor' or 'make bump patch'**)

## Finished

* [x] Add constraints to when we should run the self-tests (e.g. when CMakeLists.txt, .h, .c, .cc files have been edited)
* [ ] ~~Add Intel MKL BLAS support~~
* [x] Update GitHub self-tests to only run after pushing .h, .c, .cc, CMakeLists.txt code (see e.g. [here](https://github.com/scivision/mumps/blob/v5.5.1.7/.github/workflows/ci.yml)).
* [x] Discuss updating to C++17
  * [x] Required for Trilinos v14.4.0+ so will definitely be required for Trilinos v14.0
* [x] Add `oomph_option(...)` commands to specify the commandline options to the user
* [x] Add functionality to specify `--gmsh_command_line`
* [x] Sort out moving third-party libs after they've been built (will probably have issues with `oomphlibConfig.cmake`)
* [x] Update to C++17 again to patch issue in `demo_drivers/gzip/one_d_poisson`
  * Will need to test with Trilinos as it used a function removed in C++17
* [ ] ~~Add option to build with `ccache` support~~
* [x] Replace current BLAS/LAPACK with OpenBLAS
  * Looks quite complex... Appears to contain optimised kernels for different operating systems...
  * Could possibly be done by adding it as a submodule?
    * **Q:** Can we restrict it to a specific tag?
* [x] ~~Settle on how to handle versioning (e.g. should we define a `version.h`?)~~
* [x] Properly review to-do list and see if there's anything missing!
* [ ]
* Patch support for libraries using oomph-lib-built BLAS/LAPACK:
  * [x] Hypre
  * [x] Trilinos
  * N.B. Hypre and Trilinos rely on absolute paths to blas and lapack which doesn't work if they are built by oomph-lib...
  * **Possible solution:**
    * Construct BLAS/LAPACK paths
      * **Q:** Is this going to work after they've been installed?...
    * Make sure BLAS/LAPACK are built before Hypre/Trilinos
* [x] Convert external library version numbers to variables?
* [ ] ~~Address MUMPS 64-bit int (i.e. long long) compatibility issue~~
  * Not important for now

* [x] DISCUSS UPDATED SELF-TEST TRIGGER SYNTAX
* [x] Can I wipe the `config/` directory? or do we need to keep it around?
* [x] Discuss making `oomph_lib_third_party_libraries` a submodule so we can preserve the README badges (or add tests to main project; should only trigger those workflows if we edit external_distributions/)

* Add support for building and installing the following packages ourselves:
  * [x] GMP
  * [x] MPFR
* Add support for:
  * [x] MUMPS:
    * [x] Parallel version
    * [x] ~~Do we need Scotch/METIS?~~ **No**
    * [x] Sequential build
* [x] Decide whether to keep `OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY` around. **Removed**
* [x] Discuss option naming. E.g. `WANT_` or `ENABLE_` etc.
  * Useful for quick debugging of build with demo drivers
* [x] Discuss required base presets in `CMakePresets.json`
  * Being handled by MH
* [x] ~~Ninja multi-config support:~~
  * [x] ~~Test it out (e.g. building, installing and uninstalling)~~
    * [x] ~~If broken, patch any issues~~
  * [x] ~~Document how to use Ninja Multi-Config~~
    * Likely to create more problems than it'll solve for now. Can come back to this in the future...
* [x] Add control for setting number of procs to use in `mpiexec` command
  * Set `OOMPH_MPI_NUM_PROC`
* [x] Update installation paths to allow setting the installation prefix AFTER the build, e.g. `cmake --install . --prefix <PATH>`
  * [x] Change the `OOMPH_INSTALL_<XXX>_DIR` variables to not be absolute paths
  * [x] ~~Add warning to README.md that it won't work with `OOMPH_ENABLE_SYMBOLIC_LINKS_FOR_HEADERS=ON`~~
  * [x] Remove symbolic linking for headers
* [x] Set `<CONFIG>_POSTFIX` for each config type
  * Useful when looking at files in `install/` and trying to remember what build type they correspond to
* [x] Print additional helpful info during configure step, e.g. install destination
* [x] Get timing comparison for Accelerate.Framework BLAS and OpenBLAS on M1 Mac
  * Rupinder ran the meshing/adaptive_tet_meshes test
  * ~2m but OpenBLAS was 4s slower, so very little difference
* [x] Implement Python script to automatically convert `Makefile.am`s in `private/`
* [x] Remove `INSTALL.md` (installation instructions in main `README.md`)
* [x] Add `self_test/` directory to CMake build
* [x] Add support for M1 Macs
  * Specify `CMAKE_APPLE_SILICON_PROCESSOR=arm64`
* [x] Fix building of documentation
* [x] Let the user specify their own BLAS/LAPACK libraries
* [x] Document installing CMake (3.22)
* [x] Document downloading and installing oomph-lib
* [x] Document adding demo drivers
  * [x] Document changes to `validate.sh` (first arg will be `oomph-lib` root; absolute path!!)
  * [x] One oomph_add_test() per demo_driver directory to match each validate.sh
* [x] Update `OomphMPI.cmake` to avoid doing a global `link_libraries()`
* [x] ~~Use `PROJECT_IS_TOP_LEVEL` to enable/disable tests~~
* [x] Update `Doxyfile`s to set `QUIET = NO` again
* [x] Update `demo_drivers` to allow building from any level
* Add initial support for:
  * [x] `Hypre`
  * [x] `MUMPS`
  * [x] `Trilinos`
  * [x] Specifying the location of already-installed libraries
* [x] Add `oomph::generic` to `LIBRARIES` for every demo driver
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
