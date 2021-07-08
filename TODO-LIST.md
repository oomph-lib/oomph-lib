# Task list

## Urgent and important
1. Update oomphlibConfig.cmake.in to store CMAKE_INSTALL_PREFIX and add that to the CMAKE_PREFIX_PATH(?) so that it can be found if installed in non-standard locations (but might not make it relocatable)
2. Switch to mpic++ and mpicxx compilers for MPI-enabled code
4. Fix broken MPI-enabled demo drivers

## Urgent but not important
1. Empty.

## Not urgent but important
1. Make sure all demo driver CMake scripts have the "#---- TESTING ---" string
1. Remove set of CMAKE_OSX_DEPLOYMENT_TARGET in demo_drivers/CMakeLists.txt
2. Use glibtool.
3. Add code coverage.
4. Add a change-log.
5. Sort out .bumpversion.cfg.
7. Set up build of docs.
8. Add Docker support; see https://github.com/FEniCS/dolfinx
9. Get Git info from CMake file.
10. Add a .circleci/ci.yml for CI.
11. Find out how else to fix the issue with deployment targets causing a lot of warnings on OSX.
12. Ask Matthias about how he intended the two_d_unstructured_adaptive_poisson/ validate.sh script; should I not pass the value of MPI_VARIABLENP_RUN_COMMAND to the script? Either way, which flag is meant to be called for the "sample_point_container_flag" and where is it supplied? (Because I don't know where on earth he specifies it.)
13. Decide on whether to provide an (optional?) VALIDATE_SH_ARGS argument to oomph_add_test(...) to provide more flexibility (for the MPI_VARIABLE_NP) demo driver
14. Test the build of the demo_drivers with the main build (i.e. add_subdirectory(demo_drivers)).
15. Address the huge number of warnings from Clang.
16. Sort out external distribution build
17. Sort out external distribution build with up-to-date sources.
18. Sort out external sources build with up-to-date sources.
19. Add check_...() calls to make sure the C/C++/Fortran compiler work
20. Add (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) in conditionals where appropriate. See: https://cliutils.gitlab.io/modern-cmake/chapters/testing.html
21. Sort out building of oomph-lib with shared libraries. Currently breaks when it gets to oomph_superlu
22. Sort out sequential build of MUMPS
23. Sort out a subproject build of private/
24. Grab latest version of oomph-lib and add in new code.
25. Change all file headers, e.g. in demo drivers ("This file forms part of oomph-lib...") so that they use the correct version and "LastChangedDate". Create a script to automate this! Also update the Copyright notice at the same time.
26. Test build and self-tests with Linux build system using GitHub Actions and output test results.
27. Download CMake oomph-lib in a folder that has a space in the name.
28. Find out how users can build a Debug and Release version and easily switch between the two.
29. Remove OomphWarnings.cmake.
30. Remove OomphPrintCompilerDefinitions.cmake.

## Not urgent and not important
1. Consider replacing bash scripts with Python scripts for demo_driver tests to make them platform-independent.
2. Address MUMPS 64-bit int (i.e. long long) compatibility issue

## Finished
1. **[DONE]:** Add .clang-format to the build. Need to settle on the formatting choices. Probably worth creating clear examples of each choice.
2. **[DONE]:** Fix issue with running CMake generation step twice when FindMPI is included.
3. **[DONE]:** Create a version.h.in file. [Version info in oomph-lib-config.h]
4. **[DONE]:** Fix support for the oomph-lib-config.h.
5. **[DONE]:** Add PARANOIA and RANGE_CHECKING to possible compile options/decide on how to add it from commandline.
6. **[DONE]:** Retest new version of Scalapack and MUMPS. Patch issue with overridden install paths.
7. **[DONE]:** Sort out construction of oomph-lib-config.h.
8. **[DONE]:** Change conditional statements involving WANT_MPI to OOMPH_HAS_MPI
9. **[DONE]:** Make all "options" in main CMakeLists.txt file OOMPH_ specific so that if it is imported in, there won't be any issues with cache variables.
10. **[DONE]:** Test out demo drivers build with MPI
11. **[DONE]:** Change oomph_parmetis_* libtype to static.
12. **[DONE]:** Replace SuperLU files with an up-to-date version.
13. **[DONE]:** Address simple warning in generic/matrices.cc: "2561:2: warning: 'delete' applied to a pointer that was allocated with 'new[]'; did you mean 'delete[]'? [-Wmismatched-new-delete] delete dist_nrow_local;" (may be more than one case)
14. **[DONE]:** Sort out the mpi build! Add MPI::MPI_CXX to libraries everywhere and set MPI_RUN_COMMAND.
15. **[DONE]:** Add RESOURCE_LOCK properties to MPI-enabled drivers.
16. **[DONE]:** Update SuperLU to v5.2.2 and run demo drivers self-tests (1 fail in shell/plate; need to speak to Matthias).
17. **[DONE]:** Fix demo_drivers/biharmonic alone for now; come back and fix the issue with topologically_rectangular_domain.cc (duplicate symbols). [**A:** Needed to declare functions defined outside a class but inside the header as "inline".]
18. **[DONE]:** Think about how users can specify their own libraries to be built (but not in src/).
19. **[DONE]:** Disable SYMBOLIC_LINKS_FOR_HEADERS if not a Unix-based system.
20. **[DONE]:** For ideas on how to handle ``make check``, see [this link](https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/EmulateMakeCheck)
21. Test out whether adding WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES as a compile definition for certain executables actually produce warnings.
22. **[DONE]:** Decide on what to do with the non-mesh src/meshes sources
23. Need to make drivers default to the project C++ standard; set this in the config file.
24. **[DONE]:** Need a way to define demo driver requirements, i.e. files that it doesn't link to but it depends on, e.g. a folder of vmtk_files.
25. **[DONE]:** Sort out a subproject build of demo_drivers
26. **[DONE]:** Add a config header and export it with the library!
27. **[DONE]:** Add an "include_guard()" to CMake modules files likely to be included multiple times.
28. **[DONE]:** Specify required dependencies of the package. See: https://cmake.org/cmake/help/git-master/guide/importing-exporting/index.html#id8

