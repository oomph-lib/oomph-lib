# Task list

## Urgent and important
1. Switch to mpic++ and mpicxx compilers for MPI-enabled code
2. Fix broken MPI-enabled demo drivers

## Urgent but not important
1. Empty.

## Not urgent but important
1. Add Ninja and MPI support to Workflows
2. Add code coverage.
3. Add Docker support; see https://github.com/FEniCS/dolfinx
4. Remove warnings for external sources on a target basis; see e.g. OomphWarnings.cmake.
5. Test the build of the demo_drivers with the main build (i.e. add_subdirectory(demo_drivers)).
6. Add sanitiser support (e.g. address/memory, etc.).
7. Sort out external distribution build
8. Add check_...() calls to make sure the C/C++/Fortran compiler work
9. Add (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) in conditionals where appropriate. See: https://cliutils.gitlab.io/modern-cmake/chapters/testing.html
10. Sort out building of oomph-lib with shared libraries. Currently breaks when it gets to oomph_superlu
11. Sort out *sequential* build of MUMPS
12. Sort out a subproject build of private/
13. Grab latest version of oomph-lib and add in new code.
14. Find and doc. how users can build a Debug and Release version and easily switch between the two.
15. Remove OomphPrintCompilerDefinitions.cmake.
16. Download CMake oomph-lib in a folder that has a space in the name.

## Not urgent and not important
1. Sort out external distribution build with up-to-date sources.
2. Sort out external sources build with up-to-date sources.
3. Consider replacing bash scripts with Python scripts for demo_driver tests to make them platform-independent.
4. Address MUMPS 64-bit int (i.e. long long) compatibility issue

## Finished
1. **[DONE]:** Test build and self-tests with Linux build system using GitHub Actions and output test results.
2. **[DONE]:** Change all file headers, e.g. in demo drivers ("This file forms part of oomph-lib...") so that they use the correct version and "LastChangedDate". Create a script to automate this! Also update the Copyright notice at the same time.
3. **[DONE]:** Remove set of CMAKE_OSX_DEPLOYMENT_TARGET in demo_drivers/CMakeLists.txt. Find out how else to fix the issue with deployment targets causing a lot of warnings on OSX. [Consensus: I don't care -- it works for me with the SDK on the new MacOS.]
4. **[DONE]:** Ask Matthias about the two_d_unstructured_adaptive_poisson/ validate.sh script.
5. **[DONE]:** Add a change-log. [Already have one, but maybe it's time we actually use it...]
6. **[DONE]:** Set up build of docs. [Handled by Jon's excellent Workflow script.]
7. **[DONE]:** Make sure all demo driver CMake scripts have the "#---- TESTING ---" string
8. **[DONE]:** Sort out .bumpversion.cfg.
9. **[DONE]:** Add CI support.
10. **[DONE]:** Add .clang-format to the build. Need to settle on the formatting choices. Probably worth creating clear examples of each choice.
11. **[DONE]:** Fix issue with running CMake generation step twice when FindMPI is included.
12. **[DONE]:** Create a version.h.in file. [Version info in oomph-lib-config.h]
13. **[DONE]:** Fix support for the oomph-lib-config.h.
14. **[DONE]:** Add PARANOIA and RANGE_CHECKING to possible compile options/decide on how to add it from commandline.
15. **[DONE]:** Retest new version of Scalapack and MUMPS. Patch issue with overridden install paths.
16. **[DONE]:** Sort out construction of oomph-lib-config.h.
17. **[DONE]:** Change conditional statements involving WANT_MPI to OOMPH_HAS_MPI
18. **[DONE]:** Make all "options" in main CMakeLists.txt file OOMPH_ specific so that if it is imported in, there won't be any issues with cache variables.
19. **[DONE]:** Test out demo drivers build with MPI
20. **[DONE]:** Change oomph_parmetis_* libtype to static.
21. **[DONE]:** Replace SuperLU files with an up-to-date version.
22. **[DONE]:** Address simple warning in generic/matrices.cc: "2561:2: warning: 'delete' applied to a pointer that was allocated with 'new[]'; did you mean 'delete[]'? [-Wmismatched-new-delete] delete dist_nrow_local;" (may be more than one case)
23. **[DONE]:** Sort out the mpi build! Add MPI::MPI_CXX to libraries everywhere and set MPI_RUN_COMMAND.
24. **[DONE]:** Add RESOURCE_LOCK properties to MPI-enabled drivers.
25. **[DONE]:** Update SuperLU to v5.2.2 and run demo drivers self-tests (1 fail in shell/plate; need to speak to Matthias).
26. **[DONE]:** Fix demo_drivers/biharmonic alone for now; come back and fix the issue with topologically_rectangular_domain.cc (duplicate symbols). [**A:** Needed to declare functions defined outside a class but inside the header as "inline".]
27. **[DONE]:** Think about how users can specify their own libraries to be built (but not in src/).
28. **[DONE]:** Disable SYMBOLIC_LINKS_FOR_HEADERS if not a Unix-based system.
29. **[DONE]:** For ideas on how to handle ``make check``, see [this link](https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/EmulateMakeCheck)
30. Test out whether adding WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES as a compile definition for certain executables actually produce warnings.
31. **[DONE]:** Decide on what to do with the non-mesh src/meshes sources
32. Need to make drivers default to the project C++ standard; set this in the config file.
33. **[DONE]:** Need a way to define demo driver requirements, i.e. files that it doesn't link to but it depends on, e.g. a folder of vmtk_files.
34. **[DONE]:** Sort out a subproject build of demo_drivers
35. **[DONE]:** Add a config header and export it with the library!
36. **[DONE]:** Add an "include_guard()" to CMake modules files likely to be included multiple times.
37. **[DONE]:** Specify required dependencies of the package. See: https://cmake.org/cmake/help/git-master/guide/importing-exporting/index.html#id8

