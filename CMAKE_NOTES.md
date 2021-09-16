## Notes

This is a scrapbook for now. It will get tidied at the "end", i.e. before a beta
release.

- Important CMake files are in the `cmake/` directory. These files are "included"
  into the library scope by appending the path to the `CMAKE_MODULE_PATH` variable.

### Defining libraries

- Library directory is assumed to be the same as the library name. If not,
  set the correct directory name with `INCLUDE_SUBDIRECTORY` when calling
  `oomph_library_config()`
