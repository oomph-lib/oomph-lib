# Notes

This is a scrapbook for now. It will get tidied at the "end", i.e. before a beta
release.

- Important CMake files are in the `cmake/` directory. These files are "included"
  into the library scope by appending the path to the `CMAKE_MODULE_PATH` variable.

## Preserve previously-installed (third-party) libraries

- Once a third-party library has been installed, you may want to move it's installation from the `build/` folder to some other folder of your choice. If you wipe the `build/` folder and want to tell `oomph-lib` that you already have an installation of the third-party library that you want to use again, you'll need to specify a path to it with the variable `-DOOMPH_WITH_FOO=BAR`, where `FOO` is one of the supported values and `BAR` is the path to the installation (it should contain a `include` and `lib` subdirectory). For example, if you moved `Trilinos`
```bash
cmake -G Ninja -DOOMPH_WITH_TRILINOS="~/trilinos_install" -B build
```
In general, it is rather painful to specify this flag so we recommend you specify this value in your `CMakeUserPresets.json` file.



## Defining libraries

- Library directory is assumed to be the same as the library name. If not,
  set the correct directory name with `INCLUDE_SUBDIRECTORY` when calling
  `oomph_library_config()`


## Steps to adding the Google benchmark library

### Initial steps

- Add FetchContent recipe to download the `benchmark` library.
- Create a `benchmark/` directory.
- Add a CMakeLists.txt file to the folder with the obvious boilerplate (`cmake_minimum_required(...)` and `project(...)`).
- If the library build includes the benchmarking directory, add a `find_package(oomphlib REQUIRED)` call.
- `include(...)` the module for downloading the benchmarking library.
- To ensure the benchmarking library can find the GetGoogleBenchmark script after installation, make sure it gets distributed with the library to lib/cmake/oomphlib/.
-
### Build-time compatibility

Need
```bash
# In benchmark/generic/ directory
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=../../build -B build
cd build
ninja
```

- Add FetchContent recipe to download the benchmark library.
- Create a benchmark directory.
- Add a CMake script that defines the project.
- If the benchmarking is NOT part of the library build:
    - Add a find_package call for the oomphlib package
- Include the module for downloading the benchmarking library.
- To ensure the benchmarking library can find the GetGoogleBenchmark script after installation, make sure it gets distributed with the library to lib/cmake/oomphlib/.
- TO FIND OUT: Should the .cmake files that get installed be added to the CMAKE_MODULE_PATH? ACTUALLY, probably just needs to be the current LIST_DIR, so it works at build time and install time
- Should probably store the path to the benchmark directory and call find_dependency() after in the main config file and use the download cmake file otherwise
- Need to include the module

### Install-time compatibility
Things to check:
- [ ] Can it be used at build time?
- [ ] Can it be used after the library is installed?
Both require updating the CMake module path.

### Questions

- **When 'benchmark' is built with the library:**
  - [ ] Will it work at build-time?
    - Yep! Don't have to do anything -- piece of cake!

- **When 'benchmark' is _not_ built with the library:**
  - [ ] Will it work at build-time?
    - Yep! You need to specify the path to the `build` directory though. *sad face*
  - [ ] Will it work at install-time?
    - Oh yes! And it's as simple as just calling `find_package(oomphlib REQUIRED)`.


- **Once downloaded, will it try to redownload the 'benchmark' library?**