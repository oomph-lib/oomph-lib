# Notes

This is a scrapbook for now. It will get tidied at the "end", i.e. before a beta
release.

- Important CMake files are in the `cmake/` directory. These files are "included"
  into the library scope by appending the path to the `CMAKE_MODULE_PATH` variable.

## Overriding `CMakePresets.json` options

To override options in `CMakePresets.json` with your own, create a `CMakeUserPresets.json` file. CMake will read this file after reading `CMakePresets.json` so you can interpret such a file as being "appended" to the bottom of the `CMakePresets.json` file. However, there is one field that you must make sure to duplicate: the `version` field. Once done, you can then "inherit" the fields of configuration presets provided in `CMakePresets.json`. Below is an example of a custom configuration preset:

**`CMakeUserPresets.json`**

```json
{
  "version": 3,
  "configurePresets": [
    {
      "name": "Custom Config",
      "inherits": "base",
      // The name displated
      "displayName": "Custom Config Default",
      "generator": "Ninja",
      "binaryDir": "${sourceDir}/build",
      // Overrides of CMake configuration options
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "Debug",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_PREFIX_PATH": "/usr/local/Cellar/openblas/0.3.17;${sourceDir}/build/",
        "OOMPH_ENABLE_MPI": "ON",
        "OOMPH_USE_OPENBLAS": "ON"
      },
      // Provide your custom C/C++ compilers
      "environment": {
        "CC": "/usr/bin/clang",
        "CXX": "/usr/bin/clang++"
      },
      "warnings": {
        "unusedCli": false
      }
    }
  ],
  "buildPresets": [
    {
      "name": "Custom Config",
      "configurePreset": "Custom Config",
      // Build with 4 jobs, e.g. 'ninja -j4'/'make -j4'
      "jobs": 4
    }
  ],
  "testPresets": [
    {
      "name": "Custom Config",
      "configurePreset": "Custom Config",
      "output": {
        "outputOnFailure": true,
        "labelSummary": false
      },
      "execution": {
        "timeout": 10000,
        "noTestsAction": "error",
        "jobs": 4
      }
    }
  ]
}
```

For more information, see the [CMake presets documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

## Running a `validate.sh` script manually (without `CTest`)

- There must be at least one argument to each `validate.sh` script: the absolute path to the oomph-lib root directory. The root directory path is used to find the `fpdiff.py` script, `validate_ok_count` and sometimes other validation tools.

## Reuse installed (third-party) libraries

- Once a third-party library has been installed, you may want to move the installation from the `build/` folder **[not sure this is true, I think they get installed to the same location as the library unless we've explicitly set `<INSTALL_DIR>` with `ExternalProject`/`FetchContent`]** to some other folder of your choice. If you move the installation then wipe the `build/` folder, you can tell `oomph-lib` when you next build it that you already have an installation of the same third-party library. To do you'll need to specify a path to it with the variable `-DOOMPH_WITH_FOO=BAR`, where `FOO` is one of the supported values and `BAR` is the path to the installation (it should contain a `include` and `lib` subdirectory). For example, if you moved `Trilinos`
```bash
cmake -G Ninja -DOOMPH_WITH_TRILINOS="~/trilinos_install" -B build
```
In general, it is rather painful to specify this flag so we recommend specifying this value in your `CMakeUserPresets.json` file.



## Defining libraries

- Library directory is assumed to be the same as the library name. If not,
  set the correct directory name with `INCLUDE_SUBDIRECTORY` when calling
  `oomph_library_config()`


## Steps to adding the Google benchmark library

**FIXME:** Update these instructions.

### Initial steps

- Add `FetchContent` recipe to download the `benchmark` library.
- Create a `benchmark/` directory.
- Add a `CMakeLists.txt` file to the folder with the obvious boilerplate (`cmake_minimum_required(...)` and `project(...)`).
- If the library build includes the `benchmarking/` directory, add a `find_package(oomphlib REQUIRED)` call.
- `include(...)` the module for downloading the benchmarking library.
- To ensure the `benchmark` library can find the `OomphGetGoogleBenchmark.cmake` script after installation, make sure it gets distributed with the library to `lib/cmake/oomphlib/`.
-
### Build-time compatibility

Need
```bash
# In benchmark/generic/ directory
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=../../build -B build
cd build
ninja
```

- Add `FetchContent` recipe to download the `benchmark` library.
- Create a `benchmark/` directory.
- Add a CMake script that defines the project.
- If the benchmarking is NOT part of the library build:
    - Add a `find_package(...)` call for the `oomphlib` package
- Include the module for downloading the benchmarking library.
- To ensure the `benchmark` library can find the `OomphGetGoogleBenchmark.cmake` script after installation, make sure it gets distributed with the library to `lib/cmake/oomphlib/`.
- TO FIND OUT: Should the `.cmake` files that get installed be added to the `CMAKE_MODULE_PATH`? ACTUALLY, probably just needs to be the current `...LIST_DIR`, so it works at build time and install time
- Should probably store the path to the `benchmark/` directory and call `find_dependency()` after in the main config file and use the download cmake file otherwise
- Need to include the module

### Install-time compatibility
**Things to check:**
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