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
        "OOMPH_ENABLE_USE_OPENBLAS": "ON"
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
