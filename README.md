# GoalsARM
HIV transmission dynamic model

## Development

### Prerequisites

* [C++ compilation toolchain](https://code.visualstudio.com/docs/cpp/config-msvc) can be installed via Visual Studio, see step 3 from link
* [CMake](https://cmake.org/) (>=3.15) for packaging.
* [Boost](https://www.boost.org/users/download/) (>=1.82, <=1.85>) installed. The easiest way I find to do this, on windows, is to install one of the prebuilt binaries linked from the download page.
* [Catch2](https://github.com/catchorg/Catch2) (optional) for testing. It will be fetched by CMake if you do not have it installed locally.

### Build with CMake

First configure the project
```console
cmake -S . -B build
```

This resolves dependencies and generates build system files (Makefiles, Ninja files or VS project files) in the `build` directory. It doesn't compile any code. Then build the project

```console
cmake --build build
```

This compiles the code and produces outputs (executables, libraries). This will invoke the appropriate build tool (make, ninja or msbuild).

Note you only need to re-run the configuration step if you modify CMakeLists.txt, or want to set some CMake variable like the build type. If you have just updated source code, just running `cmake --build build` should be enough.

To exclude compiling test code add
```console
cmake -S . -B build -DGOALS_BUILD_TESTING=off --fresh
cmake --build build
```

Note that if switching between version with tests and without tests, you'll need to reconfigure with the `--fresh` command to ensure that CMake doesn't use cached files.

Note also, when building with testing turned off nothing is going to be compiled, this is a header only library and so this will only fetch dependencies and export a target so this can be compiled downstream.

### Testing

```console
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

To add a new test, you can extend the existing `tests/tests.cpp` or any new file within the `tests` directory with the `.cpp` extension will be picked up automatically.

### IDE integration

#### CLion

Should work with CMake out of the box.

#### VSCode

* Install C/C++ extension from Microsoft
* Install CMake Tools extension from Microsoft

1. Open the project in VSCode, you should be prompted to "Select a Kit for GoalsARM". A kit is the toolchain for compiling your project. You can select a specific kit or "unspecified" to let CMake decide for you.
1. (Optionally) Select a variant. VSCode will make several build variants avaialble to you, `Debug`, `Release`, `MinRelSize`, `RelWithDebInfo`. `Debug` disables optimisations and includes debug info and should be selected by default. You can select one of the others using Ctrl+Shift+P and running `CMake: Select Variant`. 
1. You should see VSCode compiling your code automatically. It will create a `build` directory with generated files and compiled code.
1. You should now be able to build and debug by using the "Build" and "Debug" buttons in the bottom bar of VSCode. Note that any compilation buttons shown in the top right won't work as these are not being compiled using the CMake configuration.
1. You should be able to run individual tests from the testing tab, but to debug from the testing tab you'll need to add a [launch.json](https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/debug-launch.md#debug-using-a-launchjson-file).
