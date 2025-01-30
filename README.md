# GoalsARM
HIV transmission dynamic model

## Development

### Prerequisites

* [C++ compilation toolchain](https://code.visualstudio.com/docs/cpp/config-msvc) can be installed via Visual Studio, see step 3 from link
* [CMake](https://cmake.org/) (>=3.15) for packaging.
* [Boost](https://www.boost.org/users/download/) (>=1.82, <=1.85>) installed. The easiest way I find to do this, on windows, is to install one of the prebuilt binaries linked from the download page.
* [Catch2](https://github.com/catchorg/Catch2) (optional) for testing. It will be fetched by CMake if you do not have it installed locally.

### Build with CMake

This project uses [CMake](https://cmake.org/) & [CMakePresets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html) to manage different builds. CMake can manage several stages of dev lifecyle. We will use it to manage configure, build & test. 
* `configure` generates build system files (Makefiles, Ninja files or VS project files) in the `build` directory, and locates dependencies. It does not compile any code.
*  `build` this compiles the code and produces outputs (executables, libraries). This will invoke the appropriate build tool (make, ninja or msbuild).
*  `test` runs tests with `ctest`.

We have the following presets configured for each stage
* Configure
    * default
    * ci
* build
    * debug
    * release
* test
    * debug
    * release

You can use these as follows

```console
# Configure
cmake --preset=default
# Build
cmake --build --preset=debug
# Test
ctest --preset=debug
```

Note you only need to re-run the configuration step if you modify CMakeLists.txt, or want to set some CMake variable like the build type. If you have just updated source code, just running `cmake --build build` should be enough.


### Testing

Testing uses [Catch2](https://github.com/catchorg/Catch2) and is run via ctest. To add a new test, you can extend the existing `tests/tests.cpp` or any new file within the `tests` directory with the `.cpp` extension will be picked up automatically.

### IDE integration

#### CLion

Should work with CMake out of the box.

#### VSCode

* Install C/C++ extension from Microsoft
* Install CMake Tools extension from Microsoft

1. At the time of writing, VSCode's CMake extension does not [pick up env variables in preset config](https://github.com/microsoft/vscode-cmake-tools/issues/4253), so we need to work around it. Copy the `cmake/CMakeUserPresets.json` from the `cmake` directory into the root. User presets are a CMake tool for individual users to write specific overrides of configuration in `CMakePresets.json`. Copy the file over, and in the `inherit` block ensure that this is inheriting from the correct presets for your system `cmake/CMakeWindowsPresets.json` on Windows.
1. Open the project in VSCode.
1. Specify the configure "preset" by opening command palette (Ctrl+Shift+P) and running "CMake: Select Configure Preset" and select "default".
1. You should see VSCode compiling your code automatically. It will create a `build` directory with generated files and compiled code. You can turn off automatic compilation and do it manually if you prefer.
1. (Optionally) Select the build & test presets by using command palette and running "CMake: Select Build Preset" and "CMake: Select Test Preset". Choose either "vscode-debug", to build with debug symbols and no optimisation, or "vscode-release" to build with optimisation and without debug symbols. For most uses you probably want the "vscode-debug" build. It will work without selecting one of these presets, and will build in debug mode.
1. You should now be able to build and debug by using the "Build" and "Debug" buttons in the bottom bar of VSCode. Note that any compilation buttons shown in the top right won't work as these are not being compiled using the CMake configuration.
1. You should be able to run individual tests from the testing tab, but to debug from the testing tab you'll need to add a [launch.json](https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/debug-launch.md#debug-using-a-launchjson-file).
