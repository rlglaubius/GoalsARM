# GoalsARM
HIV transmission dynamic model

## Development

### Prerequisites

* [CMake](https://cmake.org/) (>3.15) for packaging
* [Catch2](https://github.com/catchorg/Catch2) (optional) for testing. It will be fetched by CMake if you do not have it installed locally

### Build with CMake

```console
cmake -S . -B build
cmake --build build
```

To exclude compiling test code add
```console
cmake -S . -B build -DBUILD_TESTING=off
cmake --build build
```

Note that when building with testing off nothing is going to be compiled above, this is a header only library and so this will only fetch dependencies and export a target so this can be compiled downstream.

### Testing

```console
cmake -S . -B build
cmake --build build
ctest --test-dir build/tests
```

### IDE integration

It is recommended you setup a debug build target within your IDE. Ensure that any build target is added to the .gitignore.

#### CLion

Should work with CMake out of the box.

#### VSCode

* Install C/C++ extension for VSCode
* Install CMake Tools extension for VSCode

1. Select a kit, this is a toolchain for building your project. Using Ctrl+Shift+P and running CMake: Select a Kit.
2. Select a variant. This will make several variants avaialble to you, `Debug`, `Release`, `MinRelSize`, `RelWithDebInfo`. `Debug` disables optimisations and includes debug info. You'll probably only want to ever build with `Debug` in this repo. So pick this, using Ctrl+Shift+P and running `CMake: Select Variant`. You should see this compiling your code and creating a `build` directory with generated CMake files and compiled code.
3. You should now be able to build and debug by using the buttons "Build" and "Debug" buttons in the bottom of VSCode. Note that the options shown in the Top right won't work as these are trying to use vscode in build
4. You should be able to run individual tests from the testing tab, but to debug from the testing tab you'll need to add a [launch.json](https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/debug-launch.md#debug-using-a-launchjson-file).
