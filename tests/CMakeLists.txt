cmake_minimum_required(VERSION 3.15)

message(CHECK_START "Locating Catch2")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

find_package(Catch2 3 CONFIG)
if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found, fetching dynamically...")
    Include(FetchContent)

    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG        v3.4.0 # or a later release
    )

    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
else()
    message(STATUS "Catch2 found at ${Catch2_INCLUDE_DIRS}")
endif()
list(POP_BACK CMAKE_MESSAGE_INDENT)

# These tests can use the Catch2-provided main
add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE GoalsARM Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(tests)
