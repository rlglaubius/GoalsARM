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
    message(STATUS "Catch2 fetched succesfully")
else()
    message(STATUS "Catch2 found at ${Catch2_INCLUDE_DIRS}")
endif()
list(POP_BACK CMAKE_MESSAGE_INDENT)

file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

foreach(TEST_FILE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    file(RELATIVE_PATH FILE_NAME ${CMAKE_CURRENT_SOURCE_DIR} ${TEST_FILE})
    add_executable(${TEST_NAME} ${FILE_NAME})
    target_link_libraries(${TEST_NAME} PRIVATE GoalsARM Catch2::Catch2WithMain)
endforeach()

include(CTest)
include(Catch)
## Note it is important these are as a separate block for them to
## be built correctly
foreach(TEST_FILE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    catch_discover_tests(${TEST_NAME})
endforeach()
