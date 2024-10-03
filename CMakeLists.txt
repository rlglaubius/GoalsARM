cmake_minimum_required(VERSION 3.15)
project(GoalsARM VERSION 1.0.0
                 DESCRIPTION "HIV transmission dynamic model"
                 HOMEPAGE_URL "https://github.com/avenirHealth-org/GoalsARM"
                 LANGUAGES CXX)
add_library(GoalsARM INTERFACE)

message(CHECK_START "Locating Boost")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

find_package(Boost 1.82 CONFIG QUIET)
if(NOT Boost_FOUND)
    message(STATUS "Boost not found, fetching Boost dynamically...")

    include(FetchContent)
    FetchContent_Declare(
            boost_math
            GIT_REPOSITORY https://github.com/boostorg/math.git
            GIT_TAG boost-1.82.0
            GIT_SHALLOW TRUE
    )

    FetchContent_Declare(
            boost_multi_array
            GIT_REPOSITORY https://github.com/boostorg/multi_array.git
            GIT_TAG boost-1.82.0
            GIT_SHALLOW TRUE
    )

    FetchContent_MakeAvailable(boost_math boost_multi_array)
    target_include_directories(GoalsARM INTERFACE ${boost_math_SOURCE_DIR}/include ${boost_multi_array_SOURCE_DIR}/include)
else()
    message(STATUS "Boost found at ${Boost_INCLUDE_DIRS}")
    target_include_directories(GoalsARM INTERFACE ${Boost_INCLUDE_DIRS})
endif()
list(POP_BACK CMAKE_MESSAGE_INDENT)

# Add library
target_include_directories(GoalsARM
        INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(GoalsARM INTERFACE cxx_std_20)

# Enable testing if this project is built standalone (not included as a subproject)
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Install the library and export CMake targets
include(GNUInstallDirs)
install(TARGETS GoalsARM
        EXPORT goals-export
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(EXPORT goals-export
        FILE GoalsARMTargets.cmake
        NAMESPACE goals::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoalsARM
)