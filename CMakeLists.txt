cmake_minimum_required(VERSION 3.15)
project(GoalsARM VERSION 1.0.0
                 DESCRIPTION "HIV transmission dynamic model"
                 HOMEPAGE_URL "https://github.com/avenirHealth-org/GoalsARM"
                 LANGUAGES CXX)
add_library(GoalsARM INTERFACE)

message(CHECK_START "Locating Boost")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

if(POLICY CMP0167)
   # This policy removes FindBoost module
   # in favour of using BoostConfig.cmake files
   # but it cannot locate these as reliably, so let's
   # just use the old policy for now.
   cmake_policy(SET CMP0167 OLD)
endif()
find_package(Boost 1.82 REQUIRED)
target_include_directories(GoalsARM INTERFACE ${Boost_INCLUDE_DIRS})

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
    add_subdirectory(tests EXCLUDE_FROM_ALL)
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
