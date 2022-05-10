cmake_minimum_required(VERSION 3.14)
project(modus VERSION 1.0 LANGUAGES CXX)

################################################################################
# Library options.
################################################################################

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
  "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_MODUS_COMMON       "Build modus common library"        ON)
option(BUILD_MODUS_SYSTEM       "Build modus system library"        ON)
option(BUILD_MODUS_GEOMETRY     "Build modus geometry library"      ON)
option(BUILD_MODUS_KINEMATICS   "Build modus kinematics library"    OFF)
option(BUILD_MODUS_MODES        "Build modus enumeration library"   ON)
option(BUILD_MODUS_DYNAMICS     "Build modus dynamics library"      OFF)


option(MODUS_ASSERTS_ENABLED "Enable asserts" ON)
set(MODUS_DEFAULT_ASSERT_LEVEL "Error" CACHE STRING "Default assert level")
set_property(CACHE MODUS_DEFAULT_ASSERT_LEVEL PROPERTY STRINGS Warning Debug Error Fatal)

################################################################################
# External packages.
################################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)

# find_package(Eigen3 3.3 REQUIRED)
set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../eigen3")
# find_package(ortools CONFIG REQUIRED)
find_package(Qhull REQUIRED)


# include(FetchContent)
# FetchContent_Declare(
#   Sophus
#   GIT_REPOSITORY https://github.com/strasdat/Sophus.git
#   CMAKE_ARGS "-DBUILD_TESTS=Off" # bruh this doesn't work
# )
# FetchContent_MakeAvailable(Sophus)

# include(FindGit)
# find_package(Git)
# if (NOT Git_FOUND)
#     message(FATAL_ERROR "Git not found")
# endif()
# include(ExternalProject)

################################################################################
# Add subdirectory targets.
################################################################################

if (BUILD_MODUS_COMMON)
  add_subdirectory(src/modus_common)
endif()

if (BUILD_MODUS_SYSTEM)
  add_subdirectory(src/modus_system)
endif()

if (BUILD_MODUS_GEOMETRY)
  add_subdirectory(src/modus_geometry)
endif()

if (BUILD_MODUS_KINEMATICS)
  add_subdirectory(src/modus_kinematics)
endif()

if (BUILD_MODUS_MODES)
  add_subdirectory(src/modus_modes)
endif()

