cmake_minimum_required(VERSION 3.14)

# Source files.
set(MODUS_SYSTEM_SOURCE
  system.cpp
  body.cpp
  transform.cpp
)

add_library(modus_system SHARED ${MODUS_SYSTEM_SOURCE})

target_include_directories(modus_system
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/include
  )

# Add compiler requirements here.
target_compile_features(modus_system PUBLIC cxx_std_17)

# Link libraries.
target_link_libraries(modus_system modus_common)

# set_property(TARGET modus_system PROPERTY POSITION_INDEPENDENT_CODE ON)
