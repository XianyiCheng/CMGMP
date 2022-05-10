cmake_minimum_required(VERSION 3.14)

# Source files.
set(MODUS_KINEMATICS_SOURCE
  contact_constraint.cpp
  contact_kinematics.cpp
  jacobian_body.cpp
  jacobian_floating.cpp
  jacobian_geometric.cpp
  jacobian_spatial.cpp
)

add_library(modus_kinematics SHARED ${MODUS_KINEMATICS_SOURCE})

# Define headers for this library. PUBLIC headers are used for compiling the
# library, and will be added to install paths.
target_include_directories(modus_kinematics
  PUBLIC
  # $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  # $<INSTALL_INTERFACE:include>
  # ${EIGEN3_INCLUDE_DIR}
  # ${CMAKE_INSTALL_PREFIX}/include
  ${PROJECT_SOURCE_DIR}/include
)

# Add compiler requirements here.
target_compile_features(modus_kinematics PUBLIC cxx_std_17)

# Link libraries.
# target_link_libraries(modus_kinematics modus_common modus_system gcc_s dw)
target_link_libraries(modus_kinematics modus_common modus_system gcc_s)
