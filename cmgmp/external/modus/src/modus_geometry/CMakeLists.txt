cmake_minimum_required(VERSION 3.14)

# Source files.
set(MODUS_GEOMETRY_SOURCE
  bounds.cpp
  distance.cpp
  edge.cpp
  vertex.cpp
  face.cpp
  halfedge.cpp
  halfedgemesh.cpp
  box.cpp
  plane.cpp
)

add_library(modus_geometry SHARED ${MODUS_GEOMETRY_SOURCE})

target_include_directories(modus_geometry
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/include
  )

# Add compiler requirements here.
target_compile_features(modus_geometry PUBLIC cxx_std_17)

# Link libraries.
target_link_libraries(modus_geometry modus_common)

# set_property(TARGET modus_geometry PROPERTY POSITION_INDEPENDENT_CODE ON)

