cmake_minimum_required(VERSION 3.14)

# Source files.
set(MODUS_SOURCE
    enumerate_modes.cpp
    # preprocess.cpp
    # constraints.cpp

    # Preprocessing
    preprocess/matroid_preprocessor.cpp

    # Geometry
    geometry/arrangements.cpp
    geometry/incidence_graph.cpp
    # geometry/linear_algebra.cpp
    geometry/interior_point.cpp
    geometry/convex_hull.cpp
    geometry/halfspace_intersection.cpp
)

add_library(modus_modes SHARED ${MODUS_SOURCE})

# Define headers for this library. PUBLIC headers are used for compiling the
# library, and will be added to install paths.
target_include_directories(modus_modes 
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${EIGEN3_INCLUDE_DIR}
    ${QHULL_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/include 
    )
    #${ORTOOLS_INCLUDE_DIR})

# Add compiler requirements here.
target_compile_features(modus_modes PUBLIC cxx_std_17)

# Link libraries.
target_link_libraries(modus_modes modus_common modus_kinematics ${QHULLCPP_LIBRARY} ${QHULL_R_LIBRARY} )#ortools::ortools)
