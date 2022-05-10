cmake_minimum_required(VERSION 3.14)

# Source files.
set(MODUS_COMMON_SOURCE
  assert.cpp
  # backward.cpp
  # serialization.cpp
  convex_hull.cpp
  linear_algebra.cpp
  # logging.cpp
)

add_library(modus_common SHARED ${MODUS_COMMON_SOURCE})

# Define headers for this library. PUBLIC headers are used for compiling the
# library, and will be added to install paths.
target_include_directories(modus_common
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
  ${QHULL_INCLUDE_DIR}
  ${CMAKE_INSTALL_PREFIX}/include
  ${PROJECT_SOURCE_DIR}/include
)

# Add compiler requirements here.
target_compile_features(modus_common PUBLIC cxx_std_17)

# Set up assert behavior.
if (MODUS_ASSERTS_ENABLED)
    target_compile_definitions(modus_common PUBLIC MODUS_ASSERTS_ENABLED)
    if (${MODUS_DEFAULT_ASSERT_LEVEL} STREQUAL "Debug")
        target_compile_definitions(modus_common PUBLIC MODUS_DEFAULT_ASSERT_LEVEL_DEBUG)
    elseif(${MODUS_DEFAULT_ASSERT_LEVEL} STREQUAL "Warning")
        target_compile_definitions(modus_common PUBLIC MODUS_DEFAULT_ASSERT_LEVEL_WARNING)
    elseif(${MODUS_DEFAULT_ASSERT_LEVEL} STREQUAL "Error")
        target_compile_definitions(modus_common PUBLIC MODUS_DEFAULT_ASSERT_LEVEL_ERROR)
    elseif(${MODUS_DEFAULT_ASSERT_LEVEL} STREQUAL "Fatal")
    target_compile_definitions(modus_common PUBLIC MODUS_DEFAULT_ASSERT_LEVEL_FATAL)
    endif()
endif()

# Link libraries.
# target_link_libraries(modus_common gcc_s dw ${QHULLCPP_LIBRARY} ${QHULL_R_LIBRARY})
target_link_libraries(modus_common gcc_s ${QHULLCPP_LIBRARY} ${QHULL_R_LIBRARY})
