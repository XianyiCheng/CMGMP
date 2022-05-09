set(dart_libs dart dart-utils-urdf dart-gui)

set(example_required_libraries 
    cmgmp-utils cmgmp-contact cmgmp-world cmgmp-search 
    ${dart_libs}
    )

set(example_include_dirs 
    "${PROJECT_SOURCE_DIR}" 
    "${PROJECT_SOURCE_DIR}/cmgmp/external/eigen3" 
    "${PROJECT_SOURCE_DIR}/cmgmp/external/dartsim"
    )

set(example_link_dirs 
    "${PROJECT_BINARY_DIR}" 
    "${PROJECT_BINARY_DIR}/cmgmp/external/dartsim"
    )




include_directories(${example_include_dirs})
link_libraries(${example_required_libraries})

add_executable(kr5 KR5.cpp)
add_executable(push_along_walls push_along_walls.cpp)
add_executable(ddhand_bookself ddhand_bookshelf.cpp)
add_executable(ddhand_flip ddhand_flip.cpp)
add_executable(ddhand_slideup ddhand_slideup.cpp)
