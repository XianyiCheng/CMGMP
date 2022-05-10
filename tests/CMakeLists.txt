add_executable(test_utils test_utils.cpp)
add_executable(test_mode test_contact_mode.cpp)
add_executable(test_dart test_dartworld.cpp)
add_executable(test_kr5 test_KR5.cpp)
add_executable(test_ddhand test_ddhand.cpp)

set(example_required_libraries 
    cmgmp-utils cmgmp-contact
    )

set(example_include_dirs 
    "${PROJECT_SOURCE_DIR}" 
    "${PROJECT_SOURCE_DIR}/cmgmp/external/eigen3" 
    )

set(example_link_dirs 
    "${PROJECT_BINARY_DIR}" 
    "${PROJECT_BINARY_DIR}/cmgmp/external/dartsim"
    )

set(dart_libs dart dart-utils-urdf dart-gui) # dart-collision-bullet


target_include_directories(test_utils PRIVATE ${example_include_dirs})
# target_link_directories(test_utils PRIVATE "${PROJECT_BINARY_DIR}")
target_link_libraries(test_utils PRIVATE ${example_required_libraries})

target_include_directories(test_mode PRIVATE ${example_include_dirs})
target_link_libraries(test_mode PRIVATE ${example_required_libraries})

target_include_directories(test_dart PRIVATE ${example_include_dirs} "${PROJECT_SOURCE_DIR}/cmgmp/external/dartsim")
target_link_libraries(test_dart PRIVATE ${example_required_libraries} cmgmp-world ${dart_libs})

target_include_directories(test_kr5 PRIVATE ${example_include_dirs} "${PROJECT_SOURCE_DIR}/cmgmp/external/dartsim")
target_link_libraries(test_kr5 PRIVATE ${example_required_libraries} cmgmp-world ${dart_libs})

target_include_directories(test_ddhand PRIVATE ${example_include_dirs} "${PROJECT_SOURCE_DIR}/cmgmp/external/dartsim")
target_link_libraries(test_ddhand PRIVATE ${example_required_libraries} cmgmp-world ${dart_libs})
