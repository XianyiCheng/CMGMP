add_subdirectory(external)

find_library(GLPKLIB glpk REQUIRED)

set(cmgmp_required_libraries 
    ${GLPKLIB}
    )

set(modus_libs modus_common modus_modes)

set(dart_libs dart dart-utils-urdf dart-gui) # dart-collision-bullet

set(cmgmp_include_dirs 
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/cmgmp/external/eigen3"
    )

set(cmgmp_link_dirs 
    "${PROJECT_BINARY_DIR}" 
    # "${PROJECT_BINARY_DIR}/cmgmp/external/dartsim"
    # "${PROJECT_BINARY_DIR}/cmgmp/external/modus/modus_common"
    # "${PROJECT_BINARY_DIR}/cmgmp/external/modus/modus_modes"
    )


add_library(cmgmp-utils 
    utilities/sample.cpp 
    utilities/utilities.cpp 
    utilities/sample_grasp.cpp )
target_include_directories(cmgmp-utils PRIVATE ${cmgmp_include_dirs})
target_link_libraries(cmgmp-utils PRIVATE ${cmgmp_required_libraries})

add_library(cmgmp-contact 
    contacts/contact_kinematics.cpp 
    contacts/contact_constraints.cpp 
    contacts/contact_mode_enumeration.cpp
    contacts/contact_utils.cpp
    )

target_include_directories(cmgmp-contact PRIVATE ${cmgmp_include_dirs} "${PROJECT_BINARY_DIR}/cmgmp/external/modus/include")
# target_link_directories(cmgmp-contact PRIVATE cmgmp_link_dirs)
target_link_libraries(cmgmp-contact PRIVATE cmgmp-utils ${cmgmp_required_libraries} ${modus_libs})

add_library(cmgmp-world
    worlds/DartWorld.cpp
    manipulators/DartPointManipulator.cpp
    dart_utils/dart_utils.cpp
    )

target_include_directories(cmgmp-world PRIVATE ${cmgmp_include_dirs} "${PROJECT_BINARY_DIR}/cmgmp/external/dartsim")
target_link_libraries(cmgmp-world PRIVATE cmgmp-utils cmgmp-contact ${dart_libs})