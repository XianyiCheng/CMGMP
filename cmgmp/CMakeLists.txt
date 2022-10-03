add_subdirectory(external)

find_library(GLPKLIB glpk REQUIRED)

set(cmgmp_required_libraries 
    ${GLPKLIB}
    )

set(modus_libs modus_common modus_modes)

set(dart_libs dart dart-utils-urdf dart-gui dart-collision-bullet)

add_definitions(-DSRC_DIR="${PROJECT_SOURCE_DIR}")

set(cmgmp_include_dirs 
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/cmgmp/external/eigen3"
    "${PROJECT_BINARY_DIR}/cmgmp/external/modus/include"
    "${PROJECT_BINARY_DIR}/cmgmp/external/dartsim"
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

target_include_directories(cmgmp-contact PRIVATE ${cmgmp_include_dirs})
target_link_libraries(cmgmp-contact PRIVATE cmgmp-utils ${cmgmp_required_libraries} ${modus_libs})

add_library(cmgmp-world
    worlds/DartWorld.cpp
    worlds/PlanningWorld.cpp
    manipulators/DartPointManipulator.cpp
    # manipulators/DartParallelJawGripper.cpp
    manipulators/ddhand/two_link_kinematics.cpp
    manipulators/DartDDHandFreeY.cpp
    manipulators/DartDDHand.cpp
    manipulators/DartKR5.cpp
    dart_utils/dart_utils.cpp
    )

target_include_directories(cmgmp-world PRIVATE ${cmgmp_include_dirs})
target_link_libraries(cmgmp-world PRIVATE cmgmp-utils cmgmp-contact ${dart_libs} ${modus_libs})

add_library(cmgmp-search
    search/rrtplanner.cpp
    search/tree.cpp)

target_include_directories(cmgmp-search PRIVATE ${cmgmp_include_dirs})
target_link_libraries(cmgmp-search PRIVATE cmgmp-utils cmgmp-contact cmgmp-world ${cmgmp_required_libraries} ${dart_libs} ${modus_libs})
