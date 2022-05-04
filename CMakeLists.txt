cmake_minimum_required(VERSION 3.16.3)

project(CMGMP)



# add_subdirectory(cmgmp/external/dartsim)
find_package(Bullet COMPONENTS BulletMath BulletCollision MODULE QUIET)
if (BULLET_FOUND)
    message("find bullet")
else()
    message("cannot find bullet")
endif()

# add_subdirectory(cmgmp/external/dartsim)
add_subdirectory(cmgmp)
add_subdirectory(tests)
add_subdirectory(examples)

set(cmgmp_libs 
    cmgmp-utils 
    )

