cmake_minimum_required(VERSION 3.16.3)

project(CMGMP)

# test if we can find bullet through find_package
find_package(Bullet COMPONENTS BulletMath BulletCollision MODULE QUIET)
if (BULLET_FOUND)
    message("find bullet")
else()
    message("cannot find bullet")
endif()

add_subdirectory(cmgmp)
add_subdirectory(examples)
add_subdirectory(tests)


