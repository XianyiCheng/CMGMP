# set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no-as-needed -pg")

macro(add_benchmark TESTNAME)
add_executable(${TESTNAME} ${ARGN})
add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
target_link_libraries(${TESTNAME} modus dw ${GPERFTOOLS_PROFILER})
endmacro()

macro(add_cuda_benchmark TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    # target_link_libraries(${TESTNAME} gtest gmock gtest_main dw)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
    # set_target_properties(${TESTNAME} PROPERTIES FOLDER test)
    set_target_properties(${TESTNAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(${TESTNAME} PRIVATE mesh_geom dw)
    target_compile_options(${TESTNAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        --expt-relaxed-constexpr
        --relocatable-device-code=true
        >)
endmacro()

# add_cuda_benchmark(bench_query bench_query.cu)
# add_cuda_benchmark(bench_bvh bench_bvh.cu)

# add_benchmark(profile_arrangements profile_arrangements.cpp)
add_benchmark(bench_subspace bench_subspace.cpp)
add_benchmark(bench_arrangements bench_arrangements.cpp)