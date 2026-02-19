cmake_minimum_required(VERSION 3.10.0)
project(graviton VERSION 1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CORE_DIR src)
set(INCLUDE_DIR include)
set(TEST_DIR test)

set(CORE_SOURCES 
    ${CORE_DIR}/engine.cpp
    ${CORE_DIR}/solver.cpp
    ${CORE_DIR}/system.cpp
)

add_library(graviton_engine SHARED ${CORE_SOURCES})
target_include_directories(graviton_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIR})

set_target_properties(graviton_engine PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    SOVERSION 1
    VERSION 1.0.0
)

add_executable(${PROJECT_NAME} ${CORE_DIR}/main.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE graviton_engine)

target_compile_options(${PROJECT_NAME} PRIVATE
    -Wall -Wextra -Wpedantic
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${PROJECT_NAME} PRIVATE -O0)
    target_compile_options(graviton_engine PRIVATE -O0)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE
        -O3
        -ffast-math
        -fopenmp
        -march=native
    )
    target_compile_options(graviton_engine PRIVATE
        -O3
        -ffast-math
        -fopenmp
        -march=native
    )
    find_package(OpenMP REQUIRED)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(graviton_engine PRIVATE OpenMP::OpenMP_CXX)
endif()

enable_testing()

add_executable(graviton_test ${TEST_DIR}/test.cpp)
target_include_directories(graviton_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIR})
target_link_libraries(graviton_test PRIVATE graviton_engine)

target_compile_options(graviton_test PRIVATE -Wall -Wextra -Wpedantic)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(graviton_test PRIVATE -O3 -fopenmp -march=native)
    find_package(OpenMP REQUIRED)
    target_link_libraries(graviton_test PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME Euler COMMAND graviton_test 0)
add_test(NAME CauchyEuler COMMAND graviton_test 1)
add_test(NAME PositionVerlet COMMAND graviton_test 2)
add_test(NAME VelocityVerlet COMMAND graviton_test 3)
add_test(NAME RK4 COMMAND graviton_test 4)
add_test(NAME RK45 COMMAND graviton_test 5)
add_test(NAME PEFRL COMMAND graviton_test 6)