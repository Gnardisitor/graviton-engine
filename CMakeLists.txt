cmake_minimum_required(VERSION 3.10.0)
project(graviton VERSION 1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CORE_DIR src)
set(INCLUDE_DIR include)

set(CORE_SOURCES 
    ${CORE_DIR}/main.cpp
    ${CORE_DIR}/engine.cpp
    ${CORE_DIR}/solver.cpp
    ${CORE_DIR}/system.cpp
)

add_executable(${PROJECT_NAME} ${CORE_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND WIN32)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES)
    target_compile_options(${PROJECT_NAME} PRIVATE
        /W4
        /O2
        /EHsc
        /clang:-march=native
        /clang:-O3
        /clang:-ffast-math
        /clang:-fopenmp
    )
    
    get_filename_component(CLANG_BIN_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    set(OMP_DLL "${CLANG_BIN_DIR}/libomp.dll")
    
    target_link_libraries(${PROJECT_NAME} PRIVATE libomp)

    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OMP_DLL}"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
        COMMENT "Copying libomp.dll"
    )
else()
    target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall -Wextra -Wpedantic
        -O3
        -ffast-math
        -fopenmp
        -march=native
    )
    find_package(OpenMP REQUIRED)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()