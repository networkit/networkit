cmake_minimum_required(VERSION 3.0)
project(networkit CXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# BUILD OPTIONS
option(NETWORKIT_OPENMP "Build with OpenMP support (requires OpenMP)" ON)
option(NETWORKIT_BUILD_STATIC "Build as static library" ON)
option(NETWORKIT_BUILD_SHARED "Build as shared library" ON)
option(NETWORKIT_BUILD_TESTS "Build with tests (requires GTest); implies NETWORKIT_BUILD_SHARED" OFF)
option(NETWORKIT_CXX14 "Build with standard C++14" ON)

if (NETWORKIT_BUILD_TESTS)
    if (NOT NETWORKIT_BUILD_SHARED)
        message("Enable NETWORKIT_BUILD_SHARED as implied by NETWORKIT_BUILD_TESTS")
    endif()
    set(NETWORKIT_BUILD_SHARED ON)
endif()

if (NOT NETWORKIT_BUILD_STATIC AND NOT NETWORKIT_BUILD_SHARED)
    message(SEND_ERROR "Select at least one of NETWORKIT_BUILD_STATIC or NETWORKIT_BUILD_SHARED")
endif()

if (NOT CMAKE_BUILD_TYPE)
    message("Use Release Build Type as default")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(NETWORKIT_BUILD_TYPE "Opt")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(NETWORKIT_BUILD_TYPE "Dbg")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Profiling")
    set(NETWORKIT_BUILD_TYPE "Pro")
endif()
message("Detected NETWORKIT_BUILD_TYPE: ${NETWORKIT_BUILD_TYPE}")

if (NETWORKIT_BUILD_TESTS AND NOT "${NETWORKIT_BUILD_TYPE}" STREQUAL "Dbg")
    message(WARNING "Build tests without debug build may disable checks")
endif()

if (NETWORKIT_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message("OpenMP found -- OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    elseif()
        message(FATAL_ERROR "No CXX support for OpenMP found")
    endif()
endif()

if (NETWORKIT_CXX14)
    set (CMAKE_CXX_STANDARD 11)
else()
    set (CMAKE_CXX_STANDARD 14)
endif()

add_library(networkit "networkit/cpp/networkit.cpp")
set(NETWORKIT_LIBRARY "Networkit-Core-${NETWORKIT_BUILD_TYPE}")
set_target_properties(networkit PROPERTIES OUTPUT_NAME ${NETWORKIT_LIBRARY})

include_directories(.)
link_directories(${PROJECT_BINARY_DIR})

if (NETWORKIT_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})
    add_executable(networkit_tests "networkit/cpp/Unittests-X.cpp")
    target_link_libraries(networkit_tests
            ${GTEST_LIBRARIES}
            ${NETWORKIT_LIBRARY})
endif()

macro(networkit_add_source FILE)
        target_sources(networkit PUBLIC ${CMAKE_CURRENT_LIST_DIR}/${FILE}.cpp)
endmacro(networkit_add_source)

function(networkit_add_test FILE)
    if (NETWORKIT_BUILD_TESTS)
        target_sources(networkit_tests PUBLIC ${CMAKE_CURRENT_LIST_DIR}/${FILE}.cpp)
    endif()
endfunction(networkit_add_test)

add_subdirectory(networkit/cpp)
add_subdirectory(example)
