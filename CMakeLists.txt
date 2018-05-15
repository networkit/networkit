cmake_minimum_required(VERSION 3.9)
project(networkit CXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# BUILD OPTIONS
option(NETWORKIT_OPENMP "Build with OpenMP support (requires OpenMP)" ON)
option(NETWORKIT_LOGGING "Build with logging support" ON)
option(NETWORKIT_SANITIZE "Build with logging support" ON)

option(NETWORKIT_BUILD_TESTS "Build with tests (requires GTest); implies NETWORKIT_BUILD_SHARED" OFF)

if (NOT CMAKE_BUILD_TYPE)
    message("Use Release Build Type as default")
    set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(OpenMP REQUIRED)

function(networkit_add_test MOD NAME)
    if (NETWORKIT_BUILD_TESTS)
        set(TARGET_NAME "networkit_${MOD}_${NAME}")

        add_executable(${TARGET_NAME}
                ${CMAKE_CURRENT_LIST_DIR}/${NAME}.cpp)

        target_link_libraries(${TARGET_NAME}
            PRIVATE
                GTest::GTest
                GTest::Main
                networkit_${MOD}
            )

        foreach(dep ${ARGN})
            target_link_libraries(${TARGET_NAME} PRIVATE networkit_${dep})
        endforeach()

        add_test(${TARGET_NAME} ${TARGET_NAME})
    endif()
endfunction(networkit_add_test)

if (NETWORKIT_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
endif()

add_subdirectory(networkit/cpp)
#add_subdirectory(example)
