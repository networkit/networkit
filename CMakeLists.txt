cmake_minimum_required(VERSION 3.9)
project(networkit CXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# BUILD OPTIONS
option(NETWORKIT_OPENMP "Build with OpenMP support (requires OpenMP)" ON)
option(NETWORKIT_LOGGING "Build with logging support" ON)
option(NETWORKIT_SANITIZE "Build with logging support" ON)
option(NETWORKIT_MONOLITH "Build single library (and test)" ON)

option(NETWORKIT_BUILD_TESTS "Build with tests (requires GTest); implies NETWORKIT_BUILD_SHARED" OFF)

if (NOT CMAKE_BUILD_TYPE)
    message("Use Release Build Type as default")
    set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(OpenMP REQUIRED)

if (NETWORKIT_MONOLITH)
    add_library(networkit networkit/cpp/networkit.cpp)
    target_link_libraries(networkit PRIVATE OpenMP::OpenMP_CXX)
endif()

################################################################################
# MODULES

# Register a new NetworKit module named ${modname}
# Files additionally passed are interpreted as PUBLIC source files to this module
function(networkit_add_module modname)
    if(NETWORKIT_MONOLITH)
        # in case we are building a monolith, no submodule are registered
        # and we simple add the source file to the networkkit target
        set(MODULE_TARGET "networkit")

    else()
        set(MODULE_TARGET "networkit_${modname}")

        add_library(${MODULE_TARGET}
                    ${PROJECT_SOURCE_DIR}/networkit/cpp/networkit.cpp)


        # All tests added to this module will will also become a dependency
        # of networkit_tests_MODNAME. This target hence allows to build all
        # tests associated with this module
        if (NETWORKIT_BUILD_TESTS)
            add_custom_target(networkit_tests_${modname})
        endif()
    endif()

    foreach(file ${ARGN})
        target_sources(${MODULE_TARGET} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${file})
    endforeach()
endfunction()

function(networkit_module_link_libraries modname)
    set(options )
    set(oneValueArgs )
    set(multiValueArgs PRIVATE PUBLIC)
    cmake_parse_arguments(NMLL
            "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT NETWORKIT_MONOLITH)
        target_link_libraries(networkit_${modname}
                PRIVATE ${NMLL_PRIVATE}
                PUBLIC ${NMLL_PUBLIC})
    endif()
endfunction()

function(networkit_module_link_modules modname)
    if (NOT NETWORKIT_MONOLITH)
        foreach(dep IN LISTS ARGN)
            target_link_libraries(networkit_${modname} PUBLIC networkit_${dep})
        endforeach()
    endif()
endfunction()

################################################################################
# TESTING
if (NETWORKIT_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)

    if (NETWORKIT_MONOLITH)
        add_executable(networkit_tests networkit/cpp/networkit.cpp)
        target_link_libraries(networkit_tests
            PRIVATE
                GTest::GTest
                GTest::Main
                OpenMP::OpenMP_CXX
                networkit
        )

        add_test(networkit_tests networkit_tests)
    endif()
endif()

function(networkit_add_extra IS_TEST MOD NAME)
    if (NETWORKIT_BUILD_TESTS)
        set(TEST_SOURCE ${CMAKE_CURRENT_LIST_DIR}/${NAME}.cpp)

        if (NETWORKIT_MONOLITH)
            target_sources(networkit_tests PRIVATE ${TEST_SOURCE})

        else()
            if (NOT TARGET networkit_${MOD})
                MESSAGE(FATAL_ERROR "Unknown NetworKit module '${MOD}'")
            endif()

            set(TARGET_NAME "networkit_${MOD}_${NAME}")

            add_executable(${TARGET_NAME} ${TEST_SOURCE})

            target_link_libraries(${TARGET_NAME}
                PRIVATE
                    GTest::GTest
                    GTest::Main
                    OpenMP::OpenMP_CXX
                PRIVATE
                    networkit_${MOD}
                )

            foreach(dep IN LISTS ARGN)
                target_link_libraries(${TARGET_NAME} PRIVATE networkit_${dep})
            endforeach()

            if (${IS_TEST})
                message(${TARGET_NAME})
                add_dependencies(networkit_tests_${MOD} ${TARGET_NAME} )
                add_test(
                    NAME "${MOD}/${NAME}"
                    COMMAND ${TARGET_NAME}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    )
            endif()
        endif()
    endif()
endfunction()

function(networkit_add_test MOD NAME)
    networkit_add_extra(ON ${MOD} ${NAME} ${ARGN})
endfunction(networkit_add_test)

function(networkit_add_benchmark MOD NAME)
    networkit_add_extra(OFF ${MOD} ${NAME} ${ARGN})
endfunction(networkit_add_benchmark)




add_subdirectory(networkit/cpp)
#add_subdirectory(example)
