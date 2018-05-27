cmake_minimum_required(VERSION 3.9)
project(networkit CXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# BUILD OPTIONS
option(NETWORKIT_MONOLITH "Build single library (and tests is requested; required for shared lib)" ON)
option(NETWORKIT_LOGGING "Build with logging support" ON)
option(NETWORKIT_SANITIZE "Build with logging support" ON)

option(NETWORKIT_BUILD_TESTS "Build with tests (requires GTest)" OFF)

if (NOT CMAKE_BUILD_TYPE)
    message("Use Release Build Type as default")
    set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(OpenMP REQUIRED)

################################################################################
# MODULES

if (NETWORKIT_MONOLITH)
    add_library(networkit networkit/cpp/networkit.cpp)
    target_link_libraries(networkit PRIVATE OpenMP::OpenMP_CXX)
else()
    if (BUILD_SHARED_LIBS)
        message(WARNING
            "Currenlty only a monolithic shared library of NetworKit is supported.
             This is due to circular dependencies of modules.")
    endif()
endif()

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

    # Add source files (it's important to mark them private; otherwise
    # all targets linking to the lib, will recompile the objects from scratch)
    foreach(file ${ARGN})
        target_sources(${MODULE_TARGET}
            PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${file})
    endforeach()
endfunction()

# Analogous to target_link_libraries with KEYWORDS.
# Use the module's name (without networkit_ prefix) for target.
# In case of monolithic builds, the call is ignored.
# To link against another module use networkit_module_link_modules
# Example: networkit_module_link_libraries(graph PRIVATE OpenMP::OpenMP_CXX)
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

# Specifiy inter-module dependencies. The function expects a list of at least
# two module names (without the networkit_ prefix). The first one add all following
# ones as dependencies. In case of monolith build, the function does nothing.
# Example: networkit_module_link_modules(io graph) # io depends on graph
function(networkit_module_link_modules modname)
    if (NOT NETWORKIT_MONOLITH)
        foreach(dep IN LISTS ARGN)
            target_link_libraries(networkit_${modname} PUBLIC networkit_${dep})
        endforeach()
    endif()
endfunction()

################################################################################
# TESTING and BENCHMARKING
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

# internal use only
# IS_TEST   indicates whether add_test should be invoked for executable
#           an whether it should be assigned to the module's test target
# MOD       Name of module the test/benchmark assigned to. It will
#           join its namespace, assigned to its test target and linked
# TESTNAME  Name of the CPP File (excluding its .cpp extension). Will
#           also be used to derive the test's name
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

# Wrapper for networkit_add_extra with (IS_TEST=ON)
# Example: networkit_add_test(io SpecialIOGTest graph) compiles
#  io/test/SpecialIOGTest.cpp, registers is as an test of networkit_tests_io
#  and links it against io and graph.
function(networkit_add_test MOD NAME)
    networkit_add_extra(ON ${MOD} ${NAME} ${ARGN})
endfunction(networkit_add_test)

# Wrapper for networkit_add_extra with (IS_TEST=OFF)
function(networkit_add_benchmark MOD NAME)
    networkit_add_extra(OFF ${MOD} ${NAME} ${ARGN})
endfunction(networkit_add_benchmark)



################################################################################
# Subdirectories
add_subdirectory(networkit/cpp)
#add_subdirectory(example)
