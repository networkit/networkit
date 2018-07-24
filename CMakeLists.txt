cmake_minimum_required(VERSION 3.5)
project(networkit CXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
	message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# BUILD OPTIONS
option(NETWORKIT_BUILD_TESTS "Build with tests" OFF)
option(NETWORKIT_LOGGING "Build with logging support" ON)
option(NETWORKIT_MONOLITH "Build single library (and tests is requested; required for shared lib)" ON)
option(NETWORKIT_NATIVE "Optimize for native architecture (often better performance)" OFF)
option(NETWORKIT_WARNINGS "Issue more warnings" OFF)
option(NETWORKIT_INCLUDESYMLINK "Create symlink to cpp directory" ON)
set(NETWORKIT_PYTHON "" CACHE STRING "Directory containing Python.h. Implies MONOLITH=TRUE")
set(NETWORKIT_PYTHON_SOABI "" CACHE STRING "Platform specific file extension. Implies MONOLITH=TRUE")

set( NETWORKIT_CXX_STANDARD "11" CACHE STRING "CXX Version to compile with. Currently fixed to 11")

if (NETWORKIT_PYTHON)
	if (NOT NETWORKIT_MONOLITH)
		message(FATAL_ERROR "When building NetworKit as a Python module, NETWORKIT_MONOLITH=ON is required")
	endif()
	if(NOT NETWORKIT_PYTHON_SOABI)
		message(WARNING "No Platform specific file extension provided. Do not distribute library.")
	endif()
endif()

if (NOT CMAKE_BUILD_TYPE)
	message("Use Release Build Type as default")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################
# Compilation Flags
set(NETWORKIT_CXX_FLAGS "")
set(NETWORKIT_LINK_FLAGS "")

if (NOT NETWORKIT_LOGGING)
	set(NETWORKIT_CXX_FLAGS "-DNOLOGGING")
endif()

if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
	if (NETWORKIT_NATIVE)
		set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} -march=native")
	endif()
	if (NETWORKIT_WARNINGS)
		set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} -Wall -Wextra -Wpedantic")
	endif()
else()
	message(WARNING "Support only GCC and Clang; your compiler might work though anyways")
endif()

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	link_directories("${OPENMP_LIBRARIES}")

	set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(NETWORKIT_LINK_FLAGS "${OpenMP_CXX_FLAGS}")
else()
	message(FATAL_ERROR "OpenMP not found")
endif()

################################################################################
#include symlink
if(NETWORKIT_INCLUDESYMLINK)
	if(UNIX)
		add_custom_target(symlink_include ALL
				COMMAND ${CMAKE_COMMAND} -E make_directory "include"
				COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/networkit/cpp" "include/NetworKit")
	else()
		message(WARNING "Cannot create symlink on windows machine")
	endif()
endif()

# NETWORKIT MODULES
if (NETWORKIT_MONOLITH)
	add_library(networkit SHARED networkit/cpp/networkit.cpp)
	set_target_properties(networkit PROPERTIES
			CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
			COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
			LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")
	if(APPLE)
		set_target_properties(networkit PROPERTIES
				LINK_FLAGS "-undefined dynamic_lookup")
		target_link_libraries(networkit PRIVATE OpenMP::OpenMP_CXX)
	endif()
	if (NETWORKIT_PYTHON)
		if(EXISTS "${PROJECT_SOURCE_DIR}/networkit/_NetworKit.cpp")
			add_library(_NetworKit MODULE networkit/_NetworKit.cpp)
			target_link_libraries(_NetworKit PRIVATE networkit)
			target_include_directories(_NetworKit PRIVATE "${NETWORKIT_PYTHON}")
			set_target_properties(_NetworKit PROPERTIES
						CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
						COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
						LINK_FLAGS "${NETWORKIT_LINK_FLAGS}"
						PREFIX ""
						OUTPUT_NAME "_NetworKit.${NETWORKIT_PYTHON_SOABI}")
			if(APPLE)
				set_target_properties(_NetworKit PROPERTIES
						LINK_FLAGS "-undefined dynamic_lookup")
				target_link_libraries(_NetworKit PRIVATE OpenMP::OpenMP_CXX)
			endif()
		else()
			message(FATAL_ERROR "networkit/_NetworKit.cpp is missing. Invoke Cython manually.")
		endif()
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

		set_target_properties(${MODULE_TARGET} PROPERTIES
			CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
			COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
			LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

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
# Example: networkit_module_link_libraries(graph PRIVATE foobar_lib)
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

	if(EXISTS "${PROJECT_SOURCE_DIR}/extlibs/googletest/CMakeLists.txt")
		option(BUILD_GTEST "Builds the googletest subproject" ON)
		option(BUILD_GMOCK "Builds the googlemock subproject" OFF)
		add_subdirectory(extlibs/googletest)
	else()
		message(FATAL_ERROR
			"Missing GoogleTest and GoogleMock in extlibs/googletest. "
			"Please run `git submodule update --init` to fetch the submodule.")
	endif()

	if (NETWORKIT_MONOLITH)
		add_executable(networkit_tests networkit/cpp/Unittests-X.cpp)
		target_link_libraries(networkit_tests
			PRIVATE
				gtest_main gtest
				networkit
		)
		set_target_properties(networkit_tests PROPERTIES
			CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
			COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
			LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

		add_test(networkit_tests networkit_tests)
		if(UNIX)
			add_custom_target(symlink_input_folder ALL COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/input" "input")
			add_custom_target(symlink_output_folder ALL COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/output" "output")
		else()
			message(WARNING "Cannot create input and output symlinks on a windows machine")
		endif()
		endif()
endif()

# internal use only
# IS_TEST   indicates whether add_test should be invoked for executable
#		   an whether it should be assigned to the module's test target
# MOD	   Name of module the test/benchmark assigned to. It will
#		   join its namespace, assigned to its test target and linked
# TESTNAME  Name of the CPP File (excluding its .cpp extension). Will
#		   also be used to derive the test's name
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
					gtest_main gtest
				PRIVATE
					networkit_${MOD}
				)
			set_target_properties(${TARGET_NAME} PROPERTIES
				CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
				COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
				LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

			foreach(dep IN LISTS ARGN)
				target_link_libraries(${TARGET_NAME} PRIVATE networkit_${dep})
			endforeach()

			if (${IS_TEST})
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
