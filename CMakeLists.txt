cmake_minimum_required(VERSION 3.5)
project(networkit CXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
	message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# BUILD OPTIONS
option(NETWORKIT_BUILD_CORE "Build NetworKit core library" ON)
option(NETWORKIT_BUILD_TESTS "Build NetworKit C++ tests" OFF)
option(NETWORKIT_LOGGING "Build with logging support" ON)
option(NETWORKIT_STATIC "Build static libraries" OFF)
option(NETWORKIT_MONOLITH "Build single library (and tests is requested; required for shared lib)" ON)
option(NETWORKIT_NATIVE "Optimize for native architecture (often better performance)" OFF)
option(NETWORKIT_WARNINGS "Issue more warnings" OFF)
option(NETWORKIT_INCLUDESYMLINK "Create symlink to cpp directory" ON)
option(NETWORKIT_FLATINSTALL "Install into a flat directory structure (useful when building a Python package)" OFF)
set(NETWORKIT_PYTHON "" CACHE STRING "Directory containing Python.h. Implies MONOLITH=TRUE")
set(NETWORKIT_PYTHON_SOABI "" CACHE STRING "Platform specific file extension. Implies MONOLITH=TRUE")

set( NETWORKIT_CXX_STANDARD "11" CACHE STRING "CXX Version to compile with. Currently fixed to 11")

if (NETWORKIT_PYTHON)
	if (NOT NETWORKIT_MONOLITH)
		message(FATAL_ERROR "When building NetworKit as a Python module, NETWORKIT_MONOLITH=ON is required")
	endif()
	if(NOT NETWORKIT_PYTHON_SOABI)
		message(WARNING "No platform-specific file extension provided. Do not distribute library.")
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

if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
	("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
	("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
	if (NETWORKIT_NATIVE)
		set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} -march=native")
	endif()
	if (NETWORKIT_WARNINGS)
		set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} -Wall -Wextra -Wpedantic")
	endif()

elseif (MSVC)
	set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} /DNETWORKIT_OMP2 /MP /DNETWORKIT_WINDOWS")
	
	# TTMath requires running an ASM for MSVC which we disable here at the cost of worse performance
	set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} /DTTMATH_NOASM=1")
	 
else()
	message(WARNING "Support only GCC, Clang, MSVC and AppleClang. Your compiler may or may not work.")
endif()


# AppleClang support is added in CMake > 3.12
# the following passage checks if CMake > 3.12 is installed and manually
# links the omp lib if found
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
	if(${CMAKE_VERSION} VERSION_LESS 3.12)
 		message(WARNING "OpenMP linking with AppleClang was introduced in CMake 3.12. Current CMake Version ${CMAKE_VERSION}")
	else()
		find_library(OMP_FOUND NAMES omp)
		if(OMP_FOUND)
			set(NETWORKIT_LINK_FLAGS "-lomp ${NETWORKIT_LINK_FLAGS}")
		else()
			message(FATAL_ERROR "libomp was not found, but necessary to run NetworKit with AppleClang")
		endif()
	endif()
 endif()

# adding OpenMP flags
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(NETWORKIT_CXX_FLAGS "${NETWORKIT_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

	if (NOT MSVC)
		# In MSVC /openmp is only passed to the compiler, and must not be passed to the linker to avoid warnigns
		set(NETWORKIT_LINK_FLAGS "${NETWORKIT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
else()
	message(FATAL_ERROR "OpenMP not found")
endif()

if (CMAKE_SIZEOF_VOID_P LESS 8)
	if (MSVC)
		message(FATAL_ERROR "NetworKit supports only 64bit builds.
		                     Make sure to select a x64 target rather than x86, e.g. when invoking cmake with -G")
	else()
		message(FATAL_ERROR "NetworKit supports only 64bit builds")
	endif()
endif()

# specify linking flags for MacOS
if (APPLE)
	set(NETWORKIT_LINK_FLAGS "-undefined dynamic_lookup ${NETWORKIT_LINK_FLAGS}")
endif()

if(NOT NETWORKIT_BUILD_CORE)
	find_library(EXTERNAL_NETWORKIT_CORE NAMES networkit DOC "External NetworKit core library")
endif()

if(NETWORKIT_FLATINSTALL)
	set(NETWORKIT_LIB_DEST ".")
else()
	set(NETWORKIT_LIB_DEST "lib/")
endif()

################################################################################
# Use TLX as a CMake submodule
if(EXISTS "${PROJECT_SOURCE_DIR}/extlibs/tlx/CMakeLists.txt")
	add_subdirectory(extlibs/tlx)
else()
	message(FATAL_ERROR
			"Missing TLX library in extlibs/tlx "
			"Please run `git submodule update --init` to fetch the submodule.")
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
if(NETWORKIT_BUILD_CORE AND NETWORKIT_MONOLITH)
	if(NETWORKIT_STATIC)
		add_library(networkit networkit/cpp/networkit.cpp)
	else()
		add_library(networkit SHARED networkit/cpp/networkit.cpp)
	endif()

	set_target_properties(networkit PROPERTIES
			CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
			COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
			LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

	install(TARGETS networkit
			LIBRARY DESTINATION ${NETWORKIT_LIB_DEST}
			ARCHIVE DESTINATION ${NETWORKIT_LIB_DEST})

	target_link_libraries(networkit PRIVATE tlx)
endif()

if(NETWORKIT_PYTHON)
	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/networkit/_NetworKit.cpp")
		message(FATAL_ERROR "networkit/_NetworKit.cpp is missing. Invoke Cython manually.")
	endif()

	add_library(_NetworKit MODULE networkit/_NetworKit.cpp)
	target_include_directories(_NetworKit PRIVATE "${NETWORKIT_PYTHON}")

	if(NOT NETWORKIT_BUILD_CORE)
		if(NOT EXTERNAL_NETWORKIT_CORE)
			message(FATAL_ERROR "Core build is disabled but no external core library was found.")
		endif()
		target_link_libraries(_NetworKit PRIVATE ${EXTERNAL_NETWORKIT_CORE})
	else()
		target_link_libraries(_NetworKit PRIVATE networkit)
	endif()

	set_target_properties(_NetworKit PROPERTIES
				CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
				COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
				LINK_FLAGS "${NETWORKIT_LINK_FLAGS}"
				PREFIX ""
				OUTPUT_NAME "_NetworKit.${NETWORKIT_PYTHON_SOABI}")
	# DSOs on Apple OSes use different conventions for RPATH.
	if(APPLE)
		set_target_properties(_NetworKit PROPERTIES
				INSTALL_RPATH "@loader_path")
	else()
		set_target_properties(_NetworKit PROPERTIES
				INSTALL_RPATH "$ORIGIN")
	endif()
	install(TARGETS _NetworKit
			LIBRARY DESTINATION ${NETWORKIT_LIB_DEST})
endif()

# Register a new NetworKit module named ${modname}
# Files additionally passed are interpreted as PUBLIC source files to this module
function(networkit_add_module modname)
	if(NOT NETWORKIT_BUILD_CORE)
		return()
	endif()

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

		target_link_libraries(${MODULE_TARGET} PRIVATE tlx)

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

	if(NOT NETWORKIT_BUILD_CORE)
		return()
	endif()

	if(NOT NETWORKIT_MONOLITH)
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
	if(NOT NETWORKIT_BUILD_CORE)
		return()
	endif()

	if(NOT NETWORKIT_MONOLITH)
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
		if (MSVC)
			# While by default MSVC projects link against the shared runtime library
			# (and hence also NetworKit), GTest defaults to the static runtime lib.
			# Both must not be mix, so we request GTest here to also use the shared CRT.
			set( gtest_force_shared_crt ON CACHE BOOL "Always use msvcrt.dll" FORCE)
		endif()
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
				gtest
				networkit
				tlx
		)

		set_target_properties(networkit_tests PROPERTIES
			CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
			COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
			LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

		add_test(
			NAME networkit_tests
			COMMAND networkit_tests -t
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		)

		add_test(
			NAME networkit_tests_no_assertions
			COMMAND networkit_tests -r
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		)

	else()
		add_library(networkit_gtest_main STATIC networkit/cpp/Unittests-X.cpp)
		target_link_libraries(networkit_gtest_main
				PUBLIC gtest
				PRIVATE networkit_auxiliary)
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
					gtest networkit_gtest_main
				PRIVATE
					networkit_${MOD}
					tlx
					)
			set_target_properties(${TARGET_NAME} PROPERTIES
				CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
				COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
				LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

			foreach(dep IN LISTS ARGN)
				target_link_libraries(${TARGET_NAME} PRIVATE networkit_${dep})
			endforeach()

			if (${IS_TEST})
				add_dependencies(networkit_tests_${MOD} ${TARGET_NAME})

				add_test(
					NAME "${MOD}/${NAME}"
					COMMAND ${TARGET_NAME} -t
					WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
					)

				add_test(
					NAME "${MOD}/${NAME}_no_assertions"
					COMMAND ${TARGET_NAME} -r
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
# Benchmarks

# In case of monolithic builds we add the target networkit_benchmarks, and later add
# source files via networkit_add_gbenchmark. In case of non-monolithic builds, each
# networkit_add_test creates it own target.
if (NETWORKIT_BUILD_BENCHMARKS AND NETWORKIT_MONOLITH)
	add_executable(networkit_benchmarks networkit/cpp/networkit.cpp)

	target_link_libraries(networkit_benchmarks
			PRIVATE
			benchmark_main benchmark
			networkit
			)

	set_target_properties(networkit_tests PROPERTIES
			CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
			COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
			LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")
endif()

function(networkit_add_gbenchmark MOD NAME)
	if (NETWORKIT_BUILD_BENCHMARKS)
		set(BENCH_SOURCE ${CMAKE_CURRENT_LIST_DIR}/${NAME}.cpp)

		if (NETWORKIT_MONOLITH)
			target_sources(networkit_benchmarks PRIVATE ${BENCH_SOURCE})

		else()
			if (NOT TARGET networkit_${MOD})
				MESSAGE(FATAL_ERROR "Unknown NetworKit module '${MOD}'")
			endif()

			set(TARGET_NAME "networkit_${MOD}_${NAME}")

			add_executable(${TARGET_NAME} ${BENCH_SOURCE})
			message("${TARGET_NAME}: ${BENCH_SOURCE}")

			target_link_libraries(${TARGET_NAME}
				PRIVATE
					benchmark_main benchmark
					networkit_${MOD})

			set_target_properties(${TARGET_NAME} PROPERTIES
					CXX_STANDARD ${NETWORKIT_CXX_STANDARD}
					COMPILE_FLAGS "${NETWORKIT_CXX_FLAGS}"
					LINK_FLAGS "${NETWORKIT_LINK_FLAGS}")

			foreach(dep IN LISTS ARGN)
				target_link_libraries(${TARGET_NAME} PRIVATE networkit_${dep})
			endforeach()
		endif()
	endif()
endfunction(networkit_add_gbenchmark)

################################################################################
# Documentation

find_program(SPHINX_EXECUTABLE
	NAMES
		sphinx-build sphinx-build.exe
	PATHS
		/usr/bin
		/usr/local/bin
		/opt/local/bin
	DOC "Sphinx documentation generator")

if (NOT SPHINX_EXECUTABLE)
	message(STATUS "sphinx-build not found. Disable documentation targets")

else()
	message(STATUS "Found sphinx-build: ${SPHINX_EXECUTABLE}")
	add_custom_target(general_docs
			COMMAND rm -rf htmldocs
			COMMAND ${SPHINX_EXECUTABLE} ${CMAKE_SOURCE_DIR}/docs htmldocs)

	add_custom_target(cpp_api_docs DEPENDS general_docs
			COMMAND sed s:@BUILDDIR@:${CMAKE_CURRENT_BINARY_DIR}:
					${CMAKE_SOURCE_DIR}/docs/cpp_api/Doxyfile.in | doxygen -
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

	add_custom_target(docs DEPENDS general_docs cpp_api_docs)
endif()

################################################################################
# Subdirectories
add_subdirectory(networkit/cpp)
