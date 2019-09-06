# Project specification
project(biggles CXX C)
cmake_minimum_required(VERSION 2.8)
enable_testing()

option(BUILD_LATEX "Build LaTeX documents in doc/ directory." OFF)

# this is a bit of a hack:
# it ensures that the code it optimised in RelWithDebInfo mode but NDEBUG is not set
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message("Build type is ${CMAKE_BUILD_TYPE}")

# Enable use of pkg-config for finding modules
find_package(PkgConfig)

# Unless otherwise specified, use the RelWithDebInfo build type
if((NOT CMAKE_BUILD_TYPE) OR ("${CMAKE_BUILD_TYPE}" STREQUAL "None"))
    message(STATUS "Since CMAKE_BUILD_TYPE was not set, use RelWithDebInfo as a default.")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Uncomment to enable all warnings and position independent code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC")

# the eigen library has an memory alignment method that can cause problems. Switch it off
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_ALIGN_STATICALLY")

# Module locations
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

# Record this directory as the root of the Biggles source code
set(BIGGLES_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Use lsb-release to find out what we're building on if supported
include(CheckLSBTypes)

# git commit number
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --tags --always
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "BIGGLES_COMMIT_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "Git version: ${BIGGLES_COMMIT_VERSION}" )
  ELSE(GIT_FOUND)
    SET(BIGGLES_COMMIT_VERSION 0)
  ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

#string(TIMESTAMP "BIGGLES_BUILD_DATE" "%d.%m.%y" UTC)
MESSAGE( STATUS "build date: ${BIGGLES_BUILD_DATE}" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBIGGLES_COMMIT_VERSION=\\\"${BIGGLES_COMMIT_VERSION}\\\"")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBIGGLES_BUILD_DATE=\\\"${BIGGLES_BUILD_DATE}\\\"")

# Make use of CPack to generate packages
set(CPACK_GENERATOR "TGZ")
if(NOT "${LSB_DISTRIBUTOR_ID}" STREQUAL "unknown")
    set(CPACK_SYSTEM_NAME "${LSB_DISTRIBUTOR_ID}-${LSB_RELEASE}-${CMAKE_SYSTEM_PROCESSOR}")
endif(NOT "${LSB_DISTRIBUTOR_ID}" STREQUAL "unknown")
set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_SYSTEM_NAME}")
include(CPack)

# Make sure the compiler has OpenMP
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Try to find the Python interpreter. We need version 2 for matplotlib.
set(PythonInterp_FIND_VERSION "2.7")
set(PythonLibs_FIND_VERSION "2.7")
find_package(PythonInterp)
find_package(PythonLibs)

# Required boost libraries
find_package(Boost REQUIRED COMPONENTS filesystem program_options thread system)
find_package(Boost COMPONENTS python)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

# Optional boost libraries
set(BUILD_PYTHON_BINDINGS NO CACHED)
if(Boost_PYTHON_LIBRARY AND PYTHONLIBS_FOUND)
    set(BUILD_PYTHON_BINDINGS YES CACHED)
else(Boost_PYTHON_LIBRARY AND PYTHONLIBS_FOUND)
    message(STATUS "Boost python library was not found; Python bindings will not be built.")
endif(Boost_PYTHON_LIBRARY AND PYTHONLIBS_FOUND)

# Find the eigen linear algebra package
pkg_check_modules(EIGEN eigen3 REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

# Find the Jansson JSON library
pkg_check_modules(JANSSON jansson REQUIRED)
include_directories(${JANSSON_INCLUDE_DIRS})
link_directories(${JANSSON_LIBRARY_DIRS})

# Add the include directory to the include path
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# A 'doc' target that documentation generators can hang off
add_custom_target(doc)

# Common utility library
add_subdirectory(utility)

# Biggles library
add_subdirectory(biggles)

# Python bindings and tools
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif(BUILD_PYTHON_BINDINGS)

# Testing
add_subdirectory(test)

# Documentation
add_subdirectory(doc)

## Test data
#add_subdirectory(test-data)

# Install header files
install(DIRECTORY include/ DESTINATION include)
