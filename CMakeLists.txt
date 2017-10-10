# inspired in http://mirkokiefer.com/blog/2013/03/cmake-by-example/

cmake_minimum_required(VERSION 2.8)

project(libmla_project)

# Semantic versioning
# http://semver.org/
set (libma_VERSION_MAJOR 0)
set (libma_VERSION_MINOR 1)
set (libma_VERSION_PATCH 0)


find_library(UMFPACK_LIBRARY NAMES umfpack )

INCLUDE_DIRECTORIES(
	${PROJECT_SOURCE_DIR}
)


# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
	"${PROJECT_SOURCE_DIR}/mla/config.h.in"
	"${PROJECT_BINARY_DIR}/mla/config.h"
)


add_definitions(-std=c++11 -Wall)


list(APPEND CMAKE_C_SOURCE_FILE_EXTENSIONS 
	c++
)

# set the compiler flags
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pedantic-errors -Wextra -Wfatal-errors")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")


enable_testing()

add_subdirectory(mla)
add_subdirectory(unit_tests)


# The installation targets
#install(TARGETS mla DESTINATION lib)

