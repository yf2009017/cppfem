# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (cppfem)

enable_testing()
set(MEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind)

# Set build type
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

#add_subdirectory (triangle)
add_subdirectory (src)
