# cmake file to build PlatoSim3
# 
# Usage: ($ mkdir build)
#        $ cd build
#        $ cmake ..  (alternatively: cmake -D CMAKE_CXX_COMPILER=/usr/bin/clang++ .. on mac OSX)
#        ($ make clean)
#        $ make -j 4 
#        $ make install
#
# After a successful build you can safely remove all files and directories from build/
#
# Notes:
#   - Prefer 'cmake ..' in the build/ directory above 'cmake' in the
#     parent directory, because the former does an out-of-source build.
#


project(PlatoSim)


# cmake version shouldn't be too old

cmake_minimum_required(VERSION 2.8)



# Specify the locations of the header files

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/hdf5-1.8.16/include
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/armadillo-6.400.3/include
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/yaml-cpp/include
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/googletest/include
)


# Specify the simulator source files to be compiled, but exclude main.cpp.
# Temporarily: remove some of the unfinished source files from the compilation list.

file(GLOB platosimSourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
list(REMOVE_ITEM platosimSourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/main.cpp)

list(REMOVE_ITEM platosimSourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/camera.cpp 
                                     ${CMAKE_CURRENT_SOURCE_DIR}/source/telescope.cpp
                                     ${CMAKE_CURRENT_SOURCE_DIR}/source/timeticker.cpp)


# Define the main.cpp of platosim separately.

set(platosimMain ${CMAKE_CURRENT_SOURCE_DIR}/source/main.cpp)


# Add a platosim library. This way we need to compile the platosim source files
# only once, when using them for platosim as well as testplatosim.

add_library(platosimLibrary ${platosimSourceFiles})


# Specify the simulator unittest files to be compiled, but exclude again the main().

file(GLOB unitTestSourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
list(REMOVE_ITEM unitTestSourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/AllTests.cpp)


# Define the file with main() function for the unit tests separately. 

set(unitTestMain ${CMAKE_CURRENT_SOURCE_DIR}/tests/AllTests.cpp)


# Specify the location of the dependency libraries

link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/hdf5-1.8.16/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/armadillo-6.400.3/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/yaml-cpp/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Installs/googletest/lib
)


# Set the compiler flags
# First those common to both gcc and clang:
# -O3   = level 3 optimization; 
# -Wall = enable all compiler's warning messages

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wno-deprecated-register -Wno-unused-variable")


# Then those specific for the compiler for using C++11 support.

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++11")    
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()


# Create the platosim and testplatosim executables

add_executable(platosim ${platosimMain})
add_executable(testplatosim ${unitTestMain} ${unitTestSourceFiles})


# Link the platosim and testplatosim executable with the required libraries

target_link_libraries(platosim platosimLibrary hdf5 hdf5_cpp armadillo.6 yaml-cpp)
target_link_libraries(testplatosim platosimLibrary hdf5 hdf5_cpp armadillo.6 yaml-cpp gtest gtest_main)



