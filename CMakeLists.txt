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




# Specify the simulator source files to be compiled
# Note: using file(GLOB..) is not recommended in the CMake manual
# Temporarily: remove some of the unfinished source files from the compilation list

set(platoSimSourceFiles
#    ${CMAKE_CURRENT_SOURCE_DIR}/source/Camera.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/ConfigurationParameters.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Coordinates.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Detector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/FileUtilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/HDF5File.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/HDF5Writer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Logger.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/PointSpreadFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Simulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/StringUtilities.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/source/Telescope.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/source/TimeTicker.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/source/IOException.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/UnsupportedException.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/IllegalArgumentException.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/FileException.cpp
)

set(platoSimMain ${CMAKE_CURRENT_SOURCE_DIR}/source/main.cpp)




# Add all the platosim sources into a library
# Note: this prevents recompilation of the source files for the unit tests

add_library(platoSimLibrary ${platoSimSourceFiles})




# Specify the simulator unittest files to be compiled

set(unitTestSourceFiles
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/GTestConfigurationParameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/GTestFileUtilities.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/GTestPointSpreadFunction.h
)

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")


# Then those specific for the compiler for using C++11 support.

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++11 -Wno-deprecated-register -Wno-unused-variable")    
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()



# Create the platosim and testplatosim executables

add_executable(platosim ${platoSimMain})
add_executable(testplatosim ${unitTestMain} ${unitTestSourceFiles})



# Link the platosim and testplatosim executable with the required libraries

target_link_libraries(platosim platoSimLibrary hdf5 hdf5_cpp armadillo.6 yaml-cpp)
target_link_libraries(testplatosim platoSimLibrary hdf5 hdf5_cpp armadillo.6 yaml-cpp gtest gtest_main)



