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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/catch/include/)


# Specify the simulator source files to be compiled
              
file(GLOB sourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)


# Remove the platosim.cpp file from the list, as it unnecessary to build the library

list(REMOVE_ITEM sourceFiles platosim.cpp)


# Specify the simulator unittest files to be compiled

file(GLOB unitTestFiles ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)


# Set the compiler flags
# First those common to both gcc and clang:
# -O3   = level 3 optimization; 
# -Wall = enable all compiler's warning messages

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")


# Then those specific for the compiler for using C++11 support.

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++11 -Wno-deprecated-register")    
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()



# Create the platosim shared library

add_library(platosimlib SHARED ${sourceFiles})


# Create the platosim executable

add_executable(platosim ${CMAKE_CURRENT_SOURCE_DIR}/source/main.cpp)


# Create the testplatosim executable

add_executable(testplatosim ${unitTestFiles})


# Install the library in the lib/ folder

install(TARGETS platosimlib LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)


# Link the platosim executable with the PlatoSim library

target_link_libraries(platosim platosimlib)


# Link the testplatosim executable with the PlatoSim library

target_link_libraries(testplatosim platosimlib)



