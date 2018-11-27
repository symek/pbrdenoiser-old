# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( pbrdenoiser )

# include(ExternalProject)

# ExternalProject_Add( bcd
#   SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/bcd"
#   CMAKE_ARGS -DBCD_USE_CUDA=OFF  -DBCD_BUILD_GUI=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/lib
#  )


list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )
find_package (Threads)

# Eigen time
# set( EIGEN_INCLUDE_DIR "$ENV{EIGEN_INCLUDE_DIR}")
# if (NOT EXISTS ${EIGEN_INCLUDE_DIR})
#     message(FATAL_ERROR "Specify EIGEN location with $EIGEN_INCLUDE_DIR. CMake will exit now.")
# else()
#     message(STATUS "EIGEN_INCLUDE_DIR: ${EIGEN_INCLUDE_DIR}")
# endif()

# include_directories(${EIGEN_INCLUDE_DIR})
# # external root path 
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Can't link to hboost with standalone apps without that. 
link_directories($ENV{HFS}/dsolib)
include_directories($ENV{REZ_BOOST_ROOT}/include)
include_directories($ENV{REZ_OPENEXR_ROOT}/include/OpenEXR)
include_directories($ENV{REZ_ILMBASE_ROOT}/include/OpenEXR)
include_directories("../ext/eigen")

# Add a executable.
set( executable_name "bcd")

file(GLOB sources "*.h" "*.cpp")

set(CMAKE_CXX_STANDARD 14)
add_executable( ${executable_name}
    # external/cxxopts/include/cxxopts.hpp
    # $ENV{HT}/include/hboost/program_options.hpp 
    # src/pbrdenoiser.cpp
    # test.cpp
    ${sources}
)
target_link_libraries (${executable_name} ${CMAKE_THREAD_LIBS_INIT})
# target_compile_features(${executable_name} PRIVATE cxx_std_14)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${executable_name} Houdini )
target_link_libraries(${executable_name} IlmImf )
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${executable_name} )
