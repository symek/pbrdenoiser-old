# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( pbrdenoiser )

include(ExternalProject)

option(WITH_BCD "Build with bcd denoiser" OFF)
option(WITH_OIDN "Build with Intel Open Image Denoiser" OFF)

if (WITH_BCD)
	set(BUILD_WITH_BCD TRUE)
	ExternalProject_Add( bcd
	  SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/bcd"
	  CMAKE_ARGS -DBCD_USE_CUDA=OFF  -DBCD_BUILD_GUI=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/lib
	 )
endif()

if (WITH_OIDN)
	list( APPEND CMAKE_PREFIX_PATH "$ENV{OIDN}/lib/cmake")
	find_package( OpenImageDenoise REQUIRED )
	set(BUILD_WITH_OIDN TRUE)
	add_definitions(-DBUILD_WITH_OIDN)
endif()

list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )


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
link_directories($ENV{HDSO})

# Add a executable.
set( executable_name pbrdenoiser)



add_executable( ${executable_name}
    # external/cxxopts/include/cxxopts.hpp
    # $ENV{HT}/include/hboost/program_options.hpp 
    src/pbrdenoiser.cpp
    # src/test.cpp
)


# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${executable_name} Houdini )
if (BUILD_WITH_OIDN)
	target_link_libraries( ${executable_name} OpenImageDenoise )
endif()
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${executable_name} )
