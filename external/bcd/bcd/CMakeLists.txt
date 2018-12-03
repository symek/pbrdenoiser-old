# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( pbrdenoiser )


list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )
find_package (Threads)

# link_directories($ENV{REZ_BOOST_ROOT}/lib)
include_directories($ENV{REZ_BOOST_ROOT}/include)
include_directories($ENV{REZ_OPENEXR_ROOT}/include/OpenEXR)
include_directories($ENV{REZ_ILMBASE_ROOT}/include/OpenEXR)
include_directories("../ext/eigen")


# Add a executable.
set( executable_name "bcd")

file(GLOB sources "*.h" "*.cpp")

set(CMAKE_CXX_STANDARD 14)

# Can't link to hboost with standalone apps without that. 
link_directories($ENV{HFS}/dsolib)

add_executable( ${executable_name}
    ${sources}
)

target_link_libraries (${executable_name} ${CMAKE_THREAD_LIBS_INIT})
# target_compile_features(${executable_name} PRIVATE cxx_std_14)
# add_library( boost_program_options STATIC IMPORTED )
# set_target_properties(boost_program_options PROPERTIES IMPORTED_LOCATION $ENV{REZ_BOOST_ROOT}/lib/boost_program_options.a)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${executable_name} Houdini )
target_link_libraries( ${executable_name} IlmImf )
target_link_libraries( ${executable_name} tbb )
# target_link_libraries( ${executable_name} hboost_program_options )
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${executable_name} )
SET( CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-ggdb" )
