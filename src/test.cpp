//HDK:
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_Format.h>
#include <CMD/CMD_Args.h>
#include <UT/UT_PtrArray.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_String.h>


int 
main(int argc, char *argv[])
{   
 
    if (argc == 1)
    {
        // usage();
        return 1;
    }
 
    /// Read argumets and options:
    CMD_Args args;
    std::string plane_str = "C";            // Default working plane.
    std::vector<std::string> planes_vector; // Planes vector in case user wants to filter > 1 plane
    planes_vector.push_back(plane_str);     // dito

    args.initialize(argc, argv);
    args.stripOptions("m:b:s:k:f:p:");
    // Options opt;

    return 0;
}