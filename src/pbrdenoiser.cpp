//HDK:
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_Format.h>
#include <CMD/CMD_Args.h>
#include <UT/UT_PtrArray.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_String.h>

#ifdef BUILD_WITH_OIDN
#include <OpenImageDenoise/oidn.hpp>
#endif

#include <vector>
#include <iostream>

void usage(const char* program_name)
{
    std::cout << "USAGE: ";
    std::cout << program_name << " -i file-to-denoise.exr -o denoised-file.exr\n";
}

int 
main(int argc, char *argv[])
{   

    CMD_Args args;
    std::vector<std::string> planes_to_denoise{"C"}; 

    args.initialize(argc, argv);
    args.stripOptions("i:o:");

    const char * intput_filename;
    const char * output_filename;

    if (args.found('i') && args.found('o')) {
        intput_filename = args.argp('i');
        output_filename = args.argp('o');
    }
    else {
        std::cerr << "ERROR: We need image files to proceed.\n";
        usage(argv[0]);
        return 1;
    }

    IMG_FileParms input_parms = IMG_FileParms();
    input_parms.setDataType(IMG_FLOAT);    

    IMG_File *input_file       = IMG_File::open(intput_filename, &input_parms); 
    const IMG_Stat &input_stat = input_file->getStat();                      

    UT_PtrArray<PXL_Raster *> raster_array;
    PXL_Raster * raster = nullptr;        // we don't own it.
    const float * pixel_array = nullptr; // we don't own it.
    size_t xres = 0, yres = 0;

    const bool success = input_file->readImages(raster_array);
    if (!success) {
        std::cerr << "Integrity : Fail\n";
        return 1;
    }  else {
        // TODO: Multi raster denoise
        const std::string & plane_name = planes_to_denoise.at(0);
        const int raster_index  = input_stat.getPlaneIndex(plane_name.c_str());
        if (raster_index != -1) {
            raster = raster_array(raster_index);
            xres   = raster->getXres(); 
            yres   = raster->getYres();
            const void * pixels_ptr = raster->getPixels();
            pixel_array = static_cast<const float*>(pixels_ptr);
        }
        else
        {
            std::cerr << "No raster :" << plane_name << std::endl;
            return 1;
        }
    }  

    auto output_ptr = std::make_unique<float[]>(3*xres*yres);

   if(pixel_array && xres && yres) {
        // Create an Open Image Denoise device
        oidn::DeviceRef device = oidn::newDevice();
        device.commit();

        // Create a denoising filter
        oidn::FilterRef filter = device.newFilter("RT"); // generic ray tracing filter
        filter.setImage("color",  (void*)pixel_array,  oidn::Format::Float3, xres, yres);
        // filter.setImage("albedo", albedoPtr, oidn::Format::Float3, width, height); // optional
        // filter.setImage("normal", normalPtr, oidn::Format::Float3, width, height); // optional
        filter.setImage("output", output_ptr.get(), oidn::Format::Float3, xres, yres);
        filter.set("hdr", true); // image is HDR
        filter.commit();

        // Filter the image
        filter.execute();

        // Check for errors
        const char* errorMessage;
        if (device.getError(errorMessage) != oidn::Error::None) {
            std::cerr << "Error: " << errorMessage << std::endl;
            return 1; 
        }


         // OutputFile should be a copy of working frame with modified raster:
        IMG_Stat output_stat = IMG_Stat(input_stat);  
        output_stat.setFilename(output_filename);
        input_parms.setDataType(IMG_HALF); // Back to half when possible.
        IMG_File *output_file = IMG_File::create(output_filename, (const IMG_Stat)input_stat, &input_parms);
        raster->setRaster((void*)output_ptr.get(), false, false); // we own output_ptr so we don't give ownership to raster.

        // End
        if (output_file)
            output_file->writeImages(raster_array);
    }

    return 0;
}




