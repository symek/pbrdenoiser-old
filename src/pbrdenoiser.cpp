//HDK:
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_Format.h>
#include <CMD/CMD_Args.h>
#include <UT/UT_PtrArray.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_String.h>

#include "filtering.hpp"

#ifdef BUILD_WITH_OIDN
#include <OpenImageDenoise/oidn.hpp>
#endif

#include <vector>
#include <iostream>



void usage(const char* program_name)
{
    std::cout << "USAGE: ";
    std::cout << program_name << " [options] -i file-to-denoise.exr -o denoised-file.exr\n";
    std::cout << "Options: \n";
    std::cout << "\t-i  name of the input file\n";
    std::cout << "\t-o  name of the output file\n";
    std::cout << "\t-d  <type> name of the denoiser   ('oidn' by default)\n";
    std::cout << "\t-l  image is in low dynamic range (hdri   by default)\n";
    std::cout << "\t-s  image is in sRGB color space  (linear by default)\n";
    std::cout << "\t-p  image planes to be filtered   ('C'    by default)\n";
    std::cout << "\t-P  number of denoising passes    (1      by default)\n";
}

int 
main(int argc, char *argv[])
{   

    CMD_Args args;
    std::vector<std::string> planes_to_denoise{"C"}; 

    args.initialize(argc, argv);
    args.stripOptions("d:p:i:o:lsP:f");

    const char * intput_filename;
    const char * output_filename;
    const char * albedo_name = "basecolor";
    const char * normals_name = "N";
    bool srgb = false;
    bool hdri = true;
    size_t passes = 1;
    size_t downscale = 0;

    if (args.found('i') && args.found('o')) {
        intput_filename = args.argp('i');
        output_filename = args.argp('o');
    }
    else {
        std::cerr << "ERROR: We need image files to proceed.\n";
        usage(argv[0]);
        return 1;
    }

    if (args.found('n')) {
        normals_name = args.argp('N');
    }

    if (args.found('a')) {
        albedo_name = args.argp('a');
    }

    if (args.found('l')) {
        hdri = false;
    }

    const char * dynamic = hdri ? "high" : "low";
    std::cout << "INFO:  Color range: " << dynamic << "\n"; 

    if (args.found('s')) {
        srgb = true;

    }
    const char * gamma = srgb ? "sRGB" : "linear";
    std::cout << "INFO:  Color gamma: " << gamma << "\n"; 

    if(args.found('P')) {   
        passes = args.iargp('P');
        std::cout << "INFO: Passes to run: " << passes << "\n"; 
    }

     if(args.found('f')) {   
        downscale = args.iargp('f');
        std::cout << "INFO: Downsample: x" << passes << "\n"; 
    }


    IMG_FileParms input_parms = IMG_FileParms();
    input_parms.setDataType(IMG_FLOAT);    
    input_parms.readAlphaAsPlane(); // so we don't have to stripe it away later on.

    IMG_File *input_file       = IMG_File::open(intput_filename, &input_parms); 
    const IMG_Stat &input_stat = input_file->getStat();                      

    PXL_Raster * beauty = nullptr;
    void * beauty_ptr; 
    PXL_Raster * normals = nullptr;
    void * normals_ptr;
    PXL_Raster * albedo = nullptr;
    void * albedo_ptr;
    size_t xres = 0, yres = 0;

    UT_PtrArray<PXL_Raster *> raster_array;
    const bool success = input_file->readImages(raster_array);
    if (!success) {
        std::cerr << "Integrity : Fail\n";
        return 1;
    } else {

        // TODO: Multi raster denoise
        const std::string & plane_name = planes_to_denoise.at(0);
        int raster_index  = input_stat.getPlaneIndex(plane_name.c_str());
        if (raster_index != -1) {
            beauty = raster_array(raster_index);
            xres   = beauty->getXres(); 
            yres   = beauty->getYres();
            beauty_ptr = beauty->getPixels();
            std::cout << "INFO: Target plane: " << plane_name << "\n"; 
        }
        else {
            std::cerr << "No plane to denoise found :" << plane_name << std::endl;
            return 1;
        }
    }  
    //
    int raster_index  = input_stat.getPlaneIndex(normals_name);
    if (raster_index != -1) {
        normals = raster_array(raster_index);
        normals_ptr = normals->getPixels();
            std::cout << "INFO: Normal plane: " << normals_name << "\n"; 

    }
    //
    raster_index  = input_stat.getPlaneIndex(albedo_name);
    if (raster_index != -1) {
        albedo = raster_array(raster_index);
        albedo_ptr = albedo->getPixels();
            std::cout << "INFO: Albedo plane: " << albedo_name << "\n"; 

    }


    auto output_ptr = std::make_unique<float[]>(3*xres*yres);

    // Filtering path:
    if(downscale && beauty && xres && yres)
    {
        const size_t size = xres*yres*3;
        std::vector<float> input;
        const float * beauty_f = (const float*) beauty_ptr;
        for(size_t i=0; i<size; ++i) {
            input[i] = beauty_f[i];
        }

        MitchellNetravali kernel;
        std::vector<float> result = downsample<MitchellNetravali>(input, downscale, kernel);

        IMG_FileParms parms;
        parms.scaleByFactor(1.0f/downscale, 1.0f/downscale);
        if (!IMG_File::copyToFile(intput_filename, output_filename, parms)) {
            std::cerr << "Can't create " << output_filename << "\n";
            return 1;
        }

        IMG_FileParms downscaled_parms = IMG_FileParms();
        downscaled_parms.setDataType(IMG_FLOAT);    
        downscaled_parms.readAlphaAsPlane(); // so we don't have to stripe it away later on.
        IMG_File *downscaled = IMG_File::open(output_filename, &downscaled_parms); 
        const IMG_Stat &downscaled_stat = downscaled->getStat();   

        IMG_Stat output_stat = IMG_Stat(downscaled_stat);  
        output_stat.setFilename(output_filename);
        downscaled_parms.setDataType(IMG_HALF); // Back to half when possible.
        IMG_File *output_file = IMG_File::create(output_filename, 
            (const IMG_Stat)input_stat, &downscaled_parms);

        UT_PtrArray<PXL_Raster *> scaled_raster_array;
        PXL_Raster * scale_raster = nullptr;
        void * scaled_ptr; 
        const bool success = output_file->readImages(scaled_raster_array);
        if (!success) {
            std::cerr << "Integrity : Fail\n";
            return 1;
        } else {

            // TODO: Multi raster denoise
            const std::string & plane_name = planes_to_denoise.at(0);
            int raster_index  = input_stat.getPlaneIndex(plane_name.c_str());
            if (raster_index != -1) {
                scale_raster = scaled_raster_array(raster_index);
                xres   = scale_raster->getXres(); 
                yres   = scale_raster->getYres();
                scaled_ptr = scale_raster->getPixels();
                std::cout << "INFO: Target plane: " << plane_name << "\n"; 
            }
            else {
                std::cerr << "No plane to denoise found :" << plane_name << std::endl;
                return 1;
            }
        }  
        
        // this does
        for(size_t y=0; y<yres; ++y) {
            // const float * output_f = (const float *)out_ptr;
            const float * output_f = const_cast<float *>(&(result.front()));
            const void * data = (const void*)&output_f[y*xres*3];
            scale_raster->writeToRow(y, data);
        }
       
        // End
        if (output_file) {
            output_file->writeImages(scaled_raster_array);
            output_file->close();
        }








        return 0;   
    }


    // denoising path:
   if(beauty && xres && yres) {
        // Create an Open Image Denoise device
        oidn::DeviceRef device = oidn::newDevice();
        device.commit();

        // Create a denoising filter
        oidn::FilterRef filter = device.newFilter("RT"); // generic ray tracing filter
        if (albedo) {
            filter.setImage("albedo", albedo_ptr, oidn::Format::Float3, xres, yres); // optional
        }
        if (normals) {
            filter.setImage("normal", normals_ptr, oidn::Format::Float3, xres, yres); // optional
        }

        std::cout << "INFO: Start denoising...\n" << std::flush; 

        //we swap those to run denoiser multiply times
        void * in_ptr  = beauty_ptr;
        void * out_ptr = output_ptr.get();

        for(int pass=0; pass<passes; ++pass)
        {
            filter.setImage("color",  in_ptr,  oidn::Format::Float3, xres, yres);
            filter.setImage("output", out_ptr, oidn::Format::Float3, xres, yres);
            filter.set("hdr", hdri); // image is HDR
            filter.set("srgb", srgb); // image is srgb
            filter.commit();
            // Filter the image
            std::cout << "INFO: Pass " << pass+1<< "\n" << std::flush;
            filter.execute();

            if (pass < passes-1) {
                void *tmp = out_ptr; 
                out_ptr = in_ptr;
                in_ptr = tmp;
            }
        }


        std::cout << "INFO: done.\n"; 


        // Check for errors
        const char* errorMessage;
        if (device.getError(errorMessage) != oidn::Error::None) {
            std::cerr << "Error: " << errorMessage << std::endl;
            return 1; 
        }


        std::cout << "INFO: Saving images to: " << output_filename << "\n"; 

         // output_file should be a copy of working frame with modified raster:
        IMG_Stat output_stat = IMG_Stat(input_stat);  
        output_stat.setFilename(output_filename);
        input_parms.setDataType(IMG_HALF); // Back to half when possible.
        IMG_File *output_file = IMG_File::create(output_filename, 
            (const IMG_Stat)input_stat, &input_parms);
        // This does't work :(
        // raster->setRaster((void*)output_ptr.get(), true, true); 
        
        // this does
        for(size_t y=0; y<yres; ++y) {
            const float * output_f = (const float *)out_ptr;
            const void * data = (const void*)&output_f[y*xres*3];
            beauty->writeToRow(y, data);
        }
       
        // End
        if (output_file) {
            output_file->writeImages(raster_array);
            output_file->close();
        }
    }


    return 0;
}




