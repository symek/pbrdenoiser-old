// HDK:
#include <CMD/CMD_Args.h>
#include <IMG/IMG_File.h>
#include <IMG/IMG_Format.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_PtrArray.h>

#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "pbrdenoiser.hpp"

void usage(const char *program_name) {
    std::cout << "USAGE: ";
    std::cout << program_name
              << " [options] -i file-to-denoise.exr -o denoised-file.exr\n";
    std::cout << "Options: \n";
    std::cout << "\t-i  name of the input file\n";
    std::cout << "\t-o  name of the output file\n";
    std::cout << "\t-d  <type> name of the denoiser   ('oidn' by default)\n";
    std::cout << "\t-l  image is in low dynamic range (hdri   by default)\n";
    std::cout << "\t-s  image is in sRGB color space  (linear by default)\n";
    std::cout << "\t-p  image planes to be filtered   ('C'    by default)\n";
    std::cout << "\t-P  number of denoising passes    (1      by default)\n";
    std::cout << "\t-f  Downscaling factor            (1      by default)\n";
}

int main(int argc, char *argv[]) {

    CMD_Args args;
    pbrd::monotonic_timer timer;

    args.initialize(argc, argv);
    args.stripOptions("d:p:i:o:lsP:f:");

    using Info = pbrd::monotonic_timer::Info;
    // this might be extended by arguments:
    std::vector<std::string> planes_to_denoise{"C"};

    const char *input_filename;
    const char *output_filename;
    const char *albedo_name = "basecolor";
    const char *normal_name = "N";

    bool srgb = false;
    bool hdri = true;
    size_t passes = 1;
    size_t downscale = 0;
    size_t xres, yres;
    xres = yres = 0;

    if (args.found('i') && args.found('o')) {
        input_filename = args.argp('i');
        output_filename = args.argp('o');
    } else {
        std::cerr << "ERROR: We need image files to proceed.\n";
        usage(argv[0]);
        return 1;
    }

    if (args.found('n')) {
        normal_name = args.argp('N');
    }

    if (args.found('a')) {
        albedo_name = args.argp('a');
    }

    if (args.found('l')) {
        hdri = false;
    }

    const char *dynamic = hdri ? "high" : "low";
    std::cout << "INFO: Color  range: " << dynamic << "\n";

    if (args.found('s')) {
        srgb = true;
    }
    const char *gamma = srgb ? "sRGB" : "linear";
    std::cout << "INFO: Color  gamma: " << gamma << "\n";

    if (args.found('P')) {
        passes = args.iargp('P');
        std::cout << "INFO: Passes to run: " << passes << "\n";
    }

    if (args.found('f')) {
        downscale = args.iargp('f');
        std::cout << "INFO: Downsampling: x " << downscale << "\n";
    }

    if (!static_cast<bool>(std::ifstream(input_filename))) {
        std::cerr << "ERROR: Can't open input file " << input_filename
                  << " (It probably doesn't exist)\n";
        return 1;
    }

    timer.begin("Reading file... ");

    // Readin rasters from file:
    UT_PtrArray<PXL_Raster *> raster_array;
    auto input_file = pbrd::read_rasters_as_float(input_filename, raster_array);

    if (!input_file) {
        std::cerr << "Integrity : Fail\n";
        return 1;
    }

    xres = input_file->getStat().getXres();
    yres = input_file->getStat().getYres();

    // beauty
    const char *plane_name = planes_to_denoise.at(0).c_str();
    PXL_Raster *beauty =
        pbrd::get_raster_by_name(plane_name, input_file.get(), raster_array);

    if (!beauty) {
        std::cerr << "No plane to denoise found:" << plane_name << std::endl;
        input_file->close();
        return 1;
    }

    timer.end();
    std::cout << "INFO: Noisy  plane: " << plane_name << " (found) \n";

    void *denoise_input_buffer = beauty->getPixels();

    // normals, albedo
    PXL_Raster *normal =
        pbrd::get_raster_by_name(normal_name, input_file.get(), raster_array);
    PXL_Raster *albedo =
        pbrd::get_raster_by_name(albedo_name, input_file.get(), raster_array);
    if (normal) {
        std::cout << "INFO: Normal plane: " << normal_name << "\n";
    }
    if (albedo) {
        std::cout << "INFO: Albedo plane: " << albedo_name << "\n";
    }

#if 0
    // we need this not to evaporate too soon:
    std::vector<float> scaled_result;

    // Filtering path (will replace all above )
    if(downscale && beauty && xres && yres) {

        // Copy current beauty
        const size_t size = xres*yres*3;
        const size_t fullxres = xres;
        const size_t fullyres = yres;
        std::vector<float> full_res_input(size);
        const float * beauty_f = (const float*)beauty->getPixels();

        for(size_t i=0; i<size; ++i) {
            full_res_input[i] = beauty_f[i];
        }
        // shortcut to scale down all rasters:
        timer.begin("Down scaling aov...");
        IMG_FileParms parms;
        parms.scaleImageBy(1.0f/downscale, 1.0f/downscale);
        if (!IMG_File::copyToFile(input_filename, output_filename, &parms)) {
            std::cerr << "Can't create " << output_filename << "\n";
            return 1;
        }
        timer.end();

        // Reopen scaled image
        // we overwrite here input_file unique_ptr, raster_array, resolution x/y,
        // all rasters' pointer and denoise_input_buffer
        // Note: we assume correctness, because we've just created that file ourself.
        raster_array.clear();
        input_file = pbrd::read_rasters_as_float(output_filename, raster_array);
        beauty = pbrd::get_raster_by_name(plane_name,  input_file.get(), raster_array);
        normal = pbrd::get_raster_by_name(normal_name, input_file.get(), raster_array);
        albedo = pbrd::get_raster_by_name(albedo_name, input_file.get(), raster_array);
        xres   = input_file->getStat().getXres();
        yres   = input_file->getStat().getYres();

        // Filter an image
        gs::MitchellNetravali kernel;
        timer.begin(Info{"Down scaling ", planes_to_denoise.at(0),
            " with ",  kernel.name(),  "..." });

        auto result = gs::downsample<gs::MitchellNetravali>(full_res_input, downscale, kernel);
        scaled_result.swap(result);

        timer.end();

        denoise_input_buffer = (void*)const_cast<float*>(&scaled_result.front());

    }
#endif

    // denoising:
    if (beauty && xres && yres) {

        //
        timer.begin(Info{"Start denoising ", planes_to_denoise.at(0), "... "});

        const char *error;
        auto output = pbrd::intel::filter_with_oidn(
            denoise_input_buffer, albedo->getPixels(), normal->getPixels(),
            &error, xres, yres, hdri, srgb);

        if (!output) {
            std::cerr << "OIDN :" << error << "\n";
            return 1;
        }

        timer.end();

        // copy back to
        for (size_t y = 0; y < yres; ++y) {
            // const float * data = (const float * )denoise_input_buffer;
            const float *scanline = &output.get()[y * xres * 3];
            beauty->writeToRow(y, (const void *)scanline);
        }

        timer.begin(Info{"Saving images to ", output_filename, "... "});

        const IMG_Stat &stat = input_file->getStat();
        if (!pbrd::save_rasters_to_file(output_filename, stat, raster_array)) {
            std::cerr << "ERROR: Can't save image: " << output_filename << "\n";
            return 1;
        }

        timer.end();
    }

    return 0;
}
