// // bcd
#include "SamplesAccumulator.h"
#include "Utils.h"
#include "DeepImage.h"
#include "ImageIO.h"
#include "Denoiser.h"
#include "MultiscaleDenoiser.h"
#include "IDenoiser.h"

// HDK
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_FileTypes.h>

#include <UT/UT_Vector2.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_ParallelUtil.h>
#include <UT/UT_Options.h>
#include <CMD/CMD_Args.h>
// #include <hboost/program_options.hpp>

//own
#include "bcd.h"
#include "stopwatch.h"
// std
#include <iostream>
#include <memory>
#include <string>
#include <array>


// namespace po = hboost::program_options;

int main(const int argc, const char** argv)
{
    Stopwatch<std::chrono::seconds> watch;

    CMD_Args args;
    args.initialize(argc, argv);
    args.stripOptions("s:d:b:o:a:S:");

    if (argc < 4) {
        std::cerr << "pbrdenoiser -s subpixels.exr -b beauty.exr -o denoised.exr\n";
        return 1;
    }

    //

    if (!args.found('s') || !args.found('b') || !args.found('o')) {
        std::cerr << "pbrdenoiser -s subpixels.exr -b beauty.exr -o denoised.exr\n";
        return 1;
    }

    const auto subpixel_filename1 = std::string(args.argp('s')); 
    const auto beauty_filename    = std::string(args.argp('b')); 
    const auto output_filename    = std::string(args.argp('o')); 

    uint bcd_scales = 4;
    if(args.found('S')) {
        bcd_scales = args.iargp('S');
    }


    std::string subpixel_filename2;
    if (args.found('d')) { 
        subpixel_filename2 = std::string(args.argp('d')); 
    }

    //
    auto parms = std::make_unique<IMG_FileParms>();
    parms->readAlphaAsPlane();
    parms->setDataType(IMG_FLOAT);

    //
    std::unique_ptr<IMG_File> subpixel1(nullptr);
    subpixel1.reset(IMG_File::open(subpixel_filename1.c_str(), parms.get()));
    
    if(!subpixel1) {
        std::cerr << "Can't open " << subpixel_filename1 << '\n';
        return 1;
    }
    const uint sub1_width = subpixel1->getStat().getXres();
    const uint sub1_height= subpixel1->getStat().getYres();

    std::unique_ptr<IMG_File> beauty(nullptr);
    beauty.reset(IMG_File::open(beauty_filename.c_str(), parms.get()));
    if (!beauty) {
        std::cerr << "Can't open beauty file " << beauty_filename << '\n';
        return 1;
    } 
    //
    const uint beauty_width = beauty->getStat().getXres();
    const uint beauty_height= beauty->getStat().getYres();
    assert((beauty_width < sub1_width) && (beauty_height < sub1_height));
    std::cout << "Opening beauty(1): " << watch.restart() << " sec\n";

    std::unique_ptr<IMG_File> subpixel2(nullptr);
    if(subpixel_filename2.compare("") != 0) {
        subpixel2.reset(IMG_File::open(subpixel_filename2.c_str(), parms.get()));
        if(!subpixel2) {
            std::cerr << "Can't open " << subpixel_filename2 << '\n';
            return 1;
        }
        const uint sub2_width = subpixel2->getStat().getXres();
        const uint sub2_height= subpixel2->getStat().getYres();
        assert(sub1_width==sub1_width && sub2_height==sub2_height);
        std::cout << "Opening beauty(2): " << watch.restart() << " sec\n";

    }
    
    // 
    
    auto samples = UT_Vector2D(1,1);
    UT_SharedPtr<UT_Options> opt = subpixel1->imageTextureOptions();
    if (opt->hasOption("image:samples")) {
        samples = opt->getOptionV2("image:samples");
        std::cout << "Subpixels     (1): " << samples.x() << " x " << samples.y() << "\n";
    } else {
        std::cerr << "Can't do without subpixel metadata." << '\n';
        return 1;
    }
    

    //:
    bcd::HistogramParameters accparms;
    bcd::SamplesAccumulator accumulator(beauty_width, beauty_height, accparms);    
    //bcd::SamplesAccumulatorThreadSafe acct(width, height, accparms);
    //
    unsigned long sample_counter = 0;
    UT_Array<PXL_Raster *> images;
    if (!subpixel1->readImages(images, "C")) {
        std::cerr << "Can't read C plane from file " << subpixel_filename1 << '\n'; 
        return 1;
    }
    std::cout << "Pixels read   (1): " << watch.restart() << " sec\n";
    std::unique_ptr<PXL_Raster> C(images[0]);
    // Lets do it nethertheless 
    #if 0
    auto range = UT_BlockedRange<int>(0, beauty_height);
    accumulateThreaded<bcd::SamplesAccumulator>(range, C.get(), beauty_width, 
        samples.x(), samples.y(), &accumulator);

    std::cout << "Accumulated   (1): " << watch.restart() << " sec\n";

    if(subpixel2) {
        // This doesn't clear the memory.
        images.clear(); //?
        // ... but this does:
        C.reset(nullptr);
        if (!subpixel2->readImages(images, "C")) {
            std::cerr << "Can't read C plane from file " << subpixel_filename2 << '\n'; 
            return 1;
        }
        std::cout << "Pixels read   (2): " << watch.restart() << " sec\n";

        C.reset(images[0]);
        auto range = UT_BlockedRange<int>(0, beauty_height);
        accumulateThreaded<bcd::SamplesAccumulator>(range, C.get(), beauty_width, 
            samples.x(), samples.y(), &accumulator);
        std::cout << "Accumulated   (2): " << watch.restart() << " sec\n";
    }

    #else
    auto stat_img = IMG_Stat(beauty_width, beauty_height, IMG_FLOAT32, IMG_RGB);
    auto *plane = stat_img.addPlane("C", IMG_FLOAT32, IMG_RGB);
    auto test_img = std::unique_ptr<IMG_File>(IMG_File::create("./test_resampled.exr", stat_img));

    UT_Array<PXL_Raster *> out;
    if(!beauty->readImages(out, "C")) {
        std::cerr << "Can't read from beauty. \n";
        return 1; 
    }

    // we gonna write to beaty raster already in mem
    std::unique_ptr<PXL_Raster> raster(out[0]);
    // std::cout << beauty

    for (int y=0; y < beauty_height; ++y) {
        const uint Y = y*samples.y();                                              
        for(size_t x=0; x<beauty_width; ++x) {                                                                   
            const uint X = x*samples.x();                                              
            for(size_t sy=0; sy<samples.y(); ++sy) {                                               
                for(size_t sx=0; sx<samples.x(); ++sx) {                                            
                    float color[3] = {0,0,0};                                                             
                    C->getPixelValue(X+sx, Y+sy, color);                                                        
                    accumulator.addSample(y, x, color[0], color[1], color[2]);
                    raster->setPixelValue(x, y, color);
                }                                                                                         
            }                                                                                             
        }
    }

    test_img->writeImages(out);
    test_img->close();
    std::cout << "Test image saved: " << watch.restart() << " sec\n";
    // return 0;

    #endif
    bcd::SamplesStatisticsImages stats = accumulator.extractSamplesStatistics();
    std::cout << "Statistics       : " << watch.restart() << " sec\n";
    bcd::Deepimf histoImage = bcd::Utils::mergeHistogramAndNbOfSamples(
            stats.m_histoImage, stats.m_nbOfSamplesImage);
    
    std::cout << "Historgram       : " << watch.restart() << " sec\n";
    //
    bcd::DenoiserInputs inputs;
    bcd::DenoiserOutputs outputs;
    bcd::DenoiserParameters parameters;
    
    //
    // parameters.m_minEigenValue = 0.0001f;
    
    //
    bcd::Deepimf beauty_image;
    bcd::ImageIO::loadEXR(beauty_image, beauty_filename.c_str());
    inputs.m_pColors      = &beauty_image;
    inputs.m_pNbOfSamples = &(stats.m_nbOfSamplesImage);
    inputs.m_pHistograms  = &(stats.m_histoImage);
    inputs.m_pSampleCovariances = &(stats.m_meanImage);

    //  
    bcd::Deepimf outputDenoisedColorImage(beauty_image);
    outputs.m_pDenoisedColors = &outputDenoisedColorImage;

    //
    auto denoiser = std::make_unique<bcd::MultiscaleDenoiser>(bcd_scales);
    denoiser->setInputs(inputs);
    denoiser->setOutputs(outputs);
    denoiser->setParameters(parameters);
    
    //
    std::cout << "Denoising...  " << "\n";
    denoiser->denoise();
    std::cout << "...done in       : " << watch.restart() << " sec\n";
    //  
    bcd::ImageIO::writeMultiChannelsEXR(outputDenoisedColorImage, output_filename.c_str());
    //
    std::cout << "Image saved : " << watch.restart() << " sec\n";
    stats.m_histoImage.clearAndFreeMemory();
    stats.m_nbOfSamplesImage.clearAndFreeMemory();
    std::cout << "Stats del[] : " << watch.restart() << " sec\n";


    return 0;
}
    // return 0;
    // //const char *hname = "./histogram.exr";
    // //const char *cname = "./covariance.exr";
    
    // //
    // // bcd::ImageIO::writeMultiChannelsEXR(stats.m_meanImage, cname);
    // // bcd::ImageIO::writeMultiChannelsEXR(histoImage, hname);
    
    //
    //std::cout << "Images saved: " << watch.restart() << "sec\n";

  //return 0;

    // po::options_description _cli("bcd options");
    // _cli.add_options()
    //     ("subixel,s",po::value<std::string>()->required(),              "Unfilered super resolution image.")
    //     ("beauty,b", po::value<std::string>()->required(),                "Beauty pass to be denoised.")
    //     ("output,o", po::value<std::string>()->required(),              "Ouput file (.exr)")
    //     // ("skin,k",   po::value<StringVec>()->multitoken(),              "Input skin files  (*.bgeo)")
    //     // ("var,v",    po::value<double>(),                               "PCA Variance (if omitted, PCA won't be performed)")
    //     // ("norm,n",   po::bool_switch()->default_value(false),             "Orthonormalize PCA")
    //     // ("psd,p",    po::bool_switch()->default_value(false),           \
    //         // "Compute pose space deformation (requires tangents vectors)")
    //     ("help,h",                                                     "Prints this screen.");

    // po::variables_map options;        
    // po::store(po::parse_command_line(argc, argv, _cli), options);

    // if (options.count("help") || argc == 1) {
    //     std::cout << _cli << '\n';
    //     return 0;
    // }

    // po::notify(options);

    // bcd::DeepImage<float> input;
    // bcd::ImageIO::loadEXR(input, filename.c_str());

    // std::cout << input.getWidth() << '\n';  
    // std::cout << input.getHeight() << '\n'; 
    // std::cout << input.getDepth() << '\n';  
    // std::cout << input.getSize() << '\n';
    // return 0; 
    // Imf::DeepScanLineInputFile input2(filename.c_str());
    // auto data_window = input2.header().dataWindow();
    // width = data_window.size().x + 1;
    // height = data_window.size().y + 1;
    // // std::cout << width << '\n';  
    // // std::cout << height << '\n';
    // Imf::DeepFrameBuffer buffer;
    // auto pp = std::make_unique<float[]>(width*height*samples);

    // const char *channels[3] = {"R", "G", "B"};
    // for (size_t c=0; c<3; ++c)
    //     buffer.insert(c, 
    //         Imf::DeepSlice(Imf::FLOAT, reinterpret_cast<char*>(pp.get()), width, height, 9, 3, 3));

    // buffer.insert("sample count", 
    //         Imf::DeepSlice(Imf::INT, reinterpret_cast<char*>(pp.get()), width, height, 9, 3, 3));
    
    // input2.setFrameBuffer(buffer);
    // input2.readPixels(reinterpret_cast<char*>(pp.get()), buffer, 0, 1080);


    // const Imf::ChannelList & channel_list = input2.header().channels();
    // if (channel_list.findChannel("C"))
    //     std::cout << "C found \n";
