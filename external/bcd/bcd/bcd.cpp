// bcd
#include "SamplesAccumulator.h"
#include "Utils.h"
#include "DeepImage.h"
#include "ImageIO.h"
#include "Denoiser.h"
#include "MultiscaleDenoiser.h"
#include "IDenoiser.h"
// utils
#include "stopwatch.h"

// HDK
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_DeepShadow.h>
#include <UT/UT_Options.h>
#include <UT/UT_WorkBuffer.h>
#include <UT/UT_Vector2.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_ParallelUtil.h>
#include "bcd.h"
// std
#include <iostream>
#include <memory>
int main(const int argc, const char** argv)
{

    Stopwatch<std::chrono::seconds> watch;
    //
    if (argc < 3) {
        std::cerr << "bcd deepfile.exr beauty.exr denoised.exr\n";
        return 1;
    }

    //
    auto deep_filename   = std::string(argv[1]);
    auto beauty_filename = std::string(argv[2]); 
    auto output_filename = std::string(argv[3]);
    //
    IMG_DeepShadow deep;
    auto parms = std::make_unique<IMG_FileParms>();
    parms->readAlphaAsPlane();
    parms->setDataType(IMG_FLOAT);
    std::unique_ptr<IMG_File> flat(nullptr);
    std::unique_ptr<IMG_File> beauty(nullptr);
    int  width, height;
    //
    if (!deep.open(deep_filename.c_str())) {
        flat.reset(IMG_File::open(deep_filename.c_str()));
        if(!flat) {
            std::cerr << "Can't open " << deep_filename << '\n';
            return 1;
        }
    }

    beauty.reset(IMG_File::open(beauty_filename.c_str()));
    if (!beauty) {
        std::cerr << "Can't open beauty file " << beauty_filename << '\n';
        return 1;
    } else {
        width = beauty->getStat().getXres();
        height= beauty->getStat().getYres();
    }
    

    // 
    std::cout << "Opening file: " << watch.restart() << "sec\n";
    //
    const IMG_DeepShadowChannel *deep_channel = nullptr;
    for (int i = 0; i < deep.getChannelCount(); i++) {
        auto name = std::string(deep.getChannel(i)->getName());
        if( name == std::string("C")) {
            deep_channel = deep.getChannel(i);
            // deep.resolution(width, height);
            break;
        }
    }

    //
    if (!deep_channel && !flat) {
        std::cerr << "No {r,g,b} channel in deep file nor flat file. Quiting now.\n";
        return 1;
    }

    
    auto image_samples = UT_Vector2D(1,1);

    if (deep_channel) {
        UT_SharedPtr<UT_Options> opt = deep.getTextureOptions();
        if (opt->hasOption("image:samples")) {
            image_samples = opt->getOptionV2("image:samples");
            std::cout << "Deep samples: " << image_samples.x() \
            << ", " << image_samples.y() << "\n";
        }
    } 
    else if (flat) {
        // width = flat->getStat().getXres();
        // height= flat->getStat().getYres();
        UT_SharedPtr<UT_Options> opt = flat->imageTextureOptions();
        if (opt->hasOption("image:samples")) {
            image_samples = opt->getOptionV2("image:samples");
            std::cout << "Flat samples: " << image_samples.x() \
            << ", " << image_samples.y() << "\n";
        }
    } else {
        std::cerr << "Neither deep nor flat file could be processed." << '\n';
        return 1;
    }


    //:
    bcd::HistogramParameters accparms;
    bcd::SamplesAccumulator accumulator(width, height, accparms);    
	//bcd::SamplesAccumulatorThreadSafe acct(width, height, accparms);
    //
    unsigned long sample_counter = 0;

    if (flat){

        UT_Array<PXL_Raster *> images;
        if (!flat->readImages(images, "C")) {
            std::cerr << "Can't read C plane from file " << deep_filename << '\n'; 
            return 1;
        }
		
		std::cout << "Pixels read : " << watch.restart() << "sec\n";
        std::unique_ptr<PXL_Raster> C(images[0]);
		// What a hell? Multithread write is not implemented yet in bcd,
		// despite methods exists in headers...
		// Lets do it nethertheless 
		#if 1
		auto range = UT_BlockedRange<int>(0, height);
		accumulateThreaded<bcd::SamplesAccumulator>(range, C.get(), width, 
			image_samples.x(), image_samples.y(), &accumulator);
		std::cout << "Accumulated : " << watch.restart() << "sec\n";
        #else
		for(size_t y=0; y<height; ++y){
            for(size_t x=0; x<width; ++x) {
                for (size_t sy=0; sy<image_samples.y(); ++sy) {
                    for(size_t sx=0; sx<image_samples.x(); ++sx) {
                        float color[3] = {0,0,0};
                        const uint X = x*image_samples.x() + sx;
                        const uint Y = y*image_samples.y() + sy;
                        C->getPixelValue(X, Y, color);
                        accumulator.addSample(y, x, color[0], color[1], color[2]);
                        sample_counter++;
                    }
                }
            }
        }
		std::cout << "Accumulated : " << watch.restart() << "sec\n";
        #endif

    } else if (deep_channel) {
        IMG_DeepPixelReader pixel(deep);
        for(size_t y=height; y>0; --y) {
            for (size_t x=0; x<width; ++x) {
                pixel.open(x, y);
                const int depth = pixel.getDepth();
                for(int s = 0; s<depth; ++s) {
                    const float *data = pixel.getData(*deep_channel, s);
                    accumulator.addSample(y, x, data[0], data[1], data[2]);
                    sample_counter++;
                }
            }
        }
    } else {
       
    }

    //

    std::cout << "width,height: " << width << "," << height << '\n';
    //std::cout << "Samples num : " << sample_counter << '\n';
    //std::cout << "Need samples: " << width*height*image_samples.x()*image_samples.y() << '\n';
    //std::cout << "Is it equal?: " << bool(sample_counter >= width*height*image_samples.x()*image_samples.y()) << '\n';
    //
    //
	bcd::SamplesStatisticsImages stats = accumulator.extractSamplesStatistics();
	std::cout << "Statistics  : " << watch.restart() << "sec\n";
    bcd::Deepimf histoImage = bcd::Utils::mergeHistogramAndNbOfSamples(
            stats.m_histoImage, stats.m_nbOfSamplesImage);
    
	std::cout << "Historgram  : " << watch.restart() << "sec\n";
	//
	bcd::DenoiserInputs inputs;
	bcd::DenoiserOutputs outputs;
	bcd::DenoiserParameters parameters;
	
	//
	//parameters.m_minEigenValue = 0.0001f;
	
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
	auto denoiser = std::make_unique<bcd::MultiscaleDenoiser>(4);
	denoiser->setInputs(inputs);
	denoiser->setOutputs(outputs);
	denoiser->setParameters(parameters);
	
	//
	std::cout << "Denoising...  " << "\n";
	denoiser->denoise();
	std::cout << "...done in  : " << watch.restart() << "sec\n";
	//	
    bcd::ImageIO::writeMultiChannelsEXR(outputDenoisedColorImage, output_filename.c_str());
	//
	std::cout << "Image saved : " << watch.restart() << "sec\n";
	//const char *hname = "./histogram.exr";
    //const char *cname = "./covariance.exr";
    stats.m_histoImage.clearAndFreeMemory();
    stats.m_nbOfSamplesImage.clearAndFreeMemory();
	
	//
	std::cout << "Stats del[] : " << watch.restart() << "sec\n";
    //bcd::ImageIO::writeMultiChannelsEXR(stats.m_meanImage, cname);
    //bcd::ImageIO::writeMultiChannelsEXR(histoImage, hname);
	
	return 0;
	//
	//std::cout << "Images saved: " << watch.restart() << "sec\n";

  //return 0;
}
    // return 0;

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
