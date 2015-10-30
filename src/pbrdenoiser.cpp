//HDK:
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_Format.h>
#include <CMD/CMD_Args.h>
#include <UT/UT_PtrArray.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_String.h>

// std cmath (to consider for isnan() and isinf())
#include <cmath>
#include <fstream>

// std:
#include <string.h>

// Self:
#include "pbrdenoiser.hpp"

Image::Image( int w, int h )
{
    resize( w, h );
}

const double* Image::readable( int x, int y ) const
{
    x = std::max( std::min( x, m_width - 1 ), 0 );
    y = std::max( std::min( y, m_height - 1 ), 0 );
    return &m_data[ ( x + m_width * y ) * 3 ];
}

double* Image::writeable( int x, int y )
{
    return &m_data[ ( x + m_width * y ) * 3 ];
}

const double* Image::at( int x, int y ) const
{
    return &m_data[ ( x + m_width * y ) * 3 ];
}

SampleSet::SampleSet( const std::vector< Image > &images ) :
    m_width(0),
    m_height(0)
{
    for( unsigned int j = 0; j < images.size(); ++j )
    {
        if( m_width == 0 && m_height == 0 )
        {
            m_width = images[j].width();
            m_height = images[j].height();
        }
        else
        {
            if( m_width != images[j].width() || m_height != images[j].height() )
            {
                throw std::runtime_error( "Not all images are the same size." );
            }
        }
    }
    m_samples.resize( m_width * m_height * 3 );

    const int arraySize = m_width * m_height * 3;   
    m_mean.resize( arraySize );
    m_median.resize( arraySize );
    m_variance.resize( arraySize );
    m_deviation.resize( arraySize );
    m_min.resize( arraySize );
    m_max.resize( arraySize );

    int index = 0;
    std::vector<double> s;
    for( int y = 0; y < m_height; ++y )
    {
        for( int x = 0; x < m_width; ++x )
        {
            for( int c = 0; c < 3; ++c, ++index )
            {
                s.resize( images.size() );
                for( unsigned int i = 0; i < images.size(); ++i )
                {
                    s[i] = images[i].at( x, y )[c];
                }
            
                double mean = 0., variance = 0., norm = 1. / s.size();
                double max = std::numeric_limits<double>::min();
                double min = std::numeric_limits<double>::max();

                double areBlack = true;
                for( unsigned int i = 0; i < s.size(); ++i )
                {
                    if( s[i] != 0. )
                    {
                        areBlack = false;
                        break;
                    }
                }
                
                if( !areBlack )
                {
                    bool hasAnyBlack = false;
                    double count = 0.;
                    for( unsigned int i = 0; i < s.size(); ++i )
                    {
                        if( s[i] == 0. )
                        {
                            hasAnyBlack = true;
                            continue;
                        }

                        ++count;

                        min = ( s[i] < min ) ? s[i] : min;
                        max = ( s[i] > max ) ? s[i] : max;
                        mean += s[i];
                    }
                    mean /= count;

                    if( hasAnyBlack )
                    {
                        for( unsigned int i = 0; i < s.size(); ++i )
                        {
                            if( s[i] == 0. )
                            {
                                continue;
                            }
                            double d = s[i] - mean;
                            variance += d*d*( 1. / count );
                        }
                    
                        bool flipFlop = false;
                        double newMean = 0.;
                        variance = 0.;
                        for( unsigned int i = 0; i < s.size(); ++i )
                        {
                            if( s[i] == 0. )
                            {
                                s[i] = mean;
                                newMean += mean * norm;
                            }
                            double d = s[i] - mean;
                            variance += d*d*norm;
                        }
                        mean = newMean;
                    }
                    else
                    {
                        for( unsigned int i = 0; i < s.size(); ++i )
                        {
                            double d = s[i] - mean;
                            variance += d*d*norm;
                        }
                    }
                }
                else
                {
                    min = max = mean = variance = 0.;
                }
                
                m_min[index] = min;
                m_max[index] = max;
                m_mean[index] = mean;
                m_variance[index] = variance;
                m_deviation[index] = sqrt( variance );
                m_samples[index] = s;
                
                std::sort( s.begin(), s.end() );
                if( s.size() > 2 )
                {
                    m_median[index] = s[ std::min( int( s.size()-1 ), std::max( int( ceil( s.size() / 2. ) ), 0 ) ) ];
                }
                else
                {
                    m_median[index] = s.front() + ( s.back() - s.front() ) / 2.;
                }
            }
        }
    }
}


int read_image(const std::string &inputName, Image &image)
{
    // Optional read parameters:
    IMG_FileParms parms = IMG_FileParms();
    /// At first we just open the file to read stat:
    IMG_File *inputFile = NULL;
   
    /// Read file:
    inputFile = IMG_File::open(inputName.c_str(), &parms);
    if (!inputFile)
    { 
        std::cerr << "Can't open: " << inputName << std::endl;
        return 0;
    }
    
    /// Print basic info:
    static const IMG_Stat &stat = inputFile->getStat();
    printBasicInfo(inputFile);

    /// Under this line all routines will require myData
    UT_PtrArray<PXL_Raster *> images; // arrays of rasters
    bool       loaded        = false; // did we load the file
    const char *plane_name   = "C"  ; // default raster name
    PXL_Raster *raster       = NULL;  // working raster 
    void       *myData       = NULL;  // data ready to be casted
    int         npix         = 0;  


    loaded = inputFile->readImages(images);
    if (!loaded)
    {
        std::cout << "Integrity : Fail" << std::endl;
        return 0;
    } 
    else 
    {
        std::cout << "Integrity : Ok" << std::endl;
        // TODO: default C could not exists!
        // Look for it, and choose diffrent if neccesery
        int px = stat.getPlaneIndex(plane_name);
        if (px != -1)
        {
            raster  = images(px);
            myData  = raster->getPixels();
            npix    = raster->getNumPixels();
        }
        else
        {
            std::cerr << "No raster :" << plane_name << std::endl;
            return 0;
        }
    }  

    float *fpixels;
    fpixels = static_cast<float *>(myData);
    int xres = raster->getXres();
    int yres = raster->getYres();

    image.resize( xres, yres );
    float pixel[3];

    for( int y = 0; y < image.height(); ++y )
    {
        for( int x = 0; x < image.width(); ++x )
        {
            raster->getPixelValue(x, y, pixel);
            double *inPixel = image.writeable( x, y );
            inPixel[0] = static_cast<double>(pixel[0]);
            inPixel[1] = static_cast<double>(pixel[1]);
            inPixel[2] = static_cast<double>(pixel[2]);
        }
    }

    return 1;
}


void denoise(SampleSet &set, Image &result, Options &opt)
{
    //===================================================================
    // The algorithm.
    //===================================================================


    const int width = set.width(), height = set.height();

    int kernelRadius = opt.kernelWidth > 1 ? ( opt.kernelWidth - 1 ) / 2 : 0;
    std::vector< double > srcSamples;
    for( int y = 0; y < height; ++y )
    {
        for( int x = 0; x < width; ++x )
        {
            fprintf(stderr,"\rFiltering %5.2f%% complete.", 100. * y / ( height-1 ) );
        
            // Loop over each channel.  
            for( unsigned int c = 0; c < 3; ++c )
            {
                double destMean = set.mean( x, y, c );
                double destDeviation = set.deviation( x, y, c );
                double destVariation = set.variance( x, y, c );
                double destRange = set.max( x, y, c ) - set.min( x, y, c ); 

                // Loop over the neighbouring pixels.
                double weightedSum = 0.;
                double v = 0.;
                for( int ky = -kernelRadius; ky <= kernelRadius; ++ky )
                {
                    for( int kx = -kernelRadius; kx <= kernelRadius; ++kx )
                    {
                        // Don't include the pixel being sampled in our calculations as we are
                        // summing the deviations from it and doing so will bias our results.
                        if( ky == 0 && kx == 0 )
                        {
                            continue;
                        }

                        // Gather information on the source pixel's samples.
                        srcSamples = set.samples( x + kx, y + ky, c );
                        double srcMin = set.min( x + kx, y + ky, c );
                        double srcMax = set.max( x + kx, y + ky, c );
                        double srcMean = set.mean( x + kx, y + ky, c );
                        double srcDeviation = set.deviation( x + kx, y + ky, c );
                        double srcVariation = set.variance( x + kx, y + ky, c );
                        double srcRange = set.max( x + kx, y + ky, c ) - set.min( x + kx, y + ky, c ); 
                        
                        if( srcVariation == 0 && srcSamples[0] == 0. ) continue;
                            
                        // A gaussian falloff that weights contributing samples which are closer to the pixel being filtered higher.
                        /// \todo Intuitive falloff parameters need to be added to the distance weight or at least a suitable curve found.
                        double distanceWeight = gaussian( sqrt( kx*kx + ky*ky ) / sqrt( kernelRadius*kernelRadius + kernelRadius*kernelRadius ), 0., .7, false );
                            
                        // Similarity weight.
                        // This weight defines a measure of how similar the set of contributing samples is to the pixel being filtered.
                        // By itself it will produce a smart blur of sorts which is then attenuated by the variance of the source samples in the process of weighted offsets.
                        // Changing this value will effect how aggressive the filtering is.
                        double similarity;
                        if( opt.blurMode == Options::kAggressive )
                        {
                            similarity = ( srcMean - destMean ) * ( srcRange - destRange );
                        }
                        else
                        {
                            similarity = ( srcMean - destMean );
                        }
                        similarity *= similarity;

                        // Temporal weight.
                        // Weight the contribution using a function in the range of 0-1 which weights the importance of the
                        // contributing sample according to how close it is in time to the current time.
                        double time = 1.; // \todo: implement this! Example functions are Median, Gaussian, etc.

                        // Loop over each of the neighbouring samples.
                        for( unsigned int i = 0; i < srcSamples.size(); ++i )
                        {
                            // The contribution weight extends the range of allowed samples that can influence the pixel being filtered.
                            // It is simply a scaler that increases the width of the bell curve that the samples are weighted against.
                            double contribution = gaussian( srcSamples[i], destMean, destDeviation * ( 1 + opt.contributionStrength ) ) * gaussian( srcSamples[i], srcMean, srcDeviation );
                            contribution = contribution * ( 1. - opt.blurStrength ) + opt.blurStrength;

                            // This weight is a step function with a strong falloff close to the limits. However, it will never reach 0 so that the sample is not excluded.
                            // By using this weight the dependency on the limiting samples is much less which reduces the effect of sparkling artefacts.
                            double limitWeight = srcSamples.size() <= 2 ? 1. : softStep( srcSamples[i], srcMin, srcMax );
                        
                            // Combine the weights together and normalize to the range of 0-1.  
                            double weight = pow( M_E, -( similarity / ( contribution * srcVariation * time * distanceWeight * limitWeight ) ) );
                            weight = ( isnan( weight ) || isinf( weight ) ) ? 0. : weight;
                        
                            // Sum the offset.  
                            v += ( srcSamples[i] - destMean ) * weight;

                            // Sum the weight.
                            weightedSum += weight;
                        }
                    }
                }

                if( weightedSum == 0. || destVariation <= 0. )
                {
                    result.writeable( x, y )[c] = destMean;
                }
                else
                {
                    result.writeable( x, y )[c] = destMean + ( v / weightedSum );
                }
            }
        }
    }
}


int
main(int argc, char *argv[])
{   
    /// No options? Exit:  
    if (argc == 1)
    {
        usage();
        return 0;
    }
  
    /// Read argumets and options:
    CMD_Args args;
    args.initialize(argc, argv);
    args.stripOptions("p:c:w:L:o:ihsfSmg:b:");
    Options opt;

    /// File we work on:/
    std::vector<  std::string > image_names;
    for (int i = 1; i < argc; ++i)
        image_names.push_back(std::string(argv[i]));

    // Print info:
    for ( std::vector<std::string>::iterator i = image_names.begin();
        i != image_names.end(); i++ ) 
    {
        std::cout << *i << std::endl;
    }
    
    std::vector<Image> frames(image_names.size());

    int result = 0;
    for (int i = 0; i < image_names.size(); ++i)
    {
        result += read_image(image_names[i], frames[i]);
    }

    opt.nImages = result;
    SampleSet time_samples(frames);
    const int width = time_samples.width(), height = time_samples.height();
    // Image output( width, height );
    Image output = frames[1];

    // denoise(time_samples, output, opt);

    std::string name      = image_names[0];
    int length            = name.length();
    int lastindex         = name.find_last_of("."); 
    std::string rawname   = name.substr(0, lastindex);
    std::string extension = name.substr(lastindex, length);
    std::string outputName = rawname + ".filtered" + extension;

    std::cout << "Output to : " <<  outputName << std::endl;

    IMG_FileParms parms  = IMG_FileParms();
    IMG_File *outputFile = NULL;

   
    /// Copy settings from input:
    // IMG_File *inputFile = NULL;
    // inputFile = IMG_File::open(name.c_str());
    const IMG_Stat stat = IMG_Stat(width, height, IMG_FLOAT, IMG_RGBA);
   
    /// Read file:
    outputFile = IMG_File::create(outputName.c_str(), stat, &parms);
    if (!outputFile)
    { 
        std::cerr << "Can't write to: " << outputName << std::endl;
        return 1;
    }

    void *buffer   = outputFile->allocScanlineBuffer();
    float *fbuffer = static_cast<float *>(buffer);

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            // float pixel[3];
            const double *source = output.at(x, y);
            fbuffer[3*x]   = static_cast<float>(source[0]);
            fbuffer[3*x+1] = static_cast<float>(source[2]);
            fbuffer[3*x+2] = static_cast<float>(source[3]);
            fbuffer[3*x+3] = 1.0f;
            std::cout << source[0] << ", " << source[1]<< ", " << source[2] << std::endl;
        }

        outputFile->write(y, buffer, 0);

    }


    outputFile->close();

    delete fbuffer;
    return 0;    

}



















    
    // /// Switch working plane if requested 
    // if (args.found('p')) 
    // {
    //     plane_name = args.argp('p');
    //     std::cout << "Work plane: " << plane_name << std::endl;
    // }

    // /// 'Integrity check' basicaly means we try to load 
    // /// the file into memory: fail == it's broken.      
    // if (args.found('i') || args.found('h') || args.found('s') ||
    //     args.found('c') || args.found('w') || args.found('f'))
    // {
    //     loaded = inputFile->readImages(images);
    //     if (!loaded)
    //     {
    //         std << "Integrity : Fail" << std::endl;
    //         return 1;
    //     } 
    //     else 
    //     {
    //         cout << "Integrity : Ok" << std::endl;
    //         // TODO: default C could not exists!
    //         // Look for it, and choose diffrent if neccesery
    //         int px = stat.getPlaneIndex(plane_name);
    //         if (px != -1)
    //         {
    //             raster  = images(px);
    //             myData  = raster->getPixels();
    //             npix    = raster->getNumPixels();
    //         }
    //         else
    //         {
    //             cerr << "No raster :" << plane_name << std::endl;
    //             return 1;
    //         }
    //     }
    // }

    

    // /// Print statistics and optionally fix nans/infs: 
    // /// FIXME: fix doesn't work yet.
    // if (args.found('s') || args.found('f'))
    // {
    //     /// I'm acctually reloading the image with 32float,
    //     //float *fpixels;
    //     myData = printStats(inputName, plane_name,  args.found('f'));
    // }
    
    // /// Meta data:
    // if (args.found('m'))
    // {
    //     UT_String info = "";
    //     inputFile->getAdditionalInfo(info);
    //     cout << info.buffer()  << std::endl;
    // }

    // if (args.found('o'))
    // {
    //     const char *outputName = NULL;
    //     outputName = args.argp('o');
    //     if (outputName)
    //         inputFile->copyToFile(inputName, outputName, parms);

    // }





 //    if (args.found('L'))
 //    {
 //        lut_file = args.argp('L');
 //        if (lut_file)
 //             parms->applyLUT(lut_file, "C");
 //    }
 //    // Gamma:
 //    if (args.found('g'))
 //    {
 //        const char *gamma = args.argp('g');
 //        parms->applyGamma(atof(gamma), "C");
 //    }
 //    // Bit depth:
 //    if (args.found('b'))
 //    {
 //        const char *bitdepth = args.argp('b');
 //        parms->setDataType(getDataType(atoi(bitdepth)));
 //    }