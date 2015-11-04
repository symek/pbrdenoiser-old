//HDK:
#include <IMG/IMG_File.h>
#include <IMG/IMG_Stat.h>
#include <IMG/IMG_FileParms.h>
#include <IMG/IMG_Format.h>
#include <CMD/CMD_Args.h>
#include <UT/UT_PtrArray.h>
#include <PXL/PXL_Raster.h>
#include <UT/UT_String.h>
#include <tbb/task_scheduler_init.h>

// Houdini version info
#include <SYS/SYS_Version.h>
#if SYS_VERSION_MAJOR_INT == 14
#define HOUDINI_14
#endif

// std cmath (to consider for isnan() and isinf())
#include <cmath>
#include <fstream>

// std:
#include <string.h>
#include <time.h>

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

SampleSet::SampleSet( const std::vector< Image > &images, const std::vector< Image > &motion) :
    m_width(0),
    m_height(0),
    m_motion(motion)
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
                int mx = 0; int my = 0;
                for( unsigned int i = 0; i < images.size(); ++i )
                {
                    s[i] = images[i].readable( x + mx, y + my )[c];
                    mx += static_cast<int>(m_motion[i].at(x, y)[0]);
                    my += static_cast<int>(m_motion[i].at(x, y)[1]);
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


int read_image(const std::string &inputName, const std::string &plane_name, Image &image, bool verbose=true)
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
    if (verbose)
        printBasicInfo(inputFile);

    /// Under this line all routines will require myData
    UT_PtrArray<PXL_Raster *> images; // arrays of rasters
    bool       loaded        = false; // did we load the file
    PXL_Raster *raster       = NULL;  // working raster 
    loaded = inputFile->readImages(images);
    if (!loaded)
    {
        std::cout << "Integrity : Fail" << std::endl;
        return 0;
    } 
    else 
    {
        if(verbose)
            std::cout << "Integrity : Ok" << std::endl;
        int px = stat.getPlaneIndex(plane_name.c_str());
        if (px != -1)
        {
            raster  = images(px);
        }
        else
        {
            std::cerr << "No raster :" << plane_name << std::endl;
            return 0;
        }
    }  

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


void getPlaneNames(std::string &tokens, std::vector<std::string> &planeNames)
{
    std::istringstream stream(tokens);
    std::string token;
    while(std::getline(stream, token, ',')) 
    {
        planeNames.push_back(token);
    }
}

int
main(int argc, char *argv[])
{   
    // Profiling timer:
    Timer c = Timer();
    /// No options? Exit:  
    if (argc == 1)
    {
        usage();
        return 0;
    }
 
    int nthreads = tbb::task_scheduler_init::default_num_threads();
    std::cout << "Working with: " << nthreads << " threads." << std::endl;

    /// Read argumets and options:
    CMD_Args args;
    std::string plane_str = "C";            // Default working plane.
    std::string vel_plane_str;              // Motion vector pass.
    std::vector<std::string> planes_vector; // Planes vector in case user wants to filter > 1 plane
    planes_vector.push_back(plane_str);     // dito
    std::vector<Image> motion;              // Storage for motion vector pass 
                                            // (we will reuse for all planes thus I pref to keep it here.)

    args.initialize(argc, argv);
    args.stripOptions("m:b:s:k:f:p:v:");
    Options opt;

    // Files we work on:
    // FIXME: I dont know how to use CMD_Args to with multiply options 
    //  (argp('i', number ) doesnt work for me)
    std::vector<  std::string > image_names;
    for (int i = 1; i < argc; ++i)
    {
        if (std::ifstream(argv[i]).good()) // Ignore all non-files...
            image_names.push_back(std::string(argv[i]));
    }

    // Parse options:
    if (args.found('m'))
        opt.blurMode = args.iargp('m');
    if(args.found('b'))
        opt.blurStrength = args.fargp('b');
    if(args.found('s'))
        opt.contributionStrength = args.fargp('s');
    if (args.found('k'))
        opt.kernelWidth = args.fargp('k');
    if (args.found('f'))
        opt.startFrame = args.iargp('f');
    if (args.found('p'))
    {
        plane_str = std::string(args.argp('p'));
        planes_vector.clear();
        getPlaneNames(plane_str, planes_vector);
    }
    if(args.found('v'))
        vel_plane_str = std::string(args.argp('v'));

     // Make a name for an image to be genarated:
    std::string name      = image_names[opt.startFrame];
    int length            = name.length();
    int lastindex         = name.find_last_of("."); 
    std::string rawname   = name.substr(0, lastindex);
    std::string extension = name.substr(lastindex, length);
    std::string outputName = rawname + ".filtered" + extension;

    // Info:
    std::cout << "Output to : " <<  outputName << std::endl;

    /* ---- EXPORTING RASTER TO FILE ---- */
    IMG_FileParms parms = IMG_FileParms();
    parms.setDataType(IMG_FLOAT);                               // Readin conversion to float
    IMG_File *inputFile = IMG_File::open(name.c_str(), &parms); // 'name' is our working frame atm.
    static const IMG_Stat &stat = inputFile->getStat();         // Source stats
                 IMG_Stat ostat = IMG_Stat(stat);               // Copy source stat to be modfied.
    UT_PtrArray<PXL_Raster *> images;                           // arrays of rasters we will copy our output onto.
    PXL_Raster *raster          = NULL;                         // working raster         


    if(vel_plane_str.length() > 0)
    {
        motion.resize(image_names.size());
        int result = 0;
        for (int i = 0; i < image_names.size(); ++i)
            result += read_image(image_names[i], vel_plane_str, motion[i], false);
        if (result != image_names.size())
        {
            vel_plane_str = ""; // turn off if not all files have motion vectors.
            std::cerr << "Can't apply motion vectors. Some/all files are missing velocity plane." << std::endl;
        }
        else
            std::cout << "Applying motion vectors." << std::endl;
    }

    // for every plane in planes_vector...

    for (std::vector<std::string>::iterator i = planes_vector.begin();\
        i != planes_vector.end(); ++i)
    {
        std::string plane(*i);
        // Vector of frames to be processed:
        std::vector<Image> frames(image_names.size());
        int result = 0;
        for (int i = 0; i < image_names.size(); ++i)
            result += read_image(image_names[i], plane, frames[i]);

        // How many frames we opened:
        opt.nImages = result;
        // Start from building SampleSet:
        c.start();
        SampleSet time_samples(frames, motion);
        std::cout << "Building set for [" << plane << "] took " << c.current() << " seconds." << std::endl;
        // This will hold output to be copined back to working place later on:
        const int width = time_samples.width(), height = time_samples.height();
        Image output( width, height );
        
        // DENOISE FUNCTION:
        c.start();
        # ifdef NO_THREADING
            denoise(time_samples, output, opt); // Single thread version.
        #else
            denoiserThreaded(time_samples, opt, height, output);
        #endif

        std::cout << std::endl;
        std::cout << "Denoising of [" << plane << "] finished in: " << c.current() / nthreads << " seconds." << std::endl;


        bool loaded  = inputFile->readImages(images);        
        // this shouldn't happend as we already opened that file...
        if (!loaded) 
            return 1;
        // This also shoulnd't happen as plane was already found in read_images()...
        int px = stat.getPlaneIndex(plane.c_str());
        if (px != -1)
            raster  = images(px);
        else
        {
            std::cerr << "No raster :" << plane << std::endl;
            return 1;
        }
        

        // Copy data back to raster of intereset...
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                float pixel[3];
                const double *source = output.at(x, y);
                for(int i = 0; i < 3; ++i) pixel[i] =  static_cast<float>(source[i]);
                #ifdef HOUDINI_14 // Houdini 14 doesn't have nice setPixelValue():
                    float *target = (float*)raster->getPixel(x, y, 0);
                    for(int i = 0; i < 3; ++i) 
                        target[i] =  pixel[i];
                #else
                    raster->setPixelValue(x, y, pixel);
                #endif
            }
        }

    }

    // OutputFile should be a copy of working frame with modified raster:
    ostat.setFilename(outputName.c_str());
    parms.setDataType(IMG_HALF); // Back to half when possible.
    IMG_File *outputFile = IMG_File::create(outputName.c_str(), (const IMG_Stat)ostat, &parms);

    // End
    if (outputFile)
        outputFile->writeImages(images);
    else
        return 1;
    return 0;    
}