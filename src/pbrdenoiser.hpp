//////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2014, Luke Goddard. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining
//  a copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom
//  the Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//////////////////////////////////////////////////////////////////////////

#include <stdexcept>
#include <GA/GA_SplittableRange.h>
#include <GA/GA_Range.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>


struct Options
{
    Options() :
        blurMode( kAggressive ), 
        nImages( 5 ),
        blurStrength( .005 ),
        contributionStrength( 5 ),
        kernelWidth( 7 ),
        sequenceNumber( 0 ),
        startFrame( 0 )
    {
    }

    enum
    {
        kAggressive,
        kGentle
    };

    int blurMode;   
    int nImages;
    double blurStrength;
    double contributionStrength;
    int kernelWidth;
    int sequenceNumber;
    int startFrame;
};



inline int fromGamma22(double x)
{
    x = std::min( std::max( 0., x ), 1. );
    int v = int( pow( x, 1/2.2 ) * 255 + .5 );
    v = v > 255 ? 255 : v;
    v = v < 0 ? 0 : v;
    return v;
}

inline double toGamma22(int x)
{
    return pow( ( ( double( x ) ) / 255 ), 2.2 );
}

struct Image
{
    public :
    
        Image( int w = 1, int h = 1 );

        const double* readable( int x, int y ) const;
        double* writeable( int x, int y );
        const double* at( int x, int y ) const;
        inline int width() const { return m_width; };
        inline int height() const { return m_height; };
        
        inline void resize( int width, int height )
        {
            if( width < 1 || height < 1 )
            {
                throw std::runtime_error( "Cannot resize an image to null dimensions." );
            }

            m_width = width;
            m_height = height;
            m_data.resize( width * height * 3 );
        }

    private :
    
        int m_width, m_height;
        std::vector<double> m_data;
};

struct SampleSet
{
    public :

        SampleSet( const std::vector< Image > &i );

        inline int width() const { return m_width; };
        inline int height() const { return m_height; };
        inline const std::vector<double> &samples( int x, int y, int c ) const { return m_samples[ arrayIndex( x, y, c ) ]; }
        inline double mean( int x, int y, int c ) const { return m_mean[ arrayIndex( x, y, c ) ]; };
        inline double max( int x, int y, int c ) const { return m_max[ arrayIndex( x, y, c ) ]; };
        inline double min( int x, int y, int c ) const { return m_min[ arrayIndex( x, y, c ) ]; };
        inline double median( int x, int y, int c ) const { return m_median[ arrayIndex( x, y, c ) ]; };
        inline double variance( int x, int y, int c ) const { return m_variance[ arrayIndex( x, y, c ) ]; };
        inline double deviation( int x, int y, int c ) const { return m_deviation[ arrayIndex( x, y, c ) ]; };
        inline double midpoint( int x, int y, int c ) const
        {
            double mn = min( x, y, c );
            double mx = max( x, y, c );
            return ( mx - mn ) * .5 + mn;
         };

    private :

        inline int arrayIndex( int x, int y, int c ) const
        { 
            x = std::max( std::min( x, m_width - 1 ), 0 );
            y = std::max( std::min( y, m_height - 1 ), 0 );
            c = std::max( std::min( c, 3 ), 0 );
            return ( y * m_width + x ) * 3 + c;
        }

        int m_width, m_height;
        std::vector< std::vector< double > > m_samples;
        std::vector< double > m_mean, m_variance, m_deviation, m_min, m_max, m_median;
};

struct BmpHeader
{
    unsigned int   mFileSize;        // Size of file in bytes
    unsigned int   mReserved01;      // 2x 2 reserved bytes
    unsigned int   mDataOffset;      // Offset in bytes where data can be found (54)

    unsigned int   mHeaderSize;      // 40B
    int    mWidth;           // Width in pixels
    int    mHeight;          // Height in pixels

    short  mColorPlates;     // Must be 1
    short  mBitsPerPixel;    // We use 24bpp
    unsigned int   mCompression;     // We use BI_RGB ~ 0, uncompressed
    unsigned int   mImageSize;       // mWidth x mHeight x 3B
    unsigned int   mHorizRes;        // Pixels per meter (75dpi ~ 2953ppm)
    unsigned int   mVertRes;         // Pixels per meter (75dpi ~ 2953ppm)
    unsigned int   mPaletteColors;   // Not using palette - 0
    unsigned int   mImportantColors; // 0 - all are important
};


/// Returns the Gaussian weight for a point 'x' in a normal distribution
/// centered at the mean with the given deviation.
double gaussian( double x, double mean, double deviation, bool normalize = true )
{
    double c = deviation;
    double a = c * sqrt( 2 * M_PI );
    double b = x - mean;
    double v = a * exp( -( ( b * b ) / ( 2*c*c ) ) );
    return normalize ? v / a : v;
}

/// A step function which has a falloff that starts when the value 'x'
/// gets within 10% of a limit. The returned value will never reach 0
/// if it is within range. Values of 'x' that are out of the range 
/// are set to 0.
double softStep( double x, double min, double max )
{
    if( x < min || x > max )
    {
        return 0;
    }

    double v = ( max - min ) * .1;
    double lower = min + v;
    double upper = max - v;
    if( x < min + v )
    {
        x = ( x - min ) / ( lower - min );
    }
    else if( x > max - v )
    {
        x = 1. - ( x - upper ) / ( max - upper );
    }
    else
    {
        return 1.;
    }

    // We ensure that the weight returns a contribution of at least .0025;
    x = .05 + ( x * .95 );

    return x * x;
}


void me()
{
    std::cout<< std::endl;
    std::cout<< "\tpbrdenoiser: Denoises Monte-Carlo images.\n\tThis code is based on Luke Goddard's temporal denoiser: https://github.com/goddardl/temporalDenoise"<< std::endl;
}

void usage()
{
    me();   
    std::cout<< "\tUsage: pbrdenoiser [options] /path/image.[1-5].ext" << std::endl;
    std::cout<< "\toptions:" << std::endl;
    std::cout<< "\t -m      Blur mode (0 agressive (default), 1 gentle)." << std::endl;
    std::cout<< "\t -b      Blur strength (default 0.005)." << std::endl;
    std::cout<< "\t -s      Contribution strength (default 5)." << std::endl;
    std::cout<< "\t -k      Kernel width (default 7)." << std::endl;
    std::cout<< "\t -f      Current frame (default 0 - first frame of provided sequence)." << std::endl;
    std::cout<< "\t -p      Plane to process (default C)." << std::endl << std::endl;
}



const char * getDataTypeName(const int type)
{
    const char * typeName = "";
    switch(type)
    {
        case 8:  typeName = "32float"; break;
        case 16: typeName = "16float"; break;
        case 1:  typeName = "8Int";    break;
        case 2:  typeName = "16Int";   break;
        case 4:  typeName = "32Int";   break;
    }
    return typeName;
}

IMG_DataType getDataType(const int type)
{
    switch(type)
    {
        case 1:  return IMG_INT8;
        case 2:  return IMG_INT16;
        case 4:  return IMG_INT32;
        case 16: return IMG_HALF;
        case 32: return IMG_FLOAT;
    }
    return IMG_DT_ANY;
}

/* Prints statistics*/
void printBasicInfo(IMG_File *file)
{
    
    const IMG_Stat &stat = file->getStat();
    const int numPlanes  = stat.getNumPlanes();
    std::cout <<"Filename  : " << stat.getFilename() << std::endl;
    printf("Resolution: %ix%i\n", stat.getXres(), stat.getYres());
    printf("Rasters   : ");
    
    for (int i = 0; i < numPlanes; i++)
    {
        IMG_Plane *plane = stat.getPlane(i);
        std::cout << "";
        printf("%s/%s/", plane->getName(), getDataTypeName(plane->getDataType()));
        for (int j = 0; j < 3; j++) 
        {
            const char *compName = plane->getComponentName(j);
            if (compName) printf("%s", compName );
        }
        std::cout << ", ";
    }
    printf("\n");
}


class op_Denoiser {
public:
    op_Denoiser(const SampleSet &set, const Options &opt, Image &image)
        : set(set),  result(image), opt(opt) {};
    void    operator()(const UT_BlockedRange<int64> &range) const
            {
                //std::cout << " Progressing from " << range.begin() << " to " << range.end() << std::endl;
                const int width = set.width(), height = set.height();
                int kernelRadius = opt.kernelWidth > 1 ? ( opt.kernelWidth - 1 ) / 2 : 0;
                std::vector< double > srcSamples;
                fprintf(stderr,"\rFiltering %5.2f%% complete.", 100.0 * range.begin()  / ( height-1 )); 
                // Iterate over scanlines in image:
                for (int y = range.begin(); y != range.end(); ++y)
                {
                    for( int x = 0; x < width; ++x )
                    {
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
    private:
            const SampleSet &set;
            const Options &opt;
            Image &result;
};

void
denoiserThreaded(const SampleSet &set,  const Options &opt, const int height, Image &output)
{
   
    UTparallelFor(UT_BlockedRange<int64>(0, height), op_Denoiser(set, opt, output), 2, opt.kernelWidth*2);
}


class Timer 
{
    private:
        double begTime;
    public:
        void start()
        {
            begTime = clock();
        }

        double current() 
        {
        //int threads = UT_Thread::getNumProcessors();
            return (difftime(clock(), begTime) / CLOCKS_PER_SEC);// / (double) threads;
        }

        bool isTimeout(double seconds)
        {
            return seconds >= current();
        }
};