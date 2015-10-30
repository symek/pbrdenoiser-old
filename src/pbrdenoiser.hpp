#include <stdexcept>


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
    std::cout<< "\tpbrdenoiser: Denoises Monte-Carlo images" << std::endl << std::endl;
}

void usage()
{
    me();   
    // std::cout<< "Usage: tiit [options] picfile.ext" << std::endl;
    // std::cout<< "options:" << std::endl;
    // std::cout<< "\t -i         check file integrity." << std::endl;
    // std::cout<< "\t -h         print image sha-1 sum." << std::endl;
    // std::cout<< "\t -s         print statistics." << std::endl;
    // std::cout<< "\t -f         fix NaNs and Infs if exist." << std::endl;
    // std::cout<< "\t -p plane   replace working plane (default 'C')." << std::endl;
    // //std::cout<< "\t -w         print image wavelet signature." << std::endl;
    // //std::cout<< "\t -c file    compare wavelet sig of the image with a file." << std::endl;
    // //std::cout<< "\t -S         Suppress output (only minimal data sutable for parasing)." << std::endl;
    // std::cout<< "\t -m         Print meta-data." << std::endl;
    //    std::cout<< "\t -L lut     LUT to be applied on input." << std::endl;
    //    std::cout<< "\t -g 2.2     Apply gamma on input." << std::endl;
    //    std::cout<< "\t -b 16      Convert input bitdepth (1: 8bit, 2: 16bit, 4: 32bit, 16: half, 32: float)." << std::endl;
    //    std::cout<< "\t -o         Output file." << std::endl;

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
