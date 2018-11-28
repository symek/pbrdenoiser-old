#pragma once
//using namespace bcd;
template<typename T>
class AccumulateParallel
{
public:
	AccumulateParallel(
		const PXL_Raster *raster,
		const int width, 
		const int samplex, 
		const int sampley,
		T * accumulator)
		:
		raster(raster),
		accumulator(accumulator),
		width(width),
		samplex(samplex),
		sampley(sampley) {}

	void operator()(const UT_BlockedRange<int> &r) const {
		for (int y=r.begin(); y != r.end(); ++y) {
			for(size_t x=0; x<width; ++x) {                                                                   
				for(size_t sy=0; sy<sampley; ++sy) {                                               
					for(size_t sx=0; sx<samplex; ++sx) {                                            
						const uint X = x*samplex + sx;                                              
						const uint Y = y*sampley + sy;                                              
						float color[3] = {0,0,0};                                                             
						raster->getPixelValue(X, Y, color);                                                        
						accumulator->addSample(y, x, color[0], color[1], color[2]);
					}                                                                                         
				}                                                                                             
			}                                                                                                 
		}
	}

private:
	const PXL_Raster * raster;
	const int width;
	const int samplex;
	const int sampley;
	T * accumulator;
};

template <typename T>
void accumulateThreaded(const UT_BlockedRange<int> &range, const PXL_Raster *raster, \
    const int width, const int sx, const int sy, 
	T *accumulator)
{
    UTparallelFor(range, AccumulateParallel<T>(raster, width, sx, sy, accumulator));
}



