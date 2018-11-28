#pragma once
//using namespace bcd;
class AccumulateParallel
{
public:
	AccumulateParallel(const PXL_Raster *raster,
					   bcd::SamplesAccumulator * accumulator,
					   int width, /*int height,*/ int samplex, int sampley):
					   raster(raster),
					   accumulator(accumulator),
					   width(width),
					   //height(height),
					   samplex(samplex),
					   sampley(sampley)
					   {}

	void operator()(const UT_BlockedRange<int> &r) const {
		//GA_Offset start;
		//GA_Offset end;
		for (int y=r.begin(); y != r.end(); ++y) {
			//size_t y      = static_cast<size_t>(start);
			//size_t height = static_cast<size_t>(end); 
			//for(; y<height; ++y){                                                                       
				for(size_t x=0; x<width; ++x) {                                                                   
					for (size_t sy=0; sy<sampley; ++sy) {                                               
						for(size_t sx=0; sx<samplex; ++sx) {                                            
							const uint X = x*samplex + sx;                                              
							const uint Y = y*sampley + sy;                                              
							float color[3] = {0,0,0};                                                             
							raster->getPixelValue(X, Y, color);                                                        
							accumulator->addSample(y, x, color[0], color[1], color[2]);                            
							//sample_counter++;                                                                     
						}                                                                                         
					}                                                                                             
				}                                                                                                 
			//}  		
		}
	}

private:
	const PXL_Raster * raster;
	bcd::SamplesAccumulator * accumulator;
	int width;
	//int height;
	int samplex;
	int sampley;
};


void
accumulateThreaded(const UT_BlockedRange<int> &range, const PXL_Raster *raster, \
    const int width, const int sx, const int sy, 
	bcd::SamplesAccumulator *acc)
{
    // Create a GA_SplittableRange from the original range
    //GA_SplittableRange split_range = GA_SplittableRange(range);
    UTparallelFor(range, AccumulateParallel(raster, acc, width, sx, sy));
}



