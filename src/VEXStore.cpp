#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <VEX/VEX_VexOp.h>

#include <functional>
#include <memory>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <iostream>

#include <cstring>

#ifdef CONCURRENT_HASH_MAP
#include <tbb/concurrent_vector.h> 
#include <tbb/concurrent_hash_map.h>
#endif

#define MARGIN 5

namespace HA_HDK {

static std::mutex shaderstore_mutex;

using  SampleIntStore   = std::vector<std::vector<int>>;
using  ShaderStore      = std::unordered_map<int, SampleIntStore>;

static auto shaderStore = ShaderStore{};

template <VEX_Precision PREC>
static void vex_store_open(int argc, void *argv[], void *data) {
          VEXint<PREC> *result = (VEXint<PREC>*)       argv[0];
    const char        *channel = (const char*)         argv[1];
    const VEXint<PREC>  xres  = *(const VEXint<PREC>*) argv[2];
    const VEXint<PREC>  yres  = *(const VEXint<PREC>*) argv[3];

    const size_t long_hash  = std::hash<std::string>{}(channel);
    const int    short_hash = static_cast<int>(long_hash);

    std::lock_guard<std::mutex> guard(shaderstore_mutex);
    result[0] = short_hash;

    if (shaderStore.find(short_hash) == shaderStore.end()) {
        shaderStore.emplace(std::make_pair<int, SampleIntStore>(
            std::move((int)short_hash), SampleIntStore{}));

        auto & store = *shaderStore.find(short_hash);
        auto & samples = store.second;

        samples.resize(yres+MARGIN);
        for (auto & row: samples) {
            row.resize(xres+MARGIN);
            for (int x=0; x<xres; ++x)
                row[x] = 0;
        }

        result[0] = short_hash; 

    } else {
        const auto & samples = *shaderStore.find(short_hash);
        result[0] = samples.first; 
    }
}

template <VEX_Precision PREC>
static void vex_store_increament(int argc, void *argv[], void *data) {
          VEXint<PREC>   *result = (      VEXint<PREC>* )  argv[0];
    const VEXint<PREC>   handle = *(const VEXint<PREC>* )  argv[1];
    const VEXint<PREC>   x      = *(const VEXint<PREC>* )  argv[2];
    const VEXint<PREC>   y      = *(const VEXint<PREC>* )  argv[3];

    std::lock_guard<std::mutex> guard(shaderstore_mutex);
    if (shaderStore.find(handle) != shaderStore.end()) {
        auto & store = *shaderStore.find(handle);
        auto & samples = store.second;
        auto & row = samples.at(y);
        row.at(x) += 1;
        result[0] = row.at(x); 
    } else {
        result[0] = -1;    
    }
}

}// end of HA_HDK namespace

//
// Installation function
//
using namespace HA_HDK;
void
newVEXOp(void *)
{
     new VEX_VexOp("vex_store_open@&ISII",  // Signature
        vex_store_open<VEX_32>,      // Evaluator
        vex_store_open<VEX_64>,      // Evaluator
        VEX_ALL_CONTEXT,     // Context mask
        nullptr,nullptr,      // init function
        nullptr,nullptr,      // cleanup function
        VEX_OPTIMIZE_2       // Optimization level
        );
    new VEX_VexOp("vex_store_increament@&IIII",  // Signature
        vex_store_increament<VEX_32>,  // Evaluator
        vex_store_increament<VEX_64>,  // Evaluator
        VEX_ALL_CONTEXT,     // Context mask
        nullptr,nullptr,      // init function
        nullptr,nullptr,    // cleanup function
        VEX_OPTIMIZE_2       // Optimization level
        );         
}

