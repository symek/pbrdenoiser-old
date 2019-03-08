#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <VEX/VEX_VexOp.h>

#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

#include <cstring>

#ifdef CONCURRENT_HASH_MAP
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#endif

#define MARGIN 5

namespace HA_HDK {

static std::mutex shaderstore_mutex;

// Simple vectors of rows for pixel like storage
using SampleIntStore = std::vector<std::vector<int>>;

// Map of channels
// TODO: we could keep two structures, with hashes and names
// to avoid some .find()
using ShaderStore = std::unordered_map<int, SampleIntStore>;

// Static store holding our data
static auto shaderStore = ShaderStore{};

// This function is super inefficent, but it will be optimized
// by VEX compiler and run only once.
template <VEX_Precision PREC>
static void vex_store_open(int argc, void *argv[], void *data) {
    VEXint<PREC> *result = (VEXint<PREC> *)argv[0];
    const char *channel = (const char *)argv[1];
    const VEXint<PREC> xres = *(const VEXint<PREC> *)argv[2];
    const VEXint<PREC> yres = *(const VEXint<PREC> *)argv[3];

    const size_t long_hash = std::hash<std::string>{}(channel);
    const int short_hash = static_cast<int>(long_hash);

    std::lock_guard<std::mutex> guard(shaderstore_mutex);

    if (shaderStore.find(short_hash) == shaderStore.end()) {
        shaderStore.emplace(std::make_pair<const int, SampleIntStore>(
            std::move(short_hash), SampleIntStore{}));

        auto &store = *shaderStore.find(short_hash);
        auto &samples = store.second;

        samples.resize(static_cast<size_t>(yres + MARGIN));
        for (auto &row : samples) {
            row.resize(static_cast<size_t>(xres + MARGIN));
            for (int x = 0; x < xres; ++x)
                row[x] = 0;
        }

        result[0] = short_hash;

    } else {
        const auto &samples = *shaderStore.find(short_hash);
        result[0] = samples.first;
    }
}

template <VEX_Precision PREC>
static void vex_store_increment(int argc, void *argv[], void *data) {
    auto *result = (VEXint<PREC> *)argv[0];
    const auto handle = *(const VEXint<PREC> *)argv[1];
    const auto xcoord = *(const VEXint<PREC> *)argv[2];
    const auto ycoord = *(const VEXint<PREC> *)argv[3];

    std::lock_guard<std::mutex> guard(shaderstore_mutex);
    if (shaderStore.find(handle) != shaderStore.end()) {
        auto &store = *shaderStore.find(handle);
        auto &samples = store.second;

        assert(ycoord <= samples.size());
        assert(xcoord <= samples.at(ycoord).size());

        auto &row = samples.at(ycoord);
        row.at(xcoord) += 1;
        result[0] = row.at(xcoord);
    } else {
        result[0] = -1;
    }
}

} // namespace HA_HDK

//
// Installation function
//
using namespace HA_HDK;
void newVEXOp(void *) {
    new VEX_VexOp("vex_store_open@&ISII", // Signature
                  vex_store_open<VEX_32>, // Evaluator
                  vex_store_open<VEX_64>, // Evaluator
                  VEX_ALL_CONTEXT,        // Context mask
                  nullptr, nullptr,       // init function
                  nullptr, nullptr,       // cleanup function
                  VEX_OPTIMIZE_2          // Optimization level
    );
    new VEX_VexOp("vex_store_increment@&IIII", // Signature
                  vex_store_increment<VEX_32>, // Evaluator
                  vex_store_increment<VEX_64>, // Evaluator
                  VEX_ALL_CONTEXT,             // Context mask
                  nullptr, nullptr,            // init function
                  nullptr, nullptr,            // cleanup function
                  VEX_OPTIMIZE_2               // Optimization level
    );
}
