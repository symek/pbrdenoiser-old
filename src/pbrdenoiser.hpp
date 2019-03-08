#pragma once
#ifdef BUILD_WITH_OIDN
#include <OpenImageDenoise/oidn.hpp>
#endif
#include <atomic>
#include <chrono>
#include <ostream>

namespace pbrd {

bool save_rasters_to_file(const char *filename, const IMG_Stat &stat,
                          UT_PtrArray<PXL_Raster *> &raster_array) {

    IMG_FileParms parms = IMG_FileParms();
    IMG_Stat output_stat = IMG_Stat(stat);
    output_stat.setFilename(filename);
    parms.setDataType(IMG_HALF); // Back to half when possible.
    IMG_File *output_file =
        IMG_File::create(filename, (const IMG_Stat)output_stat, &parms);

    if (output_file) {
        output_file->writeImages(raster_array);
        output_file->close();
        return true;
    }

    return false;
}
std::unique_ptr<IMG_File>
read_rasters_as_float(const char *filename, UT_PtrArray<PXL_Raster *> raster_array)
 {
    IMG_FileParms input_parms = IMG_FileParms();
    input_parms.setDataType(IMG_FLOAT);
    input_parms.readAlphaAsPlane();
    std::unique_ptr<IMG_File> input_file;
    auto tmp = IMG_File::open(filename, &input_parms);
    input_file.reset(tmp);

    if (!input_file)
        return nullptr;

    if (!input_file->readImages(raster_array)) {
        return nullptr;
    }

    return input_file;
}

PXL_Raster *get_raster_by_name(const char *plane, const IMG_File *input_file,
                               UT_PtrArray<PXL_Raster *> raster_array) {
    const IMG_Stat &input_stat = input_file->getStat();
    const int raster_index = input_stat.getPlaneIndex(plane);
    if (raster_index == -1) {
        return nullptr;
    }
    return raster_array(raster_index);
}

/// https://codereview.stackexchange.com/questions/196245/extremely-simple-timer-class-in-c
template <typename Clock = std::chrono::high_resolution_clock> class Timer {
    typename Clock::time_point start_point;

  public:
    using Info = std::vector<std::string>;
    Timer() : start_point(Clock::now()) {}

    template <typename Rep = typename Clock::duration::rep,
              typename Units = typename Clock::duration>
    Rep elapsed_time() const {
        std::atomic_thread_fence(std::memory_order_relaxed);
        auto counted_time =
            std::chrono::duration_cast<Units>(Clock::now() - start_point)
                .count();
        std::atomic_thread_fence(std::memory_order_relaxed);
        return static_cast<Rep>(counted_time);
    }
    void restart() { start_point = Clock::now(); }

    void begin(const char *info) { print(info, false); }

    void begin(const Info &info) { print(info, false); }

    void end(const char *info = nullptr) {
        if (!info)
            print("done.", true);
        else
            print(info, true);
    }

    void print(const char *info, const bool timeit = false) {
        print(Info{info}, timeit);
    }

    void print(const Info &info, const bool timeit = false) {
        print_info(info, timeit);
        restart();
    }

  private:
    void print_info(const Info &info, const bool timeit = false) {
        if (!timeit)
            std::cout << "INFO: ";

        for (const auto &s : info)
            std::cout << s;
        std::cout << std::flush;

        if (timeit)
            std::cout
                << " ("
                << elapsed_time<unsigned int, std::chrono::milliseconds>() /
                       1000.0
                << "s)\n"
                << std::flush;
    }
};

using precise_timer = Timer<>;
using system_timer = Timer<std::chrono::system_clock>;
using monotonic_timer = Timer<std::chrono::steady_clock>;

namespace intel {
    std::unique_ptr<float[]> filter_with_oidn(
        void *color, void *albedo, void *normals, const char **error,
        const int xres, const int yres, const bool hdri, const bool srgb) {

        auto output = std::make_unique<float[]>(3 * xres * yres);

        // Create an Open Image Denoise device
        oidn::DeviceRef device = oidn::newDevice();
        device.commit();

        // Create a denoising filter
        oidn::FilterRef filter =
            device.newFilter("RT"); // generic ray tracing filter
        if (albedo) {
            filter.setImage("albedo", albedo, oidn::Format::Float3, xres,
                            yres); // optional
        }
        if (normals) {
            filter.setImage("normal", normals, oidn::Format::Float3, xres,
                            yres); // optional
        }

        filter.setImage("color", color, oidn::Format::Float3, xres, yres);
        filter.setImage("output", (void *)output.get(), oidn::Format::Float3,
                        xres, yres);
        filter.set("hdr", hdri);  // image is HDR
        filter.set("srgb", srgb); // image is srgb
        filter.commit();
        // Filter the image
        filter.execute();

        // Check for errors
        if (device.getError(*error) != oidn::Error::None) {
            return nullptr;
        }

        return output;
    }

} // namespace intel

} // namespace pbrd