#include "mmcfilters/trees/ExtremumResidualTreeByCTsBuilder.hpp"
#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mmcfilters;

namespace {

using Clock = std::chrono::steady_clock;
using ExtremumByCTsBuilder = mmcfilters::ExtremumResidualTreeByCTsBuilder<std::uint8_t>;

struct Options {
    std::string imagePath;
    int repetitions = 3;
};

Options parseOptions(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        throw std::invalid_argument("expected image path and optional repetition count");
    }
    Options options;
    options.imagePath = argv[1];
    if (argc >= 3) {
        options.repetitions = std::stoi(argv[2]);
    }
    if (options.repetitions <= 0) {
        throw std::invalid_argument("repetitions must be positive");
    }
    return options;
}

ImageUInt8Ptr loadImage(const std::string& path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    using Buffer = std::unique_ptr<unsigned char, decltype(&stbi_image_free)>;
    Buffer pixels(stbi_load(path.c_str(), &width, &height, &channels, 1), &stbi_image_free);
    if (!pixels || width <= 0 || height <= 0) {
        throw std::runtime_error("could not load image: " + path);
    }
    auto image = ImageUInt8::create(height, width);
    std::copy_n(pixels.get(), image->getSize(), image->rawData());
    return image;
}

template <class Function> double measureMilliseconds(Function&& function) {
    const auto start = Clock::now();
    const auto hierarchy = function();
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    if ((values.size() & 1U) != 0U) {
        return values[middle];
    }
    return (values[middle - 1] + values[middle]) / 2.0;
}

bool reconstructs(const NativeHierarchyOwned<std::uint8_t>& hierarchy, const ImageUInt8Ptr& image) {
    for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
        const NodeId smallestNode = hierarchy.smallestNodeMap[static_cast<std::size_t>(pixel)];
        if (hierarchy.nodeAltitudes[static_cast<std::size_t>(smallestNode)] != (*image)[pixel]) {
            return false;
        }
    }
    return true;
}

void printDomain(const Options& options, const ImageUInt8Ptr& image, const char* mode) {
    std::cout << "rows=" << image->getNumRows() << '\n'
              << "columns=" << image->getNumColumns() << '\n'
              << "pixels=" << image->getSize() << '\n'
              << "repetitions=" << options.repetitions << '\n'
              << "mode=" << mode << '\n';
}

template <class Builder>
void runIsolated(const Options& options, const ImageUInt8Ptr& image, Builder& builder, const char* mode) {
    const auto hierarchy = builder.build(image);
    if (!reconstructs(hierarchy, image)) {
        throw std::runtime_error(std::string(mode) + " min/max mode failed reconstruction");
    }
    std::vector<double> times;
    times.reserve(options.repetitions);
    for (int repetition = 0; repetition < options.repetitions; ++repetition) {
        times.push_back(measureMilliseconds([&] { return builder.build(image); }));
    }

    printDomain(options, image, mode);
    std::cout << std::fixed << std::setprecision(3)
              << mode << "_median_ms=" << median(std::move(times)) << '\n'
              << mode << "_nodes=" << hierarchy.parent.size() << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const auto image = loadImage(options.imagePath);
        const RegularGridAdjacency2D adjacency(image->getNumRows(), image->getNumColumns(), 1.0);
        ExtremumByCTsBuilder builder(adjacency);
        runIsolated(options, image, builder, "extremum");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "extremum residual tree by CTs benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
