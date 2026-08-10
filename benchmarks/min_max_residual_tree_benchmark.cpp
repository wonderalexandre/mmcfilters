#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/sdrt/SaturatedResidualTreeBuilder.hpp"
#include "mmcfilters/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp"
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
using namespace mmcfilters::sdrt;

namespace {

using Clock = std::chrono::steady_clock;
using SaturatedBuilder = SaturatedResidualTreeBuilder<std::uint8_t>;
using UnrestrictedBuilder = UnrestrictedResidualTreeBuilder<std::uint8_t>;

struct Options {
    enum class Mode { Both, Unrestricted, Saturated };

    std::string imagePath;
    int repetitions = 3;
    Mode mode = Mode::Both;
};

Options parseOptions(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        throw std::invalid_argument("expected image path, optional repetition count, and optional mode (both, unrestricted, or saturated)");
    }
    Options options;
    options.imagePath = argv[1];
    if (argc >= 3) {
        options.repetitions = std::stoi(argv[2]);
    }
    if (options.repetitions <= 0) {
        throw std::invalid_argument("repetitions must be positive");
    }
    if (argc == 4) {
        const std::string mode = argv[3];
        if (mode == "unrestricted") {
            options.mode = Options::Mode::Unrestricted;
        } else if (mode == "saturated") {
            options.mode = Options::Mode::Saturated;
        } else if (mode != "both") {
            throw std::invalid_argument("mode must be both, unrestricted, or saturated");
        }
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

template <class Builder> void build(Builder& builder, const ImageUInt8Ptr& image, const RegularGridAdjacency2D& adjacency) {
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    builder.build(image, std::move(minTree), std::move(maxTree));
}

template <class Function> double measureMilliseconds(Function&& function) {
    const auto start = Clock::now();
    function();
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

template <class Builder> bool reconstructs(const Builder& builder, const ImageUInt8Ptr& image) {
    for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
        const NodeId owner = builder.getProperPartOwner()[static_cast<std::size_t>(pixel)];
        if (builder.getAltitude()[static_cast<std::size_t>(owner)] != (*image)[pixel]) {
            return false;
        }
    }
    return true;
}

void printDomain(const Options& options, const ImageUInt8Ptr& image, const char* mode) {
    std::cout << "rows=" << image->getNumRows() << '\n'
              << "cols=" << image->getNumCols() << '\n'
              << "pixels=" << image->getSize() << '\n'
              << "repetitions=" << options.repetitions << '\n'
              << "mode=" << mode << '\n';
}

template <class Builder>
void runIsolated(const Options& options, const ImageUInt8Ptr& image, const RegularGridAdjacency2D& adjacency,
                 Builder& builder, const char* mode) {
    build(builder, image, adjacency);
    if (!reconstructs(builder, image)) {
        throw std::runtime_error(std::string(mode) + " min/max mode failed reconstruction");
    }
    std::vector<double> times;
    times.reserve(options.repetitions);
    for (int repetition = 0; repetition < options.repetitions; ++repetition) {
        times.push_back(measureMilliseconds([&] { build(builder, image, adjacency); }));
    }

    printDomain(options, image, mode);
    std::cout << std::fixed << std::setprecision(3)
              << mode << "_median_ms=" << median(std::move(times)) << '\n'
              << mode << "_nodes=" << builder.getNodeParent().size() << '\n'
              << mode << "_rejected_extrema=" << builder.getStatistics().rejectedExtrema << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const auto image = loadImage(options.imagePath);
        const RegularGridAdjacency2D adjacency(image->getNumRows(), image->getNumCols(), 1.0);
        if (options.mode == Options::Mode::Unrestricted) {
            UnrestrictedBuilder unrestricted(adjacency, UnrestrictedResidualTreeOptions{SdrtTiePolicy::ContrastInvariantSpatial});
            runIsolated(options, image, adjacency, unrestricted, "unrestricted");
            return 0;
        }
        if (options.mode == Options::Mode::Saturated) {
            SaturatedBuilder saturated(adjacency, NodeId{0}, SaturatedResidualTreeOptions{SdrtTiePolicy::ContrastInvariantSpatial});
            runIsolated(options, image, adjacency, saturated, "saturated");
            return 0;
        }

        UnrestrictedBuilder unrestricted(adjacency, UnrestrictedResidualTreeOptions{SdrtTiePolicy::ContrastInvariantSpatial});
        SaturatedBuilder saturated(adjacency, NodeId{0}, SaturatedResidualTreeOptions{SdrtTiePolicy::ContrastInvariantSpatial});

        build(unrestricted, image, adjacency);
        build(saturated, image, adjacency);
        if (!reconstructs(unrestricted, image) || !reconstructs(saturated, image)) {
            throw std::runtime_error("a min/max mode failed reconstruction");
        }

        std::vector<double> unrestrictedTimes;
        std::vector<double> saturatedTimes;
        unrestrictedTimes.reserve(options.repetitions);
        saturatedTimes.reserve(options.repetitions);
        for (int repetition = 0; repetition < options.repetitions; ++repetition) {
            if ((repetition & 1) == 0) {
                unrestrictedTimes.push_back(measureMilliseconds([&] { build(unrestricted, image, adjacency); }));
                saturatedTimes.push_back(measureMilliseconds([&] { build(saturated, image, adjacency); }));
            } else {
                saturatedTimes.push_back(measureMilliseconds([&] { build(saturated, image, adjacency); }));
                unrestrictedTimes.push_back(measureMilliseconds([&] { build(unrestricted, image, adjacency); }));
            }
        }

        const double unrestrictedMedian = median(unrestrictedTimes);
        const double saturatedMedian = median(saturatedTimes);
        const auto& unrestrictedStats = unrestricted.getStatistics();
        const auto& saturatedStats = saturated.getStatistics();
        printDomain(options, image, "both");
        std::cout << std::fixed << std::setprecision(3) << "unrestricted_median_ms=" << unrestrictedMedian << '\n'
                  << "saturated_median_ms=" << saturatedMedian << '\n'
                  << "saturated_over_unrestricted=" << saturatedMedian / unrestrictedMedian << '\n'
                  << "unrestricted_nodes=" << unrestricted.getNodeParent().size() << '\n'
                  << "saturated_nodes=" << saturated.getNodeParent().size() << '\n'
                  << "unrestricted_rejected_extrema=" << unrestrictedStats.rejectedExtrema << '\n'
                  << "saturated_rejected_extrema=" << saturatedStats.rejectedExtrema << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "min/max residual benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
