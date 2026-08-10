#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/utils/Contract.hpp"
#include "mmcfilters/utils/Image.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

mmcfilters::ImageUInt8Ptr makeInput(int rows, int cols) {
    auto image = mmcfilters::ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int radial = (row - rows / 2) * (row - rows / 2) + (col - cols / 2) * (col - cols / 2);
            (*image)[row * cols + col] = static_cast<std::uint8_t>((radial + 17 * row + 31 * col) & 0xff);
        }
    }
    return image;
}

} // namespace

int main(int argc, char** argv) {
    const int rows = argc > 1 ? std::atoi(argv[1]) : 128;
    const int cols = argc > 2 ? std::atoi(argv[2]) : 128;
    const int repetitions = argc > 3 ? std::atoi(argv[3]) : 7;
    if (rows <= 0 || cols <= 0 || repetitions <= 0) {
        std::cerr << "usage: scientific_pipeline_benchmark [rows>0] [cols>0] [repetitions>0]\n";
        return EXIT_FAILURE;
    }

    const auto image = makeInput(rows, cols);
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));
    std::uint64_t checksum = 0;

    for (int repetition = 0; repetition < repetitions; ++repetition) {
        const auto start = Clock::now();
        const auto weighted = mmcfilters::MorphologicalTreeFactory::createMaxTree(image, 1.5);
        const auto maxDist = mmcfilters::AttributeComputation::computeSingleAttribute<double>(weighted, mmcfilters::MAX_DIST);

        std::vector<bool> keep(static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()), false);
        for (mmcfilters::NodeId node = 0; node < weighted.topology().getNumInternalNodeSlots(); ++node) {
            keep[static_cast<std::size_t>(node)] = maxDist.second[static_cast<std::size_t>(node)] >= 1.0;
        }
        keep[static_cast<std::size_t>(weighted.topology().getRoot())] = true;

        mmcfilters::AttributeFilters<std::uint8_t> filters(weighted);
        const auto filtered = filters.filteringByDirectRule(keep);
        for (int pixel = 0; pixel < filtered->getSize(); ++pixel) {
            checksum += (*filtered)[pixel];
        }
        const auto finish = Clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(finish - start).count());
    }

    std::sort(samples.begin(), samples.end());
    std::cout << std::fixed << std::setprecision(3) << "contract_mode=" << (mmcfilters::contract::validationsEnabled ? "CHECKED" : "UNCHECKED") << '\n'
              << "rows=" << rows << '\n'
              << "cols=" << cols << '\n'
              << "repetitions=" << repetitions << '\n'
              << "median_pipeline_ms=" << samples[samples.size() / 2] << '\n'
              << "checksum=" << checksum << '\n';
    return EXIT_SUCCESS;
}
