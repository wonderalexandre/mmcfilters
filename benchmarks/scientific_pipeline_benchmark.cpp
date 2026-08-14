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

mmcfilters::ImageUInt8Ptr makeInput(int rows, int columns) {
    auto image = mmcfilters::ImageUInt8::create(rows, columns);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int radial = (row - rows / 2) * (row - rows / 2) + (column - columns / 2) * (column - columns / 2);
            (*image)[row * columns + column] = static_cast<std::uint8_t>((radial + 17 * row + 31 * column) & 0xff);
        }
    }
    return image;
}

} // namespace

int main(int argc, char** argv) {
    const int rows = argc > 1 ? std::atoi(argv[1]) : 128;
    const int columns = argc > 2 ? std::atoi(argv[2]) : 128;
    const int repetitions = argc > 3 ? std::atoi(argv[3]) : 7;
    if (rows <= 0 || columns <= 0 || repetitions <= 0) {
        std::cerr << "usage: scientific_pipeline_benchmark [rows>0] [columns>0] [repetitions>0]\n";
        return EXIT_FAILURE;
    }

    const auto image = makeInput(rows, columns);
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));
    std::uint64_t checksum = 0;

    for (int repetition = 0; repetition < repetitions; ++repetition) {
        const auto start = Clock::now();
        const auto valuedTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(image, 1.5);
        const auto maxDist = mmcfilters::AttributeComputation::computeSingleAttribute<double>(valuedTree, mmcfilters::MaxDist);

        std::vector<bool> keep(static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()), false);
        for (mmcfilters::NodeId node = 0; node < valuedTree.topology().numInternalNodeSlots(); ++node) {
            keep[static_cast<std::size_t>(node)] = maxDist.second[static_cast<std::size_t>(node)] >= 1.0;
        }
        keep[static_cast<std::size_t>(valuedTree.topology().root())] = true;

        mmcfilters::DirectAttributeFilter<std::uint8_t> filter(valuedTree);
        const auto filtered = filter.applyDirectAttributeFilter(mmcfilters::NodePreservationMask(std::move(keep)));
        for (mmcfilters::PixelId pixel = 0; pixel < filtered->getSize(); ++pixel) {
            checksum += (*filtered)[pixel];
        }
        const auto finish = Clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(finish - start).count());
    }

    std::sort(samples.begin(), samples.end());
    std::cout << std::fixed << std::setprecision(3) << "contract_mode=" << (mmcfilters::contract::validationsEnabled ? "CHECKED" : "UNCHECKED") << '\n'
              << "rows=" << rows << '\n'
              << "columns=" << columns << '\n'
              << "repetitions=" << repetitions << '\n'
              << "median_pipeline_ms=" << samples[samples.size() / 2] << '\n'
              << "checksum=" << checksum << '\n';
    return EXIT_SUCCESS;
}
