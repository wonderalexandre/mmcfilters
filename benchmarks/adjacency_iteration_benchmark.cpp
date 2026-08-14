#include "mmcfilters/utils/RegularGridAdjacency2D.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

volatile std::uint64_t adjacencyBenchmarkSink = 0;

template <class Callable> double medianMilliseconds(int repetitions, Callable&& callable) {
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        const auto start = Clock::now();
        callable();
        const auto finish = Clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(finish - start).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

template <class RangeFactory>
double benchmarkTraversal(const mmcfilters::RegularGridAdjacency2D& adjacency, int numPixels, int repetitions, RangeFactory&& rangeFactory) {
    return medianMilliseconds(repetitions, [&] {
        std::uint64_t checksum = 0;
        for (mmcfilters::PixelId pixel = 0; pixel < numPixels; ++pixel) {
            for (int neighbor : rangeFactory(adjacency, pixel)) {
                checksum += static_cast<std::uint64_t>(neighbor + 1);
            }
        }
        adjacencyBenchmarkSink = adjacencyBenchmarkSink + checksum;
    });
}

} // namespace

int main(int argc, char** argv) {
    const int rows = argc > 1 ? std::atoi(argv[1]) : 512;
    const int columns = argc > 2 ? std::atoi(argv[2]) : 512;
    const int repetitions = argc > 3 ? std::atoi(argv[3]) : 21;
    if (rows <= 0 || columns <= 0 || repetitions <= 0) {
        std::cerr << "usage: adjacency_iteration_benchmark "
                     "[rows>0] [columns>0] [repetitions>0]\n";
        return EXIT_FAILURE;
    }

    mmcfilters::RegularGridAdjacency2D adjacency(rows, columns, 1.5);
    const int numPixels = rows * columns;
    const double fullMilliseconds =
        benchmarkTraversal(adjacency, numPixels, repetitions,
                           [](const mmcfilters::RegularGridAdjacency2D& relation, mmcfilters::PixelId pixel) -> decltype(auto) {
                               return relation.getNeighborIndices(pixel);
                           });
    const double forwardMilliseconds =
        benchmarkTraversal(adjacency, numPixels, repetitions, [](const mmcfilters::RegularGridAdjacency2D& relation, mmcfilters::PixelId pixel) -> decltype(auto) {
            return relation.getForwardNeighborIndices(pixel);
        });

    std::cout << std::fixed << std::setprecision(3) << "rows=" << rows << '\n'
              << "columns=" << columns << '\n'
              << "full_ms=" << fullMilliseconds << '\n'
              << "forward_ms=" << forwardMilliseconds << '\n'
              << "sink=" << adjacencyBenchmarkSink << '\n';
    return EXIT_SUCCESS;
}
