/**
 * Internal profiling harness for contour trace loop materialization.
 *
 * This target compiles ContourTraceComputation with
 * MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE and reports the time spent inside
 * traceNodeLoops after edge caches have already been materialized.
 */
#define MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE

#include "mmcfilters/contours/ContourTraceComputation.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace mmcfilters;

namespace {

ImageUInt8Ptr makeBenchmarkImage(int rows, int columns) {
    auto image = ImageUInt8::create(rows, columns);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int idx = row * columns + column;
            const int radial = (row - rows / 2) * (row - rows / 2) + (column - columns / 2) * (column - columns / 2);
            (*image)[idx] = static_cast<uint8_t>((radial / 113 + row * 7 + column * 3) & 0xff);
        }
    }
    return image;
}

ImageUInt8Ptr loadGrayscaleImage(const std::filesystem::path& path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char* raw = stbi_load(path.string().c_str(), &width, &height, &channels, 1);
    if (!raw) {
        throw std::runtime_error("Failed to load image: " + path.string());
    }

    auto image = ImageUInt8::create(height, width);
    const int size = width * height;
    for (int i = 0; i < size; ++i) {
        (*image)[i] = raw[i];
    }
    stbi_image_free(raw);
    return image;
}

template <typename Fn> auto timed(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    auto result = fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return std::pair{std::move(result), static_cast<double>(micros) / 1000.0};
}

template <typename Fn> double timedVoid(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return static_cast<double>(micros) / 1000.0;
}

double nsToMs(std::int64_t ns) { return static_cast<double>(ns) / 1'000'000.0; }

double safeRatio(std::size_t numerator, std::size_t denominator) {
    if (denominator == 0) {
        return 0.0;
    }
    return static_cast<double>(numerator) / static_cast<double>(denominator);
}

using TraceProfileStats = ContourTraceComputation::IncrementalContourTraces::TraceProfileStats;

std::int64_t accountedNs(const TraceProfileStats& stats) {
    return stats.profileCountersNs + stats.buildAdjacencyNs + stats.walkLoopsNs + stats.resetOutgoingNs + stats.commitEdgesNs + stats.commitLoopsNs +
           stats.releaseScratchNs;
}

void printProfileStats(const TraceProfileStats& stats) {
    const auto accounted = accountedNs(stats);
    std::cout << "    nodes_traced=" << stats.nodesTraced << " edges_traced=" << stats.edgesTraced << " loops_traced=" << stats.loopsTraced << '\n';
    std::cout << "    outgoing_vertices=" << stats.outgoingVertices << " single_outgoing_vertices=" << stats.singleOutgoingVertices
              << " multi_outgoing_vertices=" << stats.multiOutgoingVertices << " avg_outgoing_degree=" << safeRatio(stats.edgesTraced, stats.outgoingVertices)
              << " max_outgoing_degree=" << stats.maxOutgoingDegree << '\n';
    std::cout << "    successor_steps_single=" << stats.singleSuccessorSteps << " successor_steps_ambiguous=" << stats.ambiguousSuccessorSteps
              << " avg_scan_per_ambiguous=" << safeRatio(stats.successorCandidateScans, stats.ambiguousSuccessorSteps)
              << " visited_scan_fraction=" << safeRatio(stats.successorVisitedSkips, stats.successorCandidateScans) << '\n';
    std::cout << "    successor_scans=" << stats.successorCandidateScans << " unvisited_candidates=" << stats.successorUnvisitedCandidates
              << " visited_skips=" << stats.successorVisitedSkips << " single_visited_stops=" << stats.singleSuccessorVisitedStops
              << " ambiguous_dead_ends=" << stats.ambiguousSuccessorDeadEnds << '\n';
    std::cout << "    closed_loop_stops=" << stats.closedLoopStops << " missing_outgoing_stops=" << stats.missingOutgoingStops << '\n';
    std::cout << "    profile_counters=" << nsToMs(stats.profileCountersNs) << " ms"
              << " build_adjacency=" << nsToMs(stats.buildAdjacencyNs) << " ms"
              << " walk_loops=" << nsToMs(stats.walkLoopsNs) << " ms"
              << " reset_outgoing=" << nsToMs(stats.resetOutgoingNs) << " ms\n";
    std::cout << "    commit_edges=" << nsToMs(stats.commitEdgesNs) << " ms"
              << " commit_loops=" << nsToMs(stats.commitLoopsNs) << " ms"
              << " release_scratch=" << nsToMs(stats.releaseScratchNs) << " ms"
              << " accounted=" << nsToMs(accounted) << " ms\n";
}

void profileCase(const std::string& label, const MorphologicalTree& tree, int repeats) {
    double totalExtractMs = 0.0;
    double totalEdgeMs = 0.0;
    double totalLoopWallMs = 0.0;
    TraceProfileStats totalStats;

    for (int repeat = 0; repeat < repeats; ++repeat) {
        auto [traces, extractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
        totalExtractMs += extractMs;

        totalEdgeMs += timedVoid([&]() {
            const auto edges = traces.getEdges(tree.root());
            volatile std::size_t edgeCount = edges.size();
            (void)edgeCount;
        });

        totalLoopWallMs += timedVoid([&]() {
            const TraceProfileStats stats = traces.profileMaterializeAllLoops();
            totalStats.nodesTraced += stats.nodesTraced;
            totalStats.edgesTraced += stats.edgesTraced;
            totalStats.loopsTraced += stats.loopsTraced;
            totalStats.outgoingVertices += stats.outgoingVertices;
            totalStats.singleOutgoingVertices += stats.singleOutgoingVertices;
            totalStats.multiOutgoingVertices += stats.multiOutgoingVertices;
            totalStats.maxOutgoingDegree = std::max(totalStats.maxOutgoingDegree, stats.maxOutgoingDegree);
            totalStats.closedLoopStops += stats.closedLoopStops;
            totalStats.missingOutgoingStops += stats.missingOutgoingStops;
            totalStats.singleSuccessorSteps += stats.singleSuccessorSteps;
            totalStats.singleSuccessorVisitedStops += stats.singleSuccessorVisitedStops;
            totalStats.ambiguousSuccessorSteps += stats.ambiguousSuccessorSteps;
            totalStats.ambiguousSuccessorDeadEnds += stats.ambiguousSuccessorDeadEnds;
            totalStats.successorCandidateScans += stats.successorCandidateScans;
            totalStats.successorVisitedSkips += stats.successorVisitedSkips;
            totalStats.successorUnvisitedCandidates += stats.successorUnvisitedCandidates;
            totalStats.profileCountersNs += stats.profileCountersNs;
            totalStats.buildAdjacencyNs += stats.buildAdjacencyNs;
            totalStats.walkLoopsNs += stats.walkLoopsNs;
            totalStats.resetOutgoingNs += stats.resetOutgoingNs;
            totalStats.commitEdgesNs += stats.commitEdgesNs;
            totalStats.commitLoopsNs += stats.commitLoopsNs;
            totalStats.releaseScratchNs += stats.releaseScratchNs;
        });
    }

    std::cout << "\n"
              << label << '\n'
              << "  nodes: live=" << tree.numNodes() << " slots=" << tree.numInternalNodeSlots() << " proper_parts=" << tree.numPixels()
              << '\n'
              << "  extract_avg=" << totalExtractMs / repeats << " ms"
              << " edge_materialize_avg=" << totalEdgeMs / repeats << " ms"
              << " loop_trace_wall_avg=" << totalLoopWallMs / repeats << " ms\n";

    TraceProfileStats avgStats = totalStats;
    avgStats.nodesTraced /= static_cast<std::size_t>(repeats);
    avgStats.edgesTraced /= static_cast<std::size_t>(repeats);
    avgStats.loopsTraced /= static_cast<std::size_t>(repeats);
    avgStats.outgoingVertices /= static_cast<std::size_t>(repeats);
    avgStats.singleOutgoingVertices /= static_cast<std::size_t>(repeats);
    avgStats.multiOutgoingVertices /= static_cast<std::size_t>(repeats);
    avgStats.closedLoopStops /= static_cast<std::size_t>(repeats);
    avgStats.missingOutgoingStops /= static_cast<std::size_t>(repeats);
    avgStats.singleSuccessorSteps /= static_cast<std::size_t>(repeats);
    avgStats.singleSuccessorVisitedStops /= static_cast<std::size_t>(repeats);
    avgStats.ambiguousSuccessorSteps /= static_cast<std::size_t>(repeats);
    avgStats.ambiguousSuccessorDeadEnds /= static_cast<std::size_t>(repeats);
    avgStats.successorCandidateScans /= static_cast<std::size_t>(repeats);
    avgStats.successorVisitedSkips /= static_cast<std::size_t>(repeats);
    avgStats.successorUnvisitedCandidates /= static_cast<std::size_t>(repeats);
    avgStats.profileCountersNs /= repeats;
    avgStats.buildAdjacencyNs /= repeats;
    avgStats.walkLoopsNs /= repeats;
    avgStats.resetOutgoingNs /= repeats;
    avgStats.commitEdgesNs /= repeats;
    avgStats.commitLoopsNs /= repeats;
    avgStats.releaseScratchNs /= repeats;
    printProfileStats(avgStats);
}

void runProfile(ImageUInt8Ptr image, int repeats) {
    auto [maxTree, maxBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMaxTree(image, 1.5); });
    std::cout << "max-tree build: " << maxBuildMs << " ms\n";
    profileCase("max-tree contour trace profile", maxTree.topology(), repeats);

    auto [minTree, minBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMinTree(image, 1.5); });
    std::cout << "\nmin-tree build: " << minBuildMs << " ms\n";
    profileCase("min-tree contour trace profile", minTree.topology(), repeats);
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4) {
            const int rows = std::stoi(argv[1]);
            const int columns = std::stoi(argv[2]);
            const int repeats = std::stoi(argv[3]);
            if (rows <= 0 || columns <= 0 || repeats <= 0) {
                throw std::invalid_argument("rows, columns, and repeats must be positive.");
            }
            runProfile(makeBenchmarkImage(rows, columns), repeats);
            return 0;
        }

        if (argc == 3) {
            const int repeats = std::stoi(argv[2]);
            if (repeats <= 0) {
                throw std::invalid_argument("repeats must be positive.");
            }
            runProfile(loadGrayscaleImage(argv[1]), repeats);
            return 0;
        }

        std::cerr << "usage:\n"
                  << "  " << argv[0] << " rows columns repeats\n"
                  << "  " << argv[0] << " image-path repeats\n";
        return 2;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return 1;
    }
}
