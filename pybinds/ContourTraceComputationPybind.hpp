#pragma once

#include "../mmcfilters/contours/ContourTraceComputation.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "MorphologicalTreePybind.hpp"

namespace mmcfilters {

/**
 * @brief Helper functions exposed to Python for contour trace extraction.
 */
class ContourTraceComputationPybind {
public:
    static ContourTraceComputation::IncrementalContourTraces extraction(MorphologicalTreePybindPtr tree) {
        return ContourTraceComputation::extract(*tree);
    }

    static ContourTraceComputation::IncrementalContourTraces extraction(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted) {
        return ContourTraceComputation::extract(weighted->asView());
    }
};

} // namespace mmcfilters
