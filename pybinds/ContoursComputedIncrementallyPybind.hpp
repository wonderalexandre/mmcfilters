#pragma once

#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "MorphologicalTreePybind.hpp"

namespace mmcfilters {

/**
 * @brief Helper functions exposed to Python for incremental contours.
 */
class ContoursComputedIncrementallyPybind {
public:
    static ContoursComputedIncrementally::IncrementalContours extraction(MorphologicalTreePybindPtr tree) {
        return ContoursComputedIncrementally::extractCompactContours(*tree);
    }

    static ContoursComputedIncrementally::IncrementalContours extraction(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted) {
        return ContoursComputedIncrementally::extractCompactContours(weighted->asView());
    }
};

} // namespace mmcfilters
