#pragma once

#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "MorphologicalTreePybind.hpp"

namespace mmcfilters {

/**
 * @brief Funções auxiliares expostas ao Python para contornos incrementais.
 */
class ContoursComputedIncrementallyPybind {
public:
    static ContoursComputedIncrementally::IncrementalContours extraction(MorphologicalTreePybindPtr tree) {
        return ContoursComputedIncrementally::extractCompactContours(*tree);
    }

    static ContoursComputedIncrementally::IncrementalContours extraction(std::shared_ptr<WeightedMorphologicalTree> weighted) {
        return ContoursComputedIncrementally::extractCompactContours(weighted->topology());
    }
};

} // namespace mmcfilters
