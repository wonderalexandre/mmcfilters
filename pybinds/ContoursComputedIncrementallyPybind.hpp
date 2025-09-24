#pragma once

#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "MorphologicalTreePybind.hpp"

namespace mmcfilters {

/**
 * @brief Funções auxiliares expostas ao Python para contornos incrementais.
 */
class ContoursComputedIncrementallyPybind {
public:
    static ContoursComputedIncrementally::IncrementalContours extraction(MorphologicalTreePybindPtr tree) {
        return ContoursComputedIncrementally::extractCompactContours(tree.get());
    }
};

} // namespace mmcfilters