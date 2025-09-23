#pragma once

#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"

#include "MorphologicalTreePybind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>
namespace mmcfilters {

class ContoursComputedIncrementallyPybind;
using ContoursComputedIncrementallyPybindPtr = std::shared_ptr<ContoursComputedIncrementallyPybind>;

/**
 * @brief Camada Pybind que expõe utilitários incrementais de contornos.
 */
class ContoursComputedIncrementallyPybind : public ContoursComputedIncrementally{

    public:
    using ContoursComputedIncrementally::ContoursComputedIncrementally;

	
	static py::dict extractContours(MorphologicalTreePybindPtr tree) {
		auto contours = ContoursComputedIncrementally::extractNonCompactContours(tree);  // chama o método original
	
		py::dict pyContours;
		for (size_t nodeIdx = 0; nodeIdx < contours.size(); ++nodeIdx) {
			py::set pySet;
			for (int pixel : contours[nodeIdx]) {
				pySet.add(pixel);
			}
			pyContours[py::int_(nodeIdx)] = pySet;
		}
	
		return pyContours;
	}


	static std::shared_ptr<Contours> extractCompactContours(MorphologicalTreePybindPtr tree){
			return ContoursComputedIncrementally::extractCompactContours(tree);
	}


};

} // namespace mmcfilters
