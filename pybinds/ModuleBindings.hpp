#pragma once

#include <pybind11/pybind11.h>

namespace mmcfilters::pybindings {

void initMorphologicalTree(pybind11::module_& m);
void initAttributeComputation(pybind11::module_& m);
void initContoursComputedIncrementally(pybind11::module_& m);
void initContourTraceComputation(pybind11::module_& m);
void initAttributeFilters(pybind11::module_& m);
void initDepthStableRegionComputer(pybind11::module_& m);
void initExtinctionValues(pybind11::module_& m);
void initAdjacencyRelation(pybind11::module_& m);
void initUltimateAttributeOpening(pybind11::module_& m);
void initComponentTreeAdjust(pybind11::module_& m);

} // namespace mmcfilters::pybindings
