#include "ModuleBindings.hpp"

#include <pybind11/pybind11.h>

/**
 * @brief Initializes the Python extension module.
 */
PYBIND11_MODULE(mmcfilters, m) {
    m.doc() = "Morphological tree filters with a NodeId-first Python API.";

    mmcfilters::pybindings::initMorphologicalTree(m);
    mmcfilters::pybindings::initAttributeComputation(m);
    mmcfilters::pybindings::initContoursComputedIncrementally(m);
    mmcfilters::pybindings::initContourTraceComputation(m);
    mmcfilters::pybindings::initAttributeFilters(m);
    mmcfilters::pybindings::initDepthStableRegionComputer(m);
    mmcfilters::pybindings::initExtinctionValues(m);
    mmcfilters::pybindings::initRegularGridAdjacency2D(m);
    mmcfilters::pybindings::initUltimateAttributeOpening(m);
    mmcfilters::pybindings::initComponentTreeAdjust(m);
}
