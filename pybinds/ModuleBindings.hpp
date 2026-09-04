#pragma once

#include <pybind11/pybind11.h>

namespace mmcfilters::pybindings {

/**
 * @brief Registers morphological tree bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initMorphologicalTree(pybind11::module_& m);
/**
 * @brief Registers attribute computation bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initAttributeComputation(pybind11::module_& m);
/**
 * @brief Registers contours computed incrementally bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initContourComputation(pybind11::module_& m);
/**
 * @brief Registers contour trace computation bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initContourTraceComputation(pybind11::module_& m);
/**
 * @brief Registers attribute filters bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initAttributeFilters(pybind11::module_& m);
/**
 * @brief Registers depth stable region computer bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initDepthStableRegionComputer(pybind11::module_& m);
/**
 * @brief Registers extinction values bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initExtinctionValues(pybind11::module_& m);
/**
 * @brief Registers regular grid adjacency2 d bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initRegularGridAdjacency2D(pybind11::module_& m);
/**
 * @brief Registers ultimate attribute opening bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initUltimateAttributeOpening(pybind11::module_& m);
/**
 * @brief Registers component tree adjust bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initComponentTreeAdjust(pybind11::module_& m);

} // namespace mmcfilters::pybindings
