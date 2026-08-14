#pragma once

#include "../MorphologicalTree.hpp"

#include <stdexcept>
#include <string>

namespace mmcfilters::detail {

/**
 * @brief Returns order name.
 *
 * @param order Requested hierarchy ordering convention.
 * @return Order name.
 */
inline const char* nodeAltitudeOrderName(NodeAltitudeOrder order) noexcept {
    switch (order) {
    case NodeAltitudeOrder::Increasing:
        return "Increasing";
    case NodeAltitudeOrder::Decreasing:
        return "Decreasing";
    case NodeAltitudeOrder::Unconstrained:
        return "Unconstrained";
    }
    return "UNKNOWN";
}

/**
 * @brief Tests whether global monotone altitude order holds.
 *
 * @param tree Tree topology.
 * @return True when global monotone altitude order; otherwise false.
 */
inline bool hasGlobalMonotoneAltitudeOrder(const MorphologicalTree& tree) noexcept { return tree.nodeAltitudeOrder() != NodeAltitudeOrder::Unconstrained; }

/**
 * @brief Validates global monotone altitude order.
 *
 * @param tree Tree topology.
 * @param context Operation name used in diagnostics.
 */
inline void validateGlobalMonotoneAltitudeOrder(const MorphologicalTree& tree, const char* context) {
    if (!hasGlobalMonotoneAltitudeOrder(tree)) {
        throw std::invalid_argument(std::string(context) + " requires a globally monotone altitude order; got " + nodeAltitudeOrderName(tree.nodeAltitudeOrder()) +
                                    ".");
    }
}

} // namespace mmcfilters::detail
