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
inline const char* altitudeOrderName(AltitudeOrder order) noexcept {
    switch (order) {
    case AltitudeOrder::INCREASING_FROM_ROOT:
        return "INCREASING_FROM_ROOT";
    case AltitudeOrder::DECREASING_FROM_ROOT:
        return "DECREASING_FROM_ROOT";
    case AltitudeOrder::UNCONSTRAINED:
        return "UNCONSTRAINED";
    }
    return "UNKNOWN";
}

/**
 * @brief Tests whether global monotone altitude order holds.
 *
 * @param tree Tree topology used by the operation.
 * @return True when global monotone altitude order; otherwise false.
 */
inline bool hasGlobalMonotoneAltitudeOrder(const MorphologicalTree& tree) noexcept { return tree.getAltitudeOrder() != AltitudeOrder::UNCONSTRAINED; }

/**
 * @brief Validates global monotone altitude order.
 *
 * @param tree Tree topology used by the operation.
 * @param context Operation name used in diagnostics.
 */
inline void validateGlobalMonotoneAltitudeOrder(const MorphologicalTree& tree, const char* context) {
    if (!hasGlobalMonotoneAltitudeOrder(tree)) {
        throw std::invalid_argument(std::string(context) + " requires a globally monotone altitude order; got " + altitudeOrderName(tree.getAltitudeOrder()) +
                                    ".");
    }
}

} // namespace mmcfilters::detail
