#pragma once

#include "../mmcfilters/filters/DepthStableRegionComputer.hpp"
#include "PybindConversions.hpp"

#include <concepts>
#include <cstdint>
#include <memory>
#include <variant>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Pybind11 wrapper exposing topological depth-variation computation.
 */
class DepthStableRegionComputerPybind {
    /** @brief References the weighted owner used by the component. */
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;
    /** @brief Stores the computer. */
    std::variant<DepthStableRegionComputer<float>, DepthStableRegionComputer<double>> computer_;

    /**
     * @brief Creates a depth-stability computer for the runtime attribute type.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @return Depth-stability computer variant for the array element type.
     */
    template <std::floating_point Real>
    static DepthStableRegionComputer<Real> makeComputer(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), weighted.topology(), "attribute");
        const py::buffer_info buffer = typed.request();
        const auto* data = static_cast<const Real*>(buffer.ptr);
        std::vector<Real> owned(data, data + buffer.shape[0]);
        return DepthStableRegionComputer<Real>(weighted.topology(), std::move(owned));
    }

    /**
     * @brief Creates a depth-stability computer for the runtime attribute type.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @return Depth-stability computer variant for the array element type.
     */
    static std::variant<DepthStableRegionComputer<float>, DepthStableRegionComputer<double>> makeComputer(WeightedMorphologicalTree<std::uint8_t>& weighted,
                                                                                                          py::array attribute) {
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return makeComputer<double>(weighted, std::move(attribute));
        }
        return makeComputer<float>(weighted, std::move(attribute));
    }

  public:
    /**
     * @brief Constructs `DepthStableRegionComputerPybind` from the supplied inputs.
     *
     * @param weighted Weighted tree used by the operation.
     */
    explicit DepthStableRegionComputerPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted)
        : weightedOwner_(std::move(weighted)), computer_(DepthStableRegionComputer<float>(weightedOwner_->topology())) {}

    /**
     * @brief Constructs `DepthStableRegionComputerPybind` from the supplied inputs.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     */
    DepthStableRegionComputerPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attribute)
        : weightedOwner_(std::move(weighted)), computer_(makeComputer(*weightedOwner_, std::move(attribute))) {}

    /**
     * @brief Computes by depth.
     *
     * @param depthDelta Topological-depth radius of the stability window.
     * @return Computed by depth.
     */
    py::array_t<std::uint8_t> computeByDepth(int depthDelta) {
        return std::visit(
            [depthDelta](auto& computer) -> py::array_t<std::uint8_t> {
                std::vector<std::uint8_t> selected = computer.computeByDepth(depthDelta);
                const int size = static_cast<int>(selected.size());
                return pybind_utils::toNumpyOwned(std::move(selected), size);
            },
            computer_);
    }

    /**
     * @brief Returns variation.
     *
     * @param node Node identifier used by the operation.
     * @return Variation.
     */
    double getVariation(NodeId node) {
        return std::visit([node](auto& computer) -> double { return static_cast<double>(computer.getVariation(node)); }, computer_);
    }

    /**
     * @brief Returns variations.
     *
     * @return Variations.
     */
    py::array getVariations() {
        return std::visit(
            [](auto& computer) -> py::array {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                const std::vector<Real>& values = computer.getVariations();
                return pybind_utils::toNumpyOwned(std::vector<Real>(values.begin(), values.end()), static_cast<int>(values.size()));
            },
            computer_);
    }

    /**
     * @brief Finds the node of minimum variation in the stability window.
     *
     * @param node Node identifier used by the operation.
     * @return Node identifier with minimum variation in the window.
     */
    NodeId nodeWithMinimumVariationInWindow(NodeId node) {
        return std::visit([node](auto& computer) -> NodeId { return computer.nodeWithMinimumVariationInWindow(node); }, computer_);
    }

    /**
     * @brief Finds the ascendant at the upper boundary of the stability window.
     *
     * @param node Node identifier used by the operation.
     * @return Node identifier at the upper stability-window boundary.
     */
    NodeId ascendantInStabilityWindow(NodeId node) const {
        return std::visit([node](const auto& computer) -> NodeId { return computer.ascendantInStabilityWindow(node); }, computer_);
    }

    /**
     * @brief Finds the descendant at the lower boundary of the stability window.
     *
     * @param node Node identifier used by the operation.
     * @return Node identifier at the lower stability-window boundary.
     */
    NodeId descendantInStabilityWindow(NodeId node) const {
        return std::visit([node](const auto& computer) -> NodeId { return computer.descendantInStabilityWindow(node); }, computer_);
    }

    /**
     * @brief Returns num nodes.
     *
     * @return Num nodes.
     */
    int getNumNodes() const {
        return std::visit([](const auto& computer) { return computer.getNumNodes(); }, computer_);
    }

    /**
     * @brief Sets max variation.
     *
     * @param value Value used by the operation.
     */
    void setMaxVariation(double value) {
        std::visit(
            [value](auto& computer) {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                computer.setMaxVariation(static_cast<Real>(value));
            },
            computer_);
    }

    /**
     * @brief Sets min attribute.
     *
     * @param value Value used by the operation.
     */
    void setMinAttribute(double value) {
        std::visit(
            [value](auto& computer) {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                computer.setMinAttribute(static_cast<Real>(value));
            },
            computer_);
    }

    /**
     * @brief Sets max attribute.
     *
     * @param value Value used by the operation.
     */
    void setMaxAttribute(double value) {
        std::visit(
            [value](auto& computer) {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                computer.setMaxAttribute(static_cast<Real>(value));
            },
            computer_);
    }
};

} // namespace mmcfilters
