#pragma once

#include "../mmcfilters/filters/DepthStableRegionComputer.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

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
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;
    std::variant<DepthStableRegionComputer<float>, DepthStableRegionComputer<double>> computer_;

    template <std::floating_point Real>
    static DepthStableRegionComputer<Real> makeComputer(
        WeightedMorphologicalTree<std::uint8_t>& weighted,
        py::array attribute) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(
            std::move(attribute),
            weighted.topology(),
            "attribute");
        const py::buffer_info buffer = typed.request();
        const auto* data = static_cast<const Real*>(buffer.ptr);
        std::vector<Real> owned(data, data + buffer.shape[0]);
        return DepthStableRegionComputer<Real>(weighted.topology(), std::move(owned));
    }

    static std::variant<DepthStableRegionComputer<float>, DepthStableRegionComputer<double>> makeComputer(
        WeightedMorphologicalTree<std::uint8_t>& weighted,
        py::array attribute) {
        if (PybindUtils::parseFloatingArrayDType(attribute, "attribute") == PybindUtils::FloatingDType::Float64) {
            return makeComputer<double>(weighted, std::move(attribute));
        }
        return makeComputer<float>(weighted, std::move(attribute));
    }

public:
    explicit DepthStableRegionComputerPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted)
        : weightedOwner_(std::move(weighted)),
          computer_(DepthStableRegionComputer<float>(weightedOwner_->topology())) {}

    DepthStableRegionComputerPybind(
        std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted,
        py::array attribute)
        : weightedOwner_(std::move(weighted)),
          computer_(makeComputer(*weightedOwner_, std::move(attribute))) {}

    py::array_t<std::uint8_t> computeByDepth(int depthDelta) {
        return std::visit(
            [depthDelta](auto& computer) -> py::array_t<std::uint8_t> {
                std::vector<std::uint8_t> selected = computer.computeByDepth(depthDelta);
                const int size = static_cast<int>(selected.size());
                return PybindUtils::toNumpyOwned(std::move(selected), size);
            },
            computer_);
    }

    double getVariation(NodeId node) {
        return std::visit(
            [node](auto& computer) -> double {
                return static_cast<double>(computer.getVariation(node));
            },
            computer_);
    }

    py::array getVariations() {
        return std::visit(
            [](auto& computer) -> py::array {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                const std::vector<Real>& values = computer.getVariations();
                return PybindUtils::toNumpyOwned(std::vector<Real>(values.begin(), values.end()), static_cast<int>(values.size()));
            },
            computer_);
    }

    NodeId nodeWithMinimumVariationInWindow(NodeId node) {
        return std::visit(
            [node](auto& computer) -> NodeId {
                return computer.nodeWithMinimumVariationInWindow(node);
            },
            computer_);
    }

    NodeId ascendantInStabilityWindow(NodeId node) const {
        return std::visit(
            [node](const auto& computer) -> NodeId {
                return computer.ascendantInStabilityWindow(node);
            },
            computer_);
    }

    NodeId descendantInStabilityWindow(NodeId node) const {
        return std::visit(
            [node](const auto& computer) -> NodeId {
                return computer.descendantInStabilityWindow(node);
            },
            computer_);
    }

    int getNumNodes() const {
        return std::visit(
            [](const auto& computer) {
                return computer.getNumNodes();
            },
            computer_);
    }

    void setMaxVariation(double value) {
        std::visit(
            [value](auto& computer) {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                computer.setMaxVariation(static_cast<Real>(value));
            },
            computer_);
    }

    void setMinAttribute(double value) {
        std::visit(
            [value](auto& computer) {
                using Computer = std::decay_t<decltype(computer)>;
                using Real = typename Computer::variation_value_type;
                computer.setMinAttribute(static_cast<Real>(value));
            },
            computer_);
    }

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
