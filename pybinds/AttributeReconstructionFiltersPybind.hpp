#pragma once

#include "../mmcfilters/filters/AttributeReconstructionFilters.hpp"
#include "PybindConversions.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include <concepts>
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace mmcfilters::pybindings {

namespace py = pybind11;

/**
 * @brief Validates and retains the Python-owned valued tree used by a reconstruction filter.
 * @param valuedTree Shared Python valued-tree wrapper.
 * @param context Public operation name used in validation errors.
 * @return The validated non-null wrapper.
 */
inline std::shared_ptr<PythonValuedMorphologicalTree> requireReconstructionFilterOwner(
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, const char* context) {
    if (!valuedTree) {
        throw std::invalid_argument(std::string(context) + " requires a non-null ValuedMorphologicalTree.");
    }
    return valuedTree;
}

/** @brief Type-erased Python wrapper for direct attribute reconstruction. */
class DirectAttributeFilterPybind {
  private:
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree_; ///< Python owner retained for the wrapper lifetime.
    std::size_t topologyMutationVersion_;                      ///< Topology version captured during construction.

  public:
    /**
     * @brief Creates a direct reconstruction wrapper.
     * @param valuedTree Python-owned valued tree used by subsequent reconstructions.
     */
    explicit DirectAttributeFilterPybind(std::shared_ptr<PythonValuedMorphologicalTree> valuedTree)
        : valuedTree_(requireReconstructionFilterOwner(std::move(valuedTree), "DirectAttributeFilter")),
          topologyMutationVersion_(valuedTree_->topology().getMutationVersion()) {}

    /**
     * @brief Applies direct reconstruction and returns a NumPy image.
     * @param nodePreservationMask Dense node-preservation decisions.
     * @return Newly allocated NumPy reconstruction.
     */
    [[nodiscard]] py::array applyDirectAttributeFilter(const NodePreservationMask& nodePreservationMask) const {
        valuedTree_->topology().requireMutationVersion(topologyMutationVersion_, "DirectAttributeFilter.apply_direct_attribute_filter");
        return valuedTree_->visit([&](const auto& native) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            return pybind_utils::toNumpy(DirectAttributeFilter<Altitude>(*native).applyDirectAttributeFilter(nodePreservationMask));
        });
    }
};

/** @brief Type-erased Python wrapper for hard subtractive attribute reconstruction. */
class SubtractiveAttributeFilterPybind {
  private:
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree_; ///< Python owner retained for the wrapper lifetime.
    std::size_t topologyMutationVersion_;                      ///< Topology version captured during construction.

  public:
    /**
     * @brief Creates a hard subtractive reconstruction wrapper.
     * @param valuedTree Python-owned valued tree used by subsequent reconstructions.
     */
    explicit SubtractiveAttributeFilterPybind(std::shared_ptr<PythonValuedMorphologicalTree> valuedTree)
        : valuedTree_(requireReconstructionFilterOwner(std::move(valuedTree), "SubtractiveAttributeFilter")),
          topologyMutationVersion_(valuedTree_->topology().getMutationVersion()) {}

    /**
     * @brief Applies hard subtractive reconstruction and returns a NumPy image.
     * @param nodePreservationMask Dense node-preservation decisions.
     * @return Newly allocated NumPy reconstruction in the signed altitude-difference type.
     */
    [[nodiscard]] py::array applySubtractiveAttributeFilter(const NodePreservationMask& nodePreservationMask) const {
        valuedTree_->topology().requireMutationVersion(topologyMutationVersion_,
                                                       "SubtractiveAttributeFilter.apply_subtractive_attribute_filter");
        return valuedTree_->visit([&](const auto& native) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            using OutputValue = AltitudeDifference<Altitude>;
            const auto valuedTreeView = native->asView();
            valuedTreeView.requireTopologyUnchanged("SubtractiveAttributeFilter.apply_subtractive_attribute_filter");
            const std::vector<OutputValue> nodeContributions =
                detail::attribute_filtering::computeHardNodeContributions(valuedTreeView, nodePreservationMask);
            std::vector<OutputValue> pixels = detail::tree_altitude::reconstructNodeContributionValues(
                native->topology(), std::span<const OutputValue>(nodeContributions),
                "SubtractiveAttributeFilter.apply_subtractive_attribute_filter");
            return pybind_utils::toNumpyOwned2D(std::move(pixels), native->topology().numRows(), native->topology().numColumns());
        });
    }
};

/** @brief Type-erased Python wrapper for soft subtractive attribute reconstruction. */
class SoftSubtractiveAttributeFilterPybind {
  private:
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree_; ///< Python owner retained for the wrapper lifetime.
    std::size_t topologyMutationVersion_;                      ///< Topology version captured during construction.

    /**
     * @brief Applies soft subtractive reconstruction for one floating-point score type.
     * @tparam Real Floating-point score and output type.
     * @param nodePreservationScores Dense one-dimensional NumPy score array.
     * @return Newly allocated NumPy reconstruction in `Real`.
     */
    template <std::floating_point Real> [[nodiscard]] py::array applyTyped(py::array nodePreservationScores) const {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(nodePreservationScores), valuedTree_->topology(), "node_preservation_scores");
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        const std::span<const Real> scores(values, static_cast<std::size_t>(valuedTree_->topology().numInternalNodeSlots()));
        return valuedTree_->visit([&](const auto& native) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            return pybind_utils::toNumpy(SoftSubtractiveAttributeFilter<Altitude, Real>(*native).applySoftSubtractiveAttributeFilter(scores));
        });
    }

  public:
    /**
     * @brief Creates a soft subtractive reconstruction wrapper.
     * @param valuedTree Python-owned valued tree used by subsequent reconstructions.
     */
    explicit SoftSubtractiveAttributeFilterPybind(std::shared_ptr<PythonValuedMorphologicalTree> valuedTree)
        : valuedTree_(requireReconstructionFilterOwner(std::move(valuedTree), "SoftSubtractiveAttributeFilter")),
          topologyMutationVersion_(valuedTree_->topology().getMutationVersion()) {}

    /**
     * @brief Applies soft subtractive reconstruction and returns a NumPy image.
     * @param nodePreservationScores Dense float32 or float64 node-preservation scores.
     * @return Newly allocated NumPy reconstruction with the input floating-point dtype.
     */
    [[nodiscard]] py::array applySoftSubtractiveAttributeFilter(py::array nodePreservationScores) const {
        valuedTree_->topology().requireMutationVersion(topologyMutationVersion_,
                                                       "SoftSubtractiveAttributeFilter.apply_soft_subtractive_attribute_filter");
        if (pybind_utils::parseFloatingArrayDType(nodePreservationScores, "node_preservation_scores") == pybind_utils::FloatingDType::Float64) {
            return applyTyped<double>(std::move(nodePreservationScores));
        }
        return applyTyped<float>(std::move(nodePreservationScores));
    }
};

} // namespace mmcfilters::pybindings
