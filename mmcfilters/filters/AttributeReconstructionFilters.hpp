#pragma once

#include "NodeDecisionMasks.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/ValuedMorphologicalTree.hpp"
#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../utils/CommittedImageAccess.hpp"
#include "../utils/Contract.hpp"

#include <cmath>
#include <concepts>
#include <span>
#include <stack>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace mmcfilters {

/** @brief Policy tag for parent-altitude propagation after node rejection. */
struct DirectReconstruction {};

/** @brief Policy tag for independent zero-baseline modulation of node residues. */
struct SubtractiveResidueModulation {};

namespace detail::attribute_filtering {

inline void requireNodePreservationMaskShape(const MorphologicalTree& tree, const NodePreservationMask& nodePreservationMask, const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(
        nodePreservationMask.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
        throw std::invalid_argument(std::string(context) + " nodePreservationMask size must match the internal node slot count."));
}

template <std::floating_point Real>
inline void requireNodePreservationScores(const MorphologicalTree& tree, std::span<const Real> nodePreservationScores, const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(
        nodePreservationScores.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
        throw std::invalid_argument(std::string(context) + " nodePreservationScores size must match the internal node slot count."));
    if constexpr (contract::validationsEnabled) {
        for (std::size_t index = 0; index < nodePreservationScores.size(); ++index) {
            const Real score = nodePreservationScores[index];
            if (!std::isfinite(score) || score < Real{0} || score > Real{1}) {
                throw std::invalid_argument(std::string(context) + " requires every nodePreservationScore to be finite and in [0, 1].");
            }
        }
    }
}

template <AltitudeValue T>
inline void applyDirectReconstruction(ValuedMorphologicalTreeView<T> valuedTree, const NodePreservationMask& nodePreservationMask, ImagePtr<T> output) {
    const char* context = "DirectAttributeFilter::applyDirectAttributeFilter";
    valuedTree.requireTopologyUnchanged(context);
    const MorphologicalTree& tree = valuedTree.topology();
    requireNodePreservationMaskShape(tree, nodePreservationMask, context);
    MMCFILTERS_CONTRACT_REQUIRE(nodePreservationMask[static_cast<std::size_t>(tree.root())],
                                throw std::invalid_argument(std::string(context) + " requires the root node to be preserved."));
    MMCFILTERS_CONTRACT_REQUIRE(output != nullptr, throw std::invalid_argument(std::string(context) + " requires a non-null output image."));
    MMCFILTERS_CONTRACT_REQUIRE(output->getNumRows() == tree.numRows() && output->getNumColumns() == tree.numColumns(),
                                throw std::invalid_argument(std::string(context) + " output image shape must match the tree image domain."));

    const std::span<const T> nodeAltitudes = valuedTree.nodeAltitudes();
    std::vector<T> reconstructedNodeAltitudes(static_cast<std::size_t>(tree.numInternalNodeSlots()), T{});
    const NodeId root = tree.root();
    reconstructedNodeAltitudes[static_cast<std::size_t>(root)] = nodeAltitudes[static_cast<std::size_t>(root)];

    std::stack<NodeId> pending;
    pending.push(root);
    while (!pending.empty()) {
        const NodeId nodeId = pending.top();
        pending.pop();
        for (NodeId childId : detail::CommittedTreeAccess::children(tree, nodeId)) {
            reconstructedNodeAltitudes[static_cast<std::size_t>(childId)] =
                nodePreservationMask[static_cast<std::size_t>(childId)] ? nodeAltitudes[static_cast<std::size_t>(childId)]
                                                                       : reconstructedNodeAltitudes[static_cast<std::size_t>(nodeId)];
            pending.push(childId);
        }
    }

    T* pixels = output->rawData();
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (PixelId pixel : detail::CommittedTreeAccess::properParts(tree, nodeId)) {
            pixels[pixel] = reconstructedNodeAltitudes[static_cast<std::size_t>(nodeId)];
        }
    }
}

template <AltitudeValue T>
[[nodiscard]] inline std::vector<AltitudeDifference<T>> computeHardNodeContributions(
    ValuedMorphologicalTreeView<T> valuedTree, const NodePreservationMask& nodePreservationMask) {
    const MorphologicalTree& tree = valuedTree.topology();
    requireNodePreservationMaskShape(tree, nodePreservationMask, "SubtractiveAttributeFilter::applySubtractiveAttributeFilter");
    std::vector<AltitudeDifference<T>> nodeContributions(static_cast<std::size_t>(tree.numInternalNodeSlots()), AltitudeDifference<T>{});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (nodePreservationMask[static_cast<std::size_t>(nodeId)]) {
            nodeContributions[static_cast<std::size_t>(nodeId)] = valuedTree.nodeResidue(nodeId);
        }
    }
    return nodeContributions;
}

template <AltitudeValue T, std::floating_point Real>
[[nodiscard]] inline std::vector<Real> computeSoftNodeContributions(ValuedMorphologicalTreeView<T> valuedTree,
                                                                    std::span<const Real> nodePreservationScores) {
    const MorphologicalTree& tree = valuedTree.topology();
    requireNodePreservationScores(tree, nodePreservationScores, "SoftSubtractiveAttributeFilter::applySoftSubtractiveAttributeFilter");
    std::vector<Real> nodeContributions(static_cast<std::size_t>(tree.numInternalNodeSlots()), Real{});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        nodeContributions[static_cast<std::size_t>(nodeId)] =
            static_cast<Real>(valuedTree.nodeResidue(nodeId)) * nodePreservationScores[static_cast<std::size_t>(nodeId)];
    }
    return nodeContributions;
}

} // namespace detail::attribute_filtering

/**
 * @brief Direct attribute filter using reconstructed-parent altitude propagation.
 * @tparam T Node-altitude value type.
 */
template <AltitudeValue T> class DirectAttributeFilter {
  private:
    ValuedMorphologicalTreeView<T> valuedTree_; ///< Validated non-owning valued-tree view.

  public:
    /** @brief Reconstruction policy implemented by this filter. */
    using ReconstructionPolicy = DirectReconstruction;

    /**
     * @brief Creates a direct filter over a valued-tree view.
     * @param valuedTree Valued tree whose topology and altitudes drive reconstruction.
     */
    explicit DirectAttributeFilter(ValuedMorphologicalTreeView<T> valuedTree) : valuedTree_(valuedTree) {}

    /**
     * @brief Creates a direct filter over an owning valued tree.
     * @param valuedTree Valued tree whose topology and altitudes drive reconstruction.
     */
    explicit DirectAttributeFilter(const ValuedMorphologicalTree<T>& valuedTree) : valuedTree_(valuedTree.asView()) {}

    /**
     * @brief Applies direct reconstruction under node-preservation decisions.
     * @param nodePreservationMask Dense mask where `true` preserves the node altitude.
     * @return Reconstructed image in the original altitude type.
     */
    [[nodiscard]] ImagePtr<T> applyDirectAttributeFilter(const NodePreservationMask& nodePreservationMask) const {
        ImagePtr<T> output = detail::CommittedImageAccess::create<T>(valuedTree_.topology().numRows(), valuedTree_.topology().numColumns());
        detail::attribute_filtering::applyDirectReconstruction(valuedTree_, nodePreservationMask, output);
        return output;
    }
};

/**
 * @brief Subtractive filter using independently gated zero-baseline residues.
 * @tparam T Node-altitude value type.
 */
template <AltitudeValue T> class SubtractiveAttributeFilter {
  private:
    ValuedMorphologicalTreeView<T> valuedTree_; ///< Validated non-owning valued-tree view.

  public:
    /** @brief Reconstruction policy implemented by this filter. */
    using ReconstructionPolicy = SubtractiveResidueModulation;

    /** @brief Signed value type produced by residue accumulation. */
    using OutputValue = AltitudeDifference<T>;

    /**
     * @brief Creates a hard subtractive filter over a valued-tree view.
     * @param valuedTree Valued tree whose residues drive reconstruction.
     */
    explicit SubtractiveAttributeFilter(ValuedMorphologicalTreeView<T> valuedTree) : valuedTree_(valuedTree) {}

    /**
     * @brief Creates a hard subtractive filter over an owning valued tree.
     * @param valuedTree Valued tree whose residues drive reconstruction.
     */
    explicit SubtractiveAttributeFilter(const ValuedMorphologicalTree<T>& valuedTree) : valuedTree_(valuedTree.asView()) {}

    /**
     * @brief Reconstructs from residues selected by preservation decisions.
     * @param nodePreservationMask Dense mask where `true` retains the node residue.
     * @return Zero-baseline subtractive reconstruction in the signed difference type.
     */
    [[nodiscard]] ImagePtr<OutputValue> applySubtractiveAttributeFilter(const NodePreservationMask& nodePreservationMask) const {
        valuedTree_.requireTopologyUnchanged("SubtractiveAttributeFilter::applySubtractiveAttributeFilter");
        const std::vector<OutputValue> nodeContributions = detail::attribute_filtering::computeHardNodeContributions(valuedTree_, nodePreservationMask);
        return TreeAltitudeAlgorithms::reconstructFromNodeContributions(
            valuedTree_.topology(), std::span<const OutputValue>(nodeContributions),
            "SubtractiveAttributeFilter::applySubtractiveAttributeFilter");
    }
};

/**
 * @brief Soft subtractive filter using scores in [0,1] to gate zero-baseline residues.
 * @tparam T Node-altitude value type.
 * @tparam Real Floating-point score and output type.
 */
template <AltitudeValue T, std::floating_point Real = float> class SoftSubtractiveAttributeFilter {
  private:
    ValuedMorphologicalTreeView<T> valuedTree_; ///< Validated non-owning valued-tree view.

  public:
    /** @brief Reconstruction policy implemented by this filter. */
    using ReconstructionPolicy = SubtractiveResidueModulation;

    /**
     * @brief Creates a soft subtractive filter over a valued-tree view.
     * @param valuedTree Valued tree whose residues drive reconstruction.
     */
    explicit SoftSubtractiveAttributeFilter(ValuedMorphologicalTreeView<T> valuedTree) : valuedTree_(valuedTree) {}

    /**
     * @brief Creates a soft subtractive filter over an owning valued tree.
     * @param valuedTree Valued tree whose residues drive reconstruction.
     */
    explicit SoftSubtractiveAttributeFilter(const ValuedMorphologicalTree<T>& valuedTree) : valuedTree_(valuedTree.asView()) {}

    /**
     * @brief Reconstructs from residues modulated by preservation scores.
     * @param nodePreservationScores Dense score buffer indexed by the internal node-slot domain.
     * @return Zero-baseline subtractive reconstruction in the score type.
     */
    [[nodiscard]] ImagePtr<Real> applySoftSubtractiveAttributeFilter(std::span<const Real> nodePreservationScores) const {
        valuedTree_.requireTopologyUnchanged("SoftSubtractiveAttributeFilter::applySoftSubtractiveAttributeFilter");
        const std::vector<Real> nodeContributions = detail::attribute_filtering::computeSoftNodeContributions(valuedTree_, nodePreservationScores);
        return TreeAltitudeAlgorithms::reconstructFromNodeContributions(
            valuedTree_.topology(), std::span<const Real>(nodeContributions),
            "SoftSubtractiveAttributeFilter::applySoftSubtractiveAttributeFilter");
    }
};

template <AltitudeValue T>
[[nodiscard]] inline ImagePtr<T> applyDirectAttributeFilter(ValuedMorphologicalTreeView<T> valuedTree,
                                                            const NodePreservationMask& nodePreservationMask) {
    return DirectAttributeFilter<T>(valuedTree).applyDirectAttributeFilter(nodePreservationMask);
}

template <AltitudeValue T>
[[nodiscard]] inline ImagePtr<T> applyDirectAttributeFilter(const ValuedMorphologicalTree<T>& valuedTree,
                                                            const NodePreservationMask& nodePreservationMask) {
    return DirectAttributeFilter<T>(valuedTree).applyDirectAttributeFilter(nodePreservationMask);
}

template <AltitudeValue T>
[[nodiscard]] inline ImagePtr<AltitudeDifference<T>> applySubtractiveAttributeFilter(
    ValuedMorphologicalTreeView<T> valuedTree, const NodePreservationMask& nodePreservationMask) {
    return SubtractiveAttributeFilter<T>(valuedTree).applySubtractiveAttributeFilter(nodePreservationMask);
}

template <AltitudeValue T>
[[nodiscard]] inline ImagePtr<AltitudeDifference<T>> applySubtractiveAttributeFilter(
    const ValuedMorphologicalTree<T>& valuedTree, const NodePreservationMask& nodePreservationMask) {
    return SubtractiveAttributeFilter<T>(valuedTree).applySubtractiveAttributeFilter(nodePreservationMask);
}

template <AltitudeValue T, std::floating_point Real>
[[nodiscard]] inline ImagePtr<Real> applySoftSubtractiveAttributeFilter(ValuedMorphologicalTreeView<T> valuedTree,
                                                                        std::span<const Real> nodePreservationScores) {
    return SoftSubtractiveAttributeFilter<T, Real>(valuedTree).applySoftSubtractiveAttributeFilter(nodePreservationScores);
}

template <AltitudeValue T, std::floating_point Real>
[[nodiscard]] inline ImagePtr<Real> applySoftSubtractiveAttributeFilter(const ValuedMorphologicalTree<T>& valuedTree,
                                                                        std::span<const Real> nodePreservationScores) {
    return SoftSubtractiveAttributeFilter<T, Real>(valuedTree).applySoftSubtractiveAttributeFilter(nodePreservationScores);
}

} // namespace mmcfilters
