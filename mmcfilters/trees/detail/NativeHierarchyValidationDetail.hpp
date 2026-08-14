#pragma once

#include "../NativeHierarchy.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

class NativeHierarchyValidation;
class TopologicalNativeHierarchyRecorder;

/**
 * @brief Move-only evidence that native topology buffers satisfy the generic
 * rooted partial-partition contract.
 *
 * The proof is bound to the node/pixel cardinalities and root id. It can
 * only be issued by the complete validator or by the recorder used by a
 * producer that establishes the same facts in its existing construction
 * passes.
 */
class NativeTopologyProof {
    friend class NativeHierarchyValidation;
    friend class TopologicalNativeHierarchyRecorder;

    /** @brief Number of nodes. */
    std::size_t numNodes_ = 0;
    /** @brief Number of proper parts. */
    std::size_t numPixels_ = 0;
    /** @brief Dense node identifier of the root. */
    NodeId root_ = InvalidNode;
    /** @brief Indicates whether the proof token contains validated evidence. */
    bool valid_ = false;

    /**
     * @brief Constructs `NativeTopologyProof` from the supplied inputs.
     *
     * @param numNodes Number of internal nodes.
     * @param numPixels Number of pixels.
     * @param root Root node of the operation.
     */
    NativeTopologyProof(std::size_t numNodes, std::size_t numPixels, NodeId root) noexcept
        : numNodes_(numNodes), numPixels_(numPixels), root_(root), valid_(true) {}

  public:
    /**
     * @brief Constructs a default `NativeTopologyProof`.
     */
    NativeTopologyProof() noexcept = default;
    /**
     * @brief Disables copy construction.
     */
    NativeTopologyProof(const NativeTopologyProof&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    NativeTopologyProof& operator=(const NativeTopologyProof&) = delete;
    /**
     * @brief Constructs a default `NativeTopologyProof`.
     */
    NativeTopologyProof(NativeTopologyProof&&) noexcept = default;
    /**
     * @brief Assigns the supplied object state.
     *
     * @return Reference to the assigned object.
     */
    NativeTopologyProof& operator=(NativeTopologyProof&&) noexcept = default;

    /**
     * @brief Validates matches.
     *
     * @param numNodes Number of internal nodes.
     * @param numPixels Number of pixels.
     * @param root Root node of the operation.
     */
    void requireMatches(std::size_t numNodes, std::size_t numPixels, NodeId root) const {
        if (!valid_ || numNodes_ != numNodes || numPixels_ != numPixels || root_ != root) {
            throw std::logic_error("Native topology proof does not match the materialized buffers.");
        }
    }
};

/**
 * @brief Move-only evidence that every emitted node has non-empty subtree
 * support.
 *
 * This narrower proof lets an assembler establish support inductively while
 * it emits nodes. A later validator still checks the complete parent and smallest-node
 * domains, but does not need another support-propagation pass.
 */
class NativeSubtreeSupportProof {
    friend class NativeSubtreeSupportRecorder;
    friend class NativeHierarchyValidation;

    /** @brief Number of supported nodes. */
    std::size_t numSupportedNodes_ = 0;
    /** @brief Indicates whether the proof token contains validated evidence. */
    bool valid_ = false;

    /**
     * @brief Constructs `NativeSubtreeSupportProof` from the supplied inputs.
     *
     * @param numSupportedNodes Count.
     */
    explicit NativeSubtreeSupportProof(std::size_t numSupportedNodes) noexcept : numSupportedNodes_(numSupportedNodes), valid_(true) {}

  public:
    /**
     * @brief Constructs a default `NativeSubtreeSupportProof`.
     */
    NativeSubtreeSupportProof() noexcept = default;
    /**
     * @brief Disables copy construction.
     */
    NativeSubtreeSupportProof(const NativeSubtreeSupportProof&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    NativeSubtreeSupportProof& operator=(const NativeSubtreeSupportProof&) = delete;
    /**
     * @brief Constructs a default `NativeSubtreeSupportProof`.
     */
    NativeSubtreeSupportProof(NativeSubtreeSupportProof&&) noexcept = default;
    /**
     * @brief Assigns the supplied object state.
     *
     * @return Reference to the assigned object.
     */
    NativeSubtreeSupportProof& operator=(NativeSubtreeSupportProof&&) noexcept = default;
};

/**
 * @brief Constant-space recorder for inductive subtree-support evidence.
 *
 * A producer calls `recordSupportedNode()` exactly when it creates a node from
 * at least one direct proper part or from an already supported child.
 */
class NativeSubtreeSupportRecorder {
    /** @brief Number of supported nodes. */
    std::size_t numSupportedNodes_ = 0;

  public:
    /**
     * @brief Records supported node.
     */
    void recordSupportedNode() {
        if (numSupportedNodes_ == std::numeric_limits<std::size_t>::max()) {
            throw std::length_error("Native subtree-support evidence exceeds the supported size.");
        }
        ++numSupportedNodes_;
    }

    /**
     * @brief Finalizes and returns the accumulated proof.
     *
     * @param expectedNodes Expected number of internal nodes.
     * @return Finalized and returns the accumulated proof.
     */
    [[nodiscard]] NativeSubtreeSupportProof finish(std::size_t expectedNodes) && {
        if (numSupportedNodes_ != expectedNodes) {
            throw std::logic_error("Native subtree-support evidence does not cover every node.");
        }
        return NativeSubtreeSupportProof(numSupportedNodes_);
    }
};

/**
 * @brief Constant-space proof recorder for producers that emit parents before
 * children in dense node-id order.
 *
 * `recordSupportedNode()` establishes connectivity by requiring every
 * non-root parent to have already been recorded. Its name makes the additional
 * producer obligation explicit: the recorded node must already be known to
 * have non-empty full subtree support.
 */
class TopologicalNativeHierarchyRecorder {
    /** @brief Expected nodes. */
    std::size_t expectedNodes_;
    /** @brief Expected proper parts. */
    std::size_t expectedProperParts_;
    /** @brief Dense node identifier of the expected root. */
    NodeId expectedRoot_;
    /** @brief Next node. */
    std::size_t nextNode_ = 0;
    /** @brief Next proper part. */
    std::size_t nextProperPart_ = 0;

  public:
    /**
     * @brief Constructs `TopologicalNativeHierarchyRecorder` from the supplied inputs.
     *
     * @param expectedNodes Expected number of internal nodes.
     * @param expectedProperParts Expected number of proper parts.
     * @param expectedRoot Expected root identifier.
     */
    TopologicalNativeHierarchyRecorder(std::size_t expectedNodes, std::size_t expectedProperParts, NodeId expectedRoot)
        : expectedNodes_(expectedNodes), expectedProperParts_(expectedProperParts), expectedRoot_(expectedRoot) {
        if (expectedNodes_ == 0 || expectedProperParts_ == 0) {
            throw std::invalid_argument("Native hierarchy evidence requires non-empty node and pixel domains.");
        }
        if (expectedNodes_ > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) ||
            expectedProperParts_ > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::length_error("Native hierarchy evidence exceeds the NodeId domain.");
        }
        if (expectedRoot_ != 0) {
            throw std::invalid_argument("Topological native hierarchy evidence requires the root to be emitted first.");
        }
    }

    /**
     * @brief Records supported node.
     *
     * @param node Node identifier.
     * @param parent Parent accumulator or node.
     */
    void recordSupportedNode(NodeId node, NodeId parent) {
        if (node != static_cast<NodeId>(nextNode_)) {
            throw std::logic_error("Topological native hierarchy nodes must be recorded in dense id order.");
        }
        if (node == expectedRoot_) {
            if (parent != node) {
                throw std::logic_error("The recorded native hierarchy root must be self-parented.");
            }
        } else if (parent < 0 || parent >= node) {
            throw std::logic_error("A recorded native hierarchy parent must already be connected.");
        }
        ++nextNode_;
    }

    /**
     * @brief Records proper part.
     *
     * @param pixel Proper-part identifier.
     * @param smallestNode Smallest node associated with the pixel.
     */
    void recordProperPart(PixelId pixel, NodeId smallestNode) {
        if (pixel != static_cast<PixelId>(nextProperPart_)) {
            throw std::logic_error("Native proper parts must be recorded in dense id order.");
        }
        if (smallestNode < 0 || smallestNode >= static_cast<NodeId>(expectedNodes_)) {
            throw std::logic_error("A recorded smallest node is outside the node domain.");
        }
        ++nextProperPart_;
    }

    /**
     * @brief Records one smallest node in a producer-established permutation pass.
     *
     * The producer is responsible for visiting every proper part exactly once.
     * This variant avoids a second dense-domain scan when the construction
     * algorithm already traverses a proven permutation of the proper parts.
     *
     * @param smallestNode Smallest node associated with the pixel.
     */
    void recordSmallestNode(NodeId smallestNode) {
        if (nextProperPart_ >= expectedProperParts_) {
            throw std::logic_error("Native proper-part evidence exceeds the declared domain.");
        }
        if (smallestNode < 0 || smallestNode >= static_cast<NodeId>(expectedNodes_)) {
            throw std::logic_error("A recorded smallest node is outside the node domain.");
        }
        ++nextProperPart_;
    }

    /**
     * @brief Finalizes and returns the accumulated proof.
     *
     * @return Finalized and returns the accumulated proof.
     */
    [[nodiscard]] NativeTopologyProof finish() && {
        if (nextNode_ != expectedNodes_ || nextProperPart_ != expectedProperParts_) {
            throw std::logic_error("Native hierarchy evidence does not cover every node and proper part.");
        }
        return NativeTopologyProof(expectedNodes_, expectedProperParts_, expectedRoot_);
    }
};

/**
 * @brief Owns native hierarchy buffers after structural validation.
 *
 * @tparam T Altitude type stored by the validated hierarchy.
 */
template <AltitudeValue T> struct ValidatedNativeHierarchyStorage {
    /** @brief Dense node identifier of the node parent. */
    std::vector<NodeId> parent;
    /** @brief Dense identifier of the smallest node. */
    std::vector<NodeId> smallestNodeMap;
    /** @brief Stores node altitudes. */
    std::vector<T> nodeAltitudes;
    /** @brief Dense node identifier of the root. */
    NodeId root = InvalidNode;
    /** @brief Grid domain2 d. */
    std::optional<GridDomain2D> gridDomain2D;
    /** @brief Semantics. */
    MorphologicalTreeSemantics semantics;
    /** @brief Topology proof. */
    NativeTopologyProof topologyProof;
};

/**
 * @brief Owning, move-only native hierarchy paired with structural evidence.
 *
 * Internal producers move their existing buffers into this representation.
 * Public span imports are validated first and copied once into it. In both
 * cases materialization subsequently transfers the buffers into the valuedTree
 * owner without another copy.
 */
template <AltitudeValue T> class ValidatedNativeHierarchy {
    /** @brief Storage. */
    ValidatedNativeHierarchyStorage<T> storage_;

  public:
    /**
     * @brief Takes ownership of structurally validated native hierarchy storage.
     *
     * @param storage Storage object that owns the native hierarchy arrays.
     */
    explicit ValidatedNativeHierarchy(ValidatedNativeHierarchyStorage<T>&& storage) : storage_(std::move(storage)) {
        storage_.topologyProof.requireMatches(storage_.parent.size(), storage_.smallestNodeMap.size(), storage_.root);
    }

    /**
     * @brief Disables copy construction.
     */
    ValidatedNativeHierarchy(const ValidatedNativeHierarchy&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    ValidatedNativeHierarchy& operator=(const ValidatedNativeHierarchy&) = delete;
    /**
     * @brief Transfers ownership of the validated native hierarchy.
     */
    ValidatedNativeHierarchy(ValidatedNativeHierarchy&&) noexcept = default;
    /**
     * @brief Assigns the supplied object state.
     *
     * @return Reference to the assigned object.
     */
    ValidatedNativeHierarchy& operator=(ValidatedNativeHierarchy&&) noexcept = default;

    /**
     * @brief Releases ownership of the stored native hierarchy.
     *
     * @return Owned native hierarchy storage.
     */
    [[nodiscard]] ValidatedNativeHierarchyStorage<T> release() && { return std::move(storage_); }
};

/**
 * @brief Central validator and materialization boundary for native hierarchy
 * buffers.
 */
class NativeHierarchyValidation {
    /**
     * @brief Validates descriptor.
     *
     * @param hierarchy Native hierarchy data to validate or convert.
     */
    template <AltitudeValue T> static void validateDescriptor(const NativeHierarchyView<T>& hierarchy) {
        if (hierarchy.parent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) ||
            hierarchy.smallestNodeMap.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native hierarchy domain exceeds the NodeId range.");
        }
        if (hierarchy.parent.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one internal node.");
        }
        if (hierarchy.smallestNodeMap.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one proper part.");
        }
        if (hierarchy.nodeAltitudes.size() != hierarchy.parent.size()) {
            throw std::invalid_argument("Native hierarchy altitude size must match the internal-node domain.");
        }

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        if (hierarchy.root < 0 || hierarchy.root >= numNodes) {
            throw std::invalid_argument("Native hierarchy requires a valid root node id.");
        }
        if (hierarchy.gridDomain2D && hierarchy.smallestNodeMap.size() != hierarchy.gridDomain2D->size("Native hierarchy 2D pixel domain")) {
            throw std::invalid_argument("Native hierarchy pixel domain must match the attached 2D grid.");
        }

        switch (hierarchy.semantics.kind) {
        case MorphologicalTreeKind::Generic:
        case MorphologicalTreeKind::MaxTree:
        case MorphologicalTreeKind::MinTree:
        case MorphologicalTreeKind::TreeOfShapes:
        case MorphologicalTreeKind::UnrestrictedResidualTree:
        case MorphologicalTreeKind::SaturatedResidualTree:
            break;
        default:
            throw std::invalid_argument("Native hierarchy kind is not supported.");
        }

        validateMorphologicalTreeSemantics(hierarchy.semantics);
        const auto requireMatchingGrid = [&hierarchy](const RegularGridAdjacency2D& adjacency, const char* context) {
            if (!hierarchy.gridDomain2D) {
                throw std::invalid_argument(std::string(context) + " requires an attached 2D grid domain.");
            }
            if (adjacency.getNumRows() != hierarchy.gridDomain2D->rows || adjacency.getNumColumns() != hierarchy.gridDomain2D->columns) {
                throw std::invalid_argument(std::string(context) + " must match the attached 2D grid.");
            }
        };
        if (const auto* context = std::get_if<SharedAdjacencyContext>(&hierarchy.semantics.constructionContext)) {
            requireMatchingGrid(context->adjacency, "SharedAdjacencyContext adjacency");
        } else if (const auto* context = std::get_if<SaturatedResidualContext>(&hierarchy.semantics.constructionContext)) {
            requireMatchingGrid(context->adjacency, "SaturatedResidualContext adjacency");
            if (context->infinityPixel < 0 ||
                static_cast<std::size_t>(context->infinityPixel) >= hierarchy.gridDomain2D->size("Saturated residual domain")) {
                throw std::invalid_argument("SaturatedResidualContext infinity pixel must belong to the attached 2D grid.");
            }
        } else if (const auto* convention = std::get_if<TopographicConvention>(&hierarchy.semantics.constructionContext)) {
            if (convention->infinityPixel < 0) {
                throw std::invalid_argument("TopographicConvention infinity pixel must be non-negative.");
            }
            if (!hierarchy.gridDomain2D) {
                throw std::invalid_argument("TopographicConvention requires an attached 2D grid domain.");
            }
            const std::int64_t extension = convention->domainExtension == TopographicDomainExtension::ExteriorRing ? 1 : -1;
            const std::int64_t activeRows = 2 * static_cast<std::int64_t>(hierarchy.gridDomain2D->rows) + extension;
            const std::int64_t activeColumns = 2 * static_cast<std::int64_t>(hierarchy.gridDomain2D->columns) + extension;
            if (activeRows <= 0 || activeColumns <= 0 || static_cast<std::int64_t>(convention->infinityPixel) >= activeRows * activeColumns) {
                throw std::invalid_argument("TopographicConvention infinity pixel must belong to the active topographic domain.");
            }
            if (const auto* immersion = std::get_if<ComplementaryGridImmersion>(&convention->immersion)) {
                requireMatchingGrid(immersion->complementaryAdjacencies.minAdjacency, "Topographic minimum adjacency");
                requireMatchingGrid(immersion->complementaryAdjacencies.maxAdjacency, "Topographic maximum adjacency");
            }
        }
    }

  public:
    /**
     * @brief Validates descriptor only.
     *
     * @param hierarchy Native hierarchy data to validate or convert.
     */
    template <AltitudeValue T> static void validateDescriptorOnly(const NativeHierarchyView<T>& hierarchy) { validateDescriptor(hierarchy); }

    /**
     * @brief Performs the complete structural validation required for a public
     * raw-buffer import.
     *
     * @param hierarchy Hierarchy data.
     * @return Result of performs the complete structural validation required for a public raw-buffer import.
     */
    template <AltitudeValue T> [[nodiscard]] static NativeTopologyProof validateComplete(const NativeHierarchyView<T>& hierarchy) {
        validateDescriptor(hierarchy);

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        std::vector<int> remainingChildren(hierarchy.parent.size(), 0);

        for (NodeId node = 0; node < numNodes; ++node) {
            const NodeId parent = hierarchy.parent[static_cast<std::size_t>(node)];
            if (node == hierarchy.root) {
                if (parent != node) {
                    throw std::invalid_argument("Native hierarchy root must be self-parented.");
                }
                continue;
            }
            if (parent < 0 || parent >= numNodes) {
                throw std::invalid_argument("Native hierarchy parent is outside the internal-node domain.");
            }
            if (parent == node) {
                throw std::invalid_argument("Native hierarchy contains a detached self-parented non-root node.");
            }
            ++remainingChildren[static_cast<std::size_t>(parent)];
        }

        std::vector<std::uint8_t> subtreeHasSupport(hierarchy.parent.size(), std::uint8_t{0});
        for (PixelId pixel = 0; pixel < static_cast<PixelId>(hierarchy.smallestNodeMap.size()); ++pixel) {
            const NodeId smallestNode = hierarchy.smallestNodeMap[static_cast<std::size_t>(pixel)];
            if (smallestNode < 0 || smallestNode >= numNodes) {
                throw std::invalid_argument("Native hierarchy smallest node is outside the internal-node domain.");
            }
            subtreeHasSupport[static_cast<std::size_t>(smallestNode)] = 1;
        }

        std::vector<NodeId> ready;
        ready.reserve(hierarchy.parent.size());
        for (NodeId node = 0; node < numNodes; ++node) {
            if (remainingChildren[static_cast<std::size_t>(node)] == 0) {
                ready.push_back(node);
            }
        }

        std::size_t processedNodes = 0;
        for (std::size_t index = 0; index < ready.size(); ++index) {
            const NodeId node = ready[index];
            if (subtreeHasSupport[static_cast<std::size_t>(node)] == 0) {
                throw std::invalid_argument("Native hierarchy contains a node whose subtree support is empty.");
            }
            ++processedNodes;
            if (node == hierarchy.root) {
                continue;
            }

            const NodeId parent = hierarchy.parent[static_cast<std::size_t>(node)];
            subtreeHasSupport[static_cast<std::size_t>(parent)] = 1;
            int& parentRemaining = remainingChildren[static_cast<std::size_t>(parent)];
            --parentRemaining;
            if (parentRemaining == 0) {
                ready.push_back(parent);
            }
        }
        if (processedNodes != hierarchy.parent.size()) {
            throw std::invalid_argument("Native hierarchy contains a cycle or nodes disconnected from the root.");
        }

        return NativeTopologyProof(hierarchy.parent.size(), hierarchy.smallestNodeMap.size(), hierarchy.root);
    }

    /**
     * @brief Validates parent connectivity and the smallest-node map while consuming
     * support evidence already established by the producer.
     *
     * `validateSmallestNode` runs inside the existing proper-part pass, allowing a
     * native producer to retain its reconstruction oracle without another scan.
     *
     * @param hierarchy Hierarchy data.
     * @param supportProof Proof that every node support is non-empty.
     * @param validateSmallestNode Callback that revalidates each smallest-node mapping.
     * @return Validation result for parent connectivity and the smallest-node map while consuming support evidence already established by the producer.
     */
    template <AltitudeValue T, class SmallestNodeValidator>
    [[nodiscard]] static NativeTopologyProof validateWithEstablishedSupport(const NativeHierarchyView<T>& hierarchy, NativeSubtreeSupportProof&& supportProof,
                                                                            SmallestNodeValidator&& validateSmallestNode) {
        validateDescriptor(hierarchy);
        if (!supportProof.valid_ || supportProof.numSupportedNodes_ != hierarchy.parent.size()) {
            throw std::logic_error("Native subtree-support proof does not match the node domain.");
        }

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        if (hierarchy.parent[static_cast<std::size_t>(hierarchy.root)] != hierarchy.root) {
            throw std::invalid_argument("Native hierarchy root must be self-parented.");
        }

        std::vector<std::uint8_t> parentState(hierarchy.parent.size(), std::uint8_t{0});
        parentState[static_cast<std::size_t>(hierarchy.root)] = 2;
        for (NodeId node = 0; node < numNodes; ++node) {
            if (node == hierarchy.root) {
                continue;
            }

            NodeId cursor = node;
            while (parentState[static_cast<std::size_t>(cursor)] == 0) {
                parentState[static_cast<std::size_t>(cursor)] = 1;
                const NodeId parent = hierarchy.parent[static_cast<std::size_t>(cursor)];
                if (parent < 0 || parent >= numNodes || parent == cursor) {
                    throw std::invalid_argument("Native hierarchy parent is self-referential or outside the internal-node domain.");
                }
                cursor = parent;
            }
            if (parentState[static_cast<std::size_t>(cursor)] == 1) {
                throw std::invalid_argument("Native hierarchy contains a non-root cycle.");
            }

            cursor = node;
            while (parentState[static_cast<std::size_t>(cursor)] == 1) {
                parentState[static_cast<std::size_t>(cursor)] = 2;
                cursor = hierarchy.parent[static_cast<std::size_t>(cursor)];
            }
        }

        for (PixelId pixel = 0; pixel < static_cast<PixelId>(hierarchy.smallestNodeMap.size()); ++pixel) {
            const NodeId smallestNode = hierarchy.smallestNodeMap[static_cast<std::size_t>(pixel)];
            if (smallestNode < 0 || smallestNode >= numNodes) {
                throw std::invalid_argument("Native hierarchy smallest node is outside the internal-node domain.");
            }
            validateSmallestNode(pixel, smallestNode);
        }

        return NativeTopologyProof(hierarchy.parent.size(), hierarchy.smallestNodeMap.size(), hierarchy.root);
    }
};

template<AltitudeValue T>
[[nodiscard]] ValidatedNativeHierarchy<T>
/**
 * @brief Validates a native hierarchy and copies it into owned storage.
 *
 * @param hierarchy Native hierarchy data to validate or convert.
 * @return Owning copy of the validated native hierarchy.
 */
validateAndCopyNativeHierarchy(
    NativeHierarchyView<T> hierarchy) {
    NativeTopologyProof proof = NativeHierarchyValidation::validateComplete(hierarchy);
    return ValidatedNativeHierarchy<T>(
        ValidatedNativeHierarchyStorage<T>{std::vector<NodeId>(hierarchy.parent.begin(), hierarchy.parent.end()),
                                           std::vector<NodeId>(hierarchy.smallestNodeMap.begin(), hierarchy.smallestNodeMap.end()),
                                           std::vector<T>(hierarchy.nodeAltitudes.begin(), hierarchy.nodeAltitudes.end()), hierarchy.root, hierarchy.gridDomain2D,
                                           std::move(hierarchy.semantics), std::move(proof)});
}

template<AltitudeValue T>
[[nodiscard]] ValidatedNativeHierarchy<T>
/**
 * @brief Creates validated native hierarchy.
 *
 * @param parent Parent node indexed by internal node identifier.
 * @param smallestNodeMap Smallest node indexed by pixel identifier.
 * @param nodeAltitudes Altitude values indexed by internal node identifier.
 * @param root Root node of the operation.
 * @param gridDomain2D Two-dimensional grid metadata associated with the hierarchy.
 * @param semantics Semantic metadata associated with the hierarchy.
 * @param proof Validation proof to attach to the imported hierarchy.
 * @return Created validated native hierarchy.
 */
makeValidatedNativeHierarchy(
    std::vector<NodeId>&& parent,
    std::vector<NodeId>&& smallestNodeMap,
    std::vector<T>&& nodeAltitudes,
    NodeId root,
    std::optional<GridDomain2D> gridDomain2D,
    MorphologicalTreeSemantics semantics,
    NativeTopologyProof&& proof) {
    const NativeHierarchyView<T> hierarchy{parent, smallestNodeMap, nodeAltitudes, root, gridDomain2D, semantics};
    NativeHierarchyValidation::validateDescriptorOnly(hierarchy);
    proof.requireMatches(parent.size(), smallestNodeMap.size(), root);

    return ValidatedNativeHierarchy<T>(ValidatedNativeHierarchyStorage<T>{std::move(parent), std::move(smallestNodeMap), std::move(nodeAltitudes), root,
                                                                          gridDomain2D, std::move(semantics), std::move(proof)});
}

} // namespace mmcfilters::detail
