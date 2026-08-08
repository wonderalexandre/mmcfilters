#pragma once

#include "../NativeHierarchy.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

class NativeHierarchyValidation;
class TopologicalNativeHierarchyRecorder;

/**
 * @brief Move-only evidence that native topology buffers satisfy the generic
 * rooted partial-partition contract.
 *
 * The proof is bound to the node/proper-part cardinalities and root id. It can
 * only be issued by the complete validator or by the recorder used by a
 * producer that establishes the same facts in its existing construction
 * passes.
 */
class NativeTopologyProof {
    friend class NativeHierarchyValidation;
    friend class TopologicalNativeHierarchyRecorder;

    /** @brief Stores the num nodes. */
    std::size_t numNodes_ = 0;
    /** @brief Stores the num proper parts. */
    std::size_t numProperParts_ = 0;
    /** @brief Stores the root. */
    NodeId root_ = InvalidNode;
    /** @brief Indicates whether the proof token contains validated evidence. */
    bool valid_ = false;

    /**
     * @brief Constructs `NativeTopologyProof` from the supplied inputs.
     *
     * @param numNodes Number of internal nodes.
     * @param numProperParts Number of proper parts.
     * @param root Root node of the operation.
     */
    NativeTopologyProof(std::size_t numNodes, std::size_t numProperParts, NodeId root) noexcept
        : numNodes_(numNodes), numProperParts_(numProperParts), root_(root), valid_(true) {}

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
     * @param numProperParts Number of proper parts.
     * @param root Root node of the operation.
     */
    void requireMatches(std::size_t numNodes, std::size_t numProperParts, NodeId root) const {
        if (!valid_ || numNodes_ != numNodes || numProperParts_ != numProperParts || root_ != root) {
            throw std::logic_error("Native topology proof does not match the materialized buffers.");
        }
    }
};

/**
 * @brief Move-only evidence that every emitted node has non-empty subtree
 * support.
 *
 * This narrower proof lets an assembler establish support inductively while
 * it emits nodes. A later validator still checks the complete parent and owner
 * domains, but does not need another support-propagation pass.
 */
class NativeSubtreeSupportProof {
    friend class NativeSubtreeSupportRecorder;
    friend class NativeHierarchyValidation;

    /** @brief Stores the num supported nodes. */
    std::size_t numSupportedNodes_ = 0;
    /** @brief Indicates whether the proof token contains validated evidence. */
    bool valid_ = false;

    /**
     * @brief Constructs `NativeSubtreeSupportProof` from the supplied inputs.
     *
     * @param numSupportedNodes Count represented by `numSupportedNodes`.
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
    /** @brief Stores the num supported nodes. */
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
    /** @brief Stores the expected nodes. */
    std::size_t expectedNodes_;
    /** @brief Stores the expected proper parts. */
    std::size_t expectedProperParts_;
    /** @brief Stores the expected root. */
    NodeId expectedRoot_;
    /** @brief Stores the next node. */
    std::size_t nextNode_ = 0;
    /** @brief Stores the next proper part. */
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
            throw std::invalid_argument("Native hierarchy evidence requires non-empty node and proper-part domains.");
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
     * @param node Node identifier used by the operation.
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
     * @param properPart Proper-part identifier.
     * @param owner Owner associated with the view or iterator.
     */
    void recordProperPart(NodeId properPart, NodeId owner) {
        if (properPart != static_cast<NodeId>(nextProperPart_)) {
            throw std::logic_error("Native proper parts must be recorded in dense id order.");
        }
        if (owner < 0 || owner >= static_cast<NodeId>(expectedNodes_)) {
            throw std::logic_error("A recorded proper-part owner is outside the node domain.");
        }
        ++nextProperPart_;
    }

    /**
     * @brief Records one owner in a producer-established permutation pass.
     *
     * The producer is responsible for visiting every proper part exactly once.
     * This variant avoids a second dense-domain scan when the construction
     * algorithm already traverses a proven permutation of the proper parts.
     *
     * @param owner Owner associated with the view or node.
     */
    void recordProperPartOwner(NodeId owner) {
        if (nextProperPart_ >= expectedProperParts_) {
            throw std::logic_error("Native proper-part evidence exceeds the declared domain.");
        }
        if (owner < 0 || owner >= static_cast<NodeId>(expectedNodes_)) {
            throw std::logic_error("A recorded proper-part owner is outside the node domain.");
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
    /** @brief Stores the node parent. */
    std::vector<NodeId> nodeParent;
    /** @brief Stores the proper part owner. */
    std::vector<NodeId> properPartOwner;
    /** @brief Stores the altitude. */
    std::vector<T> altitude;
    /** @brief Stores the root. */
    NodeId root = InvalidNode;
    /** @brief Stores the grid domain2 d. */
    std::optional<GridDomain2D> gridDomain2D;
    /** @brief Stores the semantics. */
    HierarchySemantics semantics;
    /** @brief Stores the topology proof. */
    NativeTopologyProof topologyProof;
};

/**
 * @brief Owning, move-only native hierarchy paired with structural evidence.
 *
 * Internal producers move their existing buffers into this representation.
 * Public span imports are validated first and copied once into it. In both
 * cases materialization subsequently transfers the buffers into the weighted
 * owner without another copy.
 */
template <AltitudeValue T> class ValidatedNativeHierarchy {
    /** @brief Stores the storage. */
    ValidatedNativeHierarchyStorage<T> storage_;

  public:
    /**
     * @brief Takes ownership of structurally validated native hierarchy storage.
     *
     * @param storage Storage object that owns the native hierarchy arrays.
     */
    explicit ValidatedNativeHierarchy(ValidatedNativeHierarchyStorage<T>&& storage) : storage_(std::move(storage)) {
        storage_.topologyProof.requireMatches(storage_.nodeParent.size(), storage_.properPartOwner.size(), storage_.root);
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
        if (hierarchy.nodeParent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) ||
            hierarchy.properPartOwner.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native hierarchy domain exceeds the NodeId range.");
        }
        if (hierarchy.nodeParent.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one internal node.");
        }
        if (hierarchy.properPartOwner.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one proper part.");
        }
        if (hierarchy.altitude.size() != hierarchy.nodeParent.size()) {
            throw std::invalid_argument("Native hierarchy altitude size must match the internal-node domain.");
        }

        const NodeId numNodes = static_cast<NodeId>(hierarchy.nodeParent.size());
        if (hierarchy.root < 0 || hierarchy.root >= numNodes) {
            throw std::invalid_argument("Native hierarchy requires a valid root node id.");
        }
        if (hierarchy.gridDomain2D && hierarchy.properPartOwner.size() != hierarchy.gridDomain2D->size("Native hierarchy 2D proper-part domain")) {
            throw std::invalid_argument("Native hierarchy proper-part domain must match the attached 2D grid.");
        }

        switch (hierarchy.semantics.descriptiveKind) {
        case MorphologicalTreeKind::GENERIC:
        case MorphologicalTreeKind::MAX_TREE:
        case MorphologicalTreeKind::MIN_TREE:
        case MorphologicalTreeKind::TREE_OF_SHAPES:
        case MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE:
            break;
        default:
            throw std::invalid_argument("Native hierarchy kind is not supported.");
        }

        const auto* uniformAdjacency = std::get_if<UniformGridAdjacency2D>(&hierarchy.semantics.adjacency);
        if (uniformAdjacency && !hierarchy.gridDomain2D) {
            throw std::invalid_argument("Native hierarchy grid adjacency requires an attached 2D grid domain.");
        }
        if (uniformAdjacency && (uniformAdjacency->relation.getNumRows() != hierarchy.gridDomain2D->rows ||
                                 uniformAdjacency->relation.getNumCols() != hierarchy.gridDomain2D->cols)) {
            throw std::invalid_argument("Native hierarchy adjacency must match the attached 2D grid.");
        }

        const auto* directionalAdjacency = std::get_if<DirectionalGridAdjacency2D>(&hierarchy.semantics.adjacency);
        if (directionalAdjacency && !hierarchy.gridDomain2D) {
            throw std::invalid_argument("Native hierarchy directional grid adjacency requires an attached 2D grid domain.");
        }
        if (directionalAdjacency && (directionalAdjacency->decreasing.getNumRows() != hierarchy.gridDomain2D->rows ||
                                     directionalAdjacency->decreasing.getNumCols() != hierarchy.gridDomain2D->cols ||
                                     directionalAdjacency->increasing.getNumRows() != hierarchy.gridDomain2D->rows ||
                                     directionalAdjacency->increasing.getNumCols() != hierarchy.gridDomain2D->cols)) {
            throw std::invalid_argument("Native hierarchy directional adjacency context must match the attached 2D grid.");
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
     * @param hierarchy Hierarchy data used by the operation.
     * @return Result of performs the complete structural validation required for a public raw-buffer import.
     */
    template <AltitudeValue T> [[nodiscard]] static NativeTopologyProof validateComplete(const NativeHierarchyView<T>& hierarchy) {
        validateDescriptor(hierarchy);

        const NodeId numNodes = static_cast<NodeId>(hierarchy.nodeParent.size());
        std::vector<int> remainingChildren(hierarchy.nodeParent.size(), 0);

        for (NodeId node = 0; node < numNodes; ++node) {
            const NodeId parent = hierarchy.nodeParent[static_cast<std::size_t>(node)];
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

        std::vector<std::uint8_t> subtreeHasSupport(hierarchy.nodeParent.size(), std::uint8_t{0});
        for (NodeId properPart = 0; properPart < static_cast<NodeId>(hierarchy.properPartOwner.size()); ++properPart) {
            const NodeId owner = hierarchy.properPartOwner[static_cast<std::size_t>(properPart)];
            if (owner < 0 || owner >= numNodes) {
                throw std::invalid_argument("Native hierarchy proper-part owner is outside the internal-node domain.");
            }
            subtreeHasSupport[static_cast<std::size_t>(owner)] = 1;
        }

        std::vector<NodeId> ready;
        ready.reserve(hierarchy.nodeParent.size());
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

            const NodeId parent = hierarchy.nodeParent[static_cast<std::size_t>(node)];
            subtreeHasSupport[static_cast<std::size_t>(parent)] = 1;
            int& parentRemaining = remainingChildren[static_cast<std::size_t>(parent)];
            --parentRemaining;
            if (parentRemaining == 0) {
                ready.push_back(parent);
            }
        }
        if (processedNodes != hierarchy.nodeParent.size()) {
            throw std::invalid_argument("Native hierarchy contains a cycle or nodes disconnected from the root.");
        }

        return NativeTopologyProof(hierarchy.nodeParent.size(), hierarchy.properPartOwner.size(), hierarchy.root);
    }

    /**
     * @brief Validates parent connectivity and ownership while consuming
     * support evidence already established by the producer.
     *
     * `validateOwner` runs inside the existing proper-part pass, allowing an
     * native producer to retain its reconstruction oracle without another scan.
     *
     * @param hierarchy Hierarchy data used by the operation.
     * @param supportProof Proof that support ownership is valid.
     * @param validateOwner Whether proper-part ownership must be revalidated.
     * @return Validation result for parent connectivity and ownership while consuming support evidence already established by the producer.
     */
    template <AltitudeValue T, class OwnerValidator>
    [[nodiscard]] static NativeTopologyProof validateWithEstablishedSupport(const NativeHierarchyView<T>& hierarchy, NativeSubtreeSupportProof&& supportProof,
                                                                            OwnerValidator&& validateOwner) {
        validateDescriptor(hierarchy);
        if (!supportProof.valid_ || supportProof.numSupportedNodes_ != hierarchy.nodeParent.size()) {
            throw std::logic_error("Native subtree-support proof does not match the node domain.");
        }

        const NodeId numNodes = static_cast<NodeId>(hierarchy.nodeParent.size());
        if (hierarchy.nodeParent[static_cast<std::size_t>(hierarchy.root)] != hierarchy.root) {
            throw std::invalid_argument("Native hierarchy root must be self-parented.");
        }

        std::vector<std::uint8_t> parentState(hierarchy.nodeParent.size(), std::uint8_t{0});
        parentState[static_cast<std::size_t>(hierarchy.root)] = 2;
        for (NodeId node = 0; node < numNodes; ++node) {
            if (node == hierarchy.root) {
                continue;
            }

            NodeId cursor = node;
            while (parentState[static_cast<std::size_t>(cursor)] == 0) {
                parentState[static_cast<std::size_t>(cursor)] = 1;
                const NodeId parent = hierarchy.nodeParent[static_cast<std::size_t>(cursor)];
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
                cursor = hierarchy.nodeParent[static_cast<std::size_t>(cursor)];
            }
        }

        for (NodeId properPart = 0; properPart < static_cast<NodeId>(hierarchy.properPartOwner.size()); ++properPart) {
            const NodeId owner = hierarchy.properPartOwner[static_cast<std::size_t>(properPart)];
            if (owner < 0 || owner >= numNodes) {
                throw std::invalid_argument("Native hierarchy proper-part owner is outside the internal-node domain.");
            }
            validateOwner(properPart, owner);
        }

        return NativeTopologyProof(hierarchy.nodeParent.size(), hierarchy.properPartOwner.size(), hierarchy.root);
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
        ValidatedNativeHierarchyStorage<T>{std::vector<NodeId>(hierarchy.nodeParent.begin(), hierarchy.nodeParent.end()),
                                           std::vector<NodeId>(hierarchy.properPartOwner.begin(), hierarchy.properPartOwner.end()),
                                           std::vector<T>(hierarchy.altitude.begin(), hierarchy.altitude.end()), hierarchy.root, hierarchy.gridDomain2D,
                                           std::move(hierarchy.semantics), std::move(proof)});
}

template<AltitudeValue T>
[[nodiscard]] ValidatedNativeHierarchy<T>
/**
 * @brief Creates validated native hierarchy.
 *
 * @param nodeParent Parent-related value represented by `nodeParent`.
 * @param properPartOwner Owner-related value represented by `properPartOwner`.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param root Root node of the operation.
 * @param gridDomain2D Two-dimensional grid metadata associated with the hierarchy.
 * @param semantics Semantic metadata associated with the hierarchy.
 * @param proof Validation proof to attach to the imported hierarchy.
 * @return Created validated native hierarchy.
 */
makeValidatedNativeHierarchy(
    std::vector<NodeId>&& nodeParent,
    std::vector<NodeId>&& properPartOwner,
    std::vector<T>&& altitude,
    NodeId root,
    std::optional<GridDomain2D> gridDomain2D,
    HierarchySemantics semantics,
    NativeTopologyProof&& proof) {
    const NativeHierarchyView<T> hierarchy{nodeParent, properPartOwner, altitude, root, gridDomain2D, semantics};
    NativeHierarchyValidation::validateDescriptorOnly(hierarchy);
    proof.requireMatches(nodeParent.size(), properPartOwner.size(), root);

    return ValidatedNativeHierarchy<T>(ValidatedNativeHierarchyStorage<T>{std::move(nodeParent), std::move(properPartOwner), std::move(altitude), root,
                                                                          gridDomain2D, std::move(semantics), std::move(proof)});
}

} // namespace mmcfilters::detail
