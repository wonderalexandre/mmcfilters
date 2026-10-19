#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeAlgorithms.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Contract.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>

namespace mmcfilters {

template <AltitudeValue T> class ValuedMorphologicalTree;

/**
 * @brief Non-owning view pairing a topology with an external altitude span.
 *
 * `ValuedMorphologicalTreeView` is the lightweight C++ representation of the Higra-like
 * model used by generic-altitude kernels: topology storage stays elsewhere,
 * and node altitudes are borrowed as a contiguous read-only span indexed by the dense
 * internal `NodeId` domain. The caller must keep both the topology and the
 * altitude storage alive and unchanged for the lifetime of the view. Public
 * constructors validate external altitudes against the declared topology order.
 * Views returned by `ValuedMorphologicalTree::asView()` reuse the guarantees
 * maintained by the owner without scanning the altitude buffer. The view intentionally
 * carries no ownership or canonical-tree identity; operations that need owned
 * valued state must receive a `ValuedMorphologicalTree<T>` instead.
 */
template <AltitudeValue T> class ValuedMorphologicalTreeView {
    friend class ValuedMorphologicalTree<T>;

  private:
    /** @brief Topology. */
    const MorphologicalTree* topology_ = nullptr;
    /** @brief Node altitudes. */
    NodeAltitudeSpan<T> nodeAltitudes_;
    /** @brief Topology mutation version used to detect stale derived state. */
    std::size_t topologyMutationVersion_ = 0;

    /**
     * @brief Validates that the altitude span covers every internal node slot.
     */
    void validateNodeAltitudeBufferShape() const {
        MMCFILTERS_CONTRACT_REQUIRE(nodeAltitudes_.size() == static_cast<std::size_t>(topology_->numInternalNodeSlots()),
                                    throw std::runtime_error("ValuedMorphologicalTreeView altitude size must match the dense internal-node domain."));
    }

    /**
     * @brief Borrows the established valued state without repeating validation.
     * @param owner Owner that maintains topology and altitude invariants.
     */
    explicit ValuedMorphologicalTreeView(const ValuedMorphologicalTree<T>& owner)
        : topology_(&owner.topology()), nodeAltitudes_(owner.nodeAltitudeSpan()), topologyMutationVersion_(owner.topology().getMutationVersion()) {}

  public:
    /// Altitude scalar type borrowed by this view.
    using AltitudeType = T;

    /**
     * @brief Builds a view over caller-owned topology and altitude storage.
     *
     * The span must be indexed by the dense internal `NodeId` domain and remain
     * alive and unchanged for the full lifetime of the view. Construction checks
     * its shape, finite values, and the declared parent-child altitude order.
     *
     * @param topology Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     */
    explicit ValuedMorphologicalTreeView(const MorphologicalTree& topology, NodeAltitudeSpan<T> altitude)
        : topology_(&topology), nodeAltitudes_(altitude), topologyMutationVersion_(topology.getMutationVersion()) {
        validateNodeAltitudeBufferShape();
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(altitude, "ValuedMorphologicalTreeView");
        TreeAltitudeAlgorithms::validateMonotoneNodeAltitudes(topology, altitude);
    }

    /**
     * @brief Builds a view over a caller-owned altitude buffer.
     *
     * This overload is equivalent to passing `std::span<const T>(altitude)` and
     * exists so `ValuedMorphologicalTreeView(tree, buffer)` deduces the altitude type.
     *
     * @param topology Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     */
    explicit ValuedMorphologicalTreeView(const MorphologicalTree& topology, const NodeAltitudeBuffer<T>& altitude)
        : ValuedMorphologicalTreeView(topology, std::span<const T>(altitude)) {}

    /**
     * @brief Rejects temporary altitude buffers that would leave a dangling span.
     */
    ValuedMorphologicalTreeView(const MorphologicalTree&, NodeAltitudeBuffer<T>&&) = delete;

    /**
     * @brief Returns the borrowed tree topology.
     *
     * @return The borrowed tree topology.
     */
    [[nodiscard]] const MorphologicalTree& topology() const noexcept { return *topology_; }

    /**
     * @brief Returns the borrowed altitude span indexed by internal `NodeId`.
     *
     * @return The borrowed altitude span indexed by internal NodeId.
     */
    [[nodiscard]] NodeAltitudeSpan<T> nodeAltitudes() const noexcept { return nodeAltitudes_; }

    /**
     * @brief Throws if the borrowed topology changed since view construction.
     *
     * @param context Operation context or diagnostic label.
     */
    void requireTopologyUnchanged(const char* context) const { topology().requireMutationVersion(topologyMutationVersion_, context); }

    /**
     * @brief Returns the altitude associated with an internal node id.
     *
     * @param nodeId Dense internal node identifier.
     * @return The altitude associated with an internal node id.
     */
    [[nodiscard]] T nodeAltitude(NodeId nodeId) const {
        MMCFILTERS_CONTRACT_REQUIRE(nodeId >= 0 && static_cast<std::size_t>(nodeId) < nodeAltitudes_.size(),
                                    throw std::invalid_argument("ValuedMorphologicalTreeView::nodeAltitude requires a valid internal NodeId."));
        return nodeAltitudes_[static_cast<std::size_t>(nodeId)];
    }

    /**
     * @brief Returns the altitude difference between a live node and its parent.
     *
     * The project adopts a fixed zero reconstruction baseline. Root-like nodes,
     * whose parent is invalid or themselves, therefore use their own altitude
     * as residue.
     *
     * @param nodeId Dense internal node identifier.
     * @return The altitude difference between a live node and its parent.
     */
    [[nodiscard]] AltitudeDifference<T> nodeResidue(NodeId nodeId) const {
        const MorphologicalTree& tree = topology();
        MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId), throw std::invalid_argument("ValuedMorphologicalTreeView::nodeResidue requires a live internal NodeId."));
        const NodeId parentNodeId = tree.parent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            return static_cast<AltitudeDifference<T>>(nodeAltitude(nodeId));
        }
        return static_cast<AltitudeDifference<T>>(nodeAltitude(nodeId)) - static_cast<AltitudeDifference<T>>(nodeAltitude(parentNodeId));
    }
};

} // namespace mmcfilters
