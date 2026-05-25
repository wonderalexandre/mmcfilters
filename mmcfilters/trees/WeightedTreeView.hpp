#pragma once

#include "MorphologicalTree.hpp"
#include "../utils/Altitude.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>

namespace mmcfilters {

/**
 * @brief Non-owning view pairing a topology with an external altitude span.
 *
 * `WeightedTreeView` is the lightweight C++ representation of the Higra-like
 * model used by generic-altitude kernels: topology ownership stays elsewhere,
 * and altitude is borrowed as a contiguous read-only span indexed by the dense
 * internal `NodeId` domain. The caller must keep both the topology and the
 * altitude storage alive for the lifetime of the view. The view intentionally
 * carries no ownership or canonical-tree identity; operations that need owned
 * weighted state must receive a `WeightedMorphologicalTree<T>` instead.
 */
template<AltitudeValue T>
class WeightedTreeView {
private:
    const MorphologicalTree* topology_ = nullptr;
    AltitudeSpan<T> altitude_;
    std::size_t topologyMutationVersion_ = 0;

    /**
     * @brief Validates that the altitude span covers every internal node slot.
     */
    void validateAltitudeBufferShape() const {
        if (altitude_.size() != static_cast<std::size_t>(topology_->getNumInternalNodeSlots())) {
            throw std::runtime_error("WeightedTreeView altitude size must match the dense internal-node domain.");
        }
    }

public:
    /// Altitude scalar type borrowed by this view.
    using altitude_type = T;

    /**
     * @brief Builds a view over caller-owned topology and altitude storage.
     *
     * The span must be indexed by the dense internal `NodeId` domain and remain
     * alive for the full lifetime of the view.
     */
    explicit WeightedTreeView(const MorphologicalTree& topology, AltitudeSpan<T> altitude)
        : topology_(&topology), altitude_(altitude), topologyMutationVersion_(topology.getMutationVersion()) {
        validateAltitudeBufferShape();
    }

    /**
     * @brief Builds a view over a caller-owned altitude buffer.
     *
     * This overload is equivalent to passing `std::span<const T>(altitude)` and
     * exists so `WeightedTreeView(tree, buffer)` deduces the altitude type.
     */
    explicit WeightedTreeView(const MorphologicalTree& topology, const AltitudeBuffer<T>& altitude)
        : topology_(&topology), altitude_(std::span<const T>(altitude)), topologyMutationVersion_(topology.getMutationVersion()) {
        validateAltitudeBufferShape();
    }

    /**
     * @brief Rejects temporary altitude buffers that would leave a dangling span.
     */
    WeightedTreeView(const MorphologicalTree&, AltitudeBuffer<T>&&) = delete;

    /**
     * @brief Returns the borrowed tree topology.
     */
    [[nodiscard]] const MorphologicalTree& topology() const noexcept {
        return *topology_;
    }

    /**
     * @brief Returns the borrowed altitude span indexed by internal `NodeId`.
     */
    [[nodiscard]] AltitudeSpan<T> altitude() const noexcept {
        return altitude_;
    }

    /**
     * @brief Throws if the borrowed topology changed since view construction.
     */
    void requireTopologyUnchanged(const char* context) const {
        topology().requireMutationVersion(topologyMutationVersion_, context);
    }

    /**
     * @brief Returns the altitude associated with an internal node id.
     */
    [[nodiscard]] T getAltitude(NodeId nodeId) const {
        if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedTreeView::getAltitude requires a valid internal NodeId.");
        }
        return altitude_[static_cast<std::size_t>(nodeId)];
    }

    /**
     * @brief Returns the altitude difference between a live node and its parent.
     *
     * Root-like nodes, whose parent is invalid or themselves, use their own
     * altitude as residue.
     */
    [[nodiscard]] AltitudeDiff<T> getNodeResidue(NodeId nodeId) const {
        const MorphologicalTree& tree = topology();
        if (!tree.isAlive(nodeId)) {
            throw std::invalid_argument("WeightedTreeView::getNodeResidue requires a live internal NodeId.");
        }
        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            return static_cast<AltitudeDiff<T>>(getAltitude(nodeId));
        }
        return static_cast<AltitudeDiff<T>>(getAltitude(nodeId)) -
               static_cast<AltitudeDiff<T>>(getAltitude(parentNodeId));
    }
};

} // namespace mmcfilters
