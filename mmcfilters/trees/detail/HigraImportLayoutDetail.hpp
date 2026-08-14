#pragma once

#include "../../utils/Common.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>

namespace mmcfilters::detail {

/**
 * @brief Leaf-first Higra id layout used by the import adapter.
 *
 * Higra stores proper-part leaves in `[0, numLeaves)` and internal nodes in
 * `[numLeaves, numVertices)`. `MorphologicalTree` stores internal nodes in a
 * separate dense domain, so conversion is a constant offset.
 *
 * This helper owns no buffers and performs no traversal. Import loops use it
 * in place, preserving their pass count, allocation behavior, and node-id
 * assignment.
 */
class HigraImportLayout {
  private:
    /** @brief Number of leaves. */
    int numLeaves_ = 0;
    /** @brief Number of vertices. */
    int numVertices_ = 0;
    /** @brief Number of internal nodes. */
    int numInternalNodes_ = 0;

  public:
    /**
     * @brief Constructs `HigraImportLayout` from the supplied inputs.
     *
     * @param numLeaves Count.
     * @param numVertices Count.
     */
    HigraImportLayout(int numLeaves, std::size_t numVertices) : numLeaves_(numLeaves) {
        if (numVertices > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Higra parent array exceeds the supported node-id range.");
        }

        numVertices_ = static_cast<NodeId>(numVertices);
        if (numLeaves_ <= 0 || numLeaves_ >= numVertices_) {
            throw std::invalid_argument("Higra parent array must contain leaves followed by at least one internal node.");
        }
        numInternalNodes_ = numVertices_ - numLeaves_;
    }

    /**
     * @brief Returns the number of leaves in the imported Higra hierarchy.
     *
     * @return Number of leaf vertices.
     */
    [[nodiscard]] int numLeaves() const noexcept { return numLeaves_; }

    /**
     * @brief Returns the total number of hierarchy vertices.
     *
     * @return Total number of hierarchy vertices.
     */
    [[nodiscard]] int numVertices() const noexcept { return numVertices_; }

    /**
     * @brief Returns the number of internal hierarchy nodes.
     *
     * @return Number of leaf vertices.
     */
    [[nodiscard]] int numInternalNodes() const noexcept { return numInternalNodes_; }

    /**
     * @brief Maps one internal Higra id to the dense internal-node domain.
     *
     * @param higraNodeId Node identifier.
     * @param invalidIdMessage Diagnostic message used for invalid identifiers.
     * @return The mapped internal Higra id to the dense internal-node domain.
     */
    [[nodiscard]] NodeId internalNodeId(NodeId higraNodeId, const char* invalidIdMessage) const {
        if (higraNodeId < numLeaves_ || higraNodeId >= numVertices_) {
            throw std::invalid_argument(invalidIdMessage);
        }
        return higraNodeId - numLeaves_;
    }

    /**
     * @brief Maps one dense internal node id to its imported Higra id.
     *
     * @param internalNodeId Node identifier.
     * @return The mapped dense internal node id to its imported Higra id.
     */
    [[nodiscard]] NodeId higraNodeId(NodeId internalNodeId) const noexcept { return numLeaves_ + internalNodeId; }

    /**
     * @brief Tests the self-parent encoding used for the unique Higra root.
     *
     * @param higraNodeId Node identifier.
     * @param parentHigraNodeId Node identifier.
     * @return True when the documented condition holds; otherwise false.
     */
    [[nodiscard]] static bool isRootLink(NodeId higraNodeId, NodeId parentHigraNodeId) noexcept { return higraNodeId == parentHigraNodeId; }
};

} // namespace mmcfilters::detail
