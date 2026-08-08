#pragma once

/**
 * @file UnionFindResidualTreeAssembler.hpp
 * @brief Sequential replay of union-find contraction events into tree buffers.
 */

#include "UnionFindRegionTypes.hpp"
#include "../../detail/NativeHierarchyValidationDetail.hpp"
#include "../../../utils/Common.hpp"

#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Replayable residual-tree assembler independent of region adjacency.
 *
 * Each initial union-find region owns an intrusive list of pixels whose proper
 * part has not yet been emitted and a list of residual roots that still need a
 * parent. A contraction concatenates these lists in constant time per absorbed
 * region; it neither inspects nor mutates the implicit contact-edge structure.
 *
 * @tparam T Residual-node altitude type.
 */
template <AltitudeValue T> class UnionFindResidualTreeAssembler {
  public:
    /** @brief Fully materialized native buffers returned at finalization. */
    struct Output {
        /// Parent indexed by residual node id; the root is self-parented.
        std::vector<NodeId> nodeParent;
        /// Direct residual-node owner indexed by row-major pixel id.
        std::vector<NodeId> properPartOwner;
        /// Valuation indexed by residual node id.
        std::vector<T> altitude;
        /// Inductive evidence that every emitted node has subtree support.
        mmcfilters::detail::NativeSubtreeSupportProof subtreeSupportProof;
    };

    /**
     * @brief Initializes one intrusive pixel list per initial flat zone.
     * @param numRegions Size of the dense initial region-id domain.
     * @param regionByPixel Initial region id indexed by row-major pixel id.
     * @throws std::length_error if ids cannot represent the supplied domain.
     * @throws std::out_of_range if a pixel refers to an invalid region.
     */
    UnionFindResidualTreeAssembler(std::size_t numRegions, const std::vector<RegionId>& regionByPixel)
        : regions_(numRegions), properPartOwner_(regionByPixel.size(), InvalidNode), nextPixel_(regionByPixel.size(), InvalidNode), nodeParent_(1, InvalidNode),
          altitude_(1, T{}), nextOpenRoot_(1, InvalidNode) {
        if (numRegions > static_cast<std::size_t>(std::numeric_limits<RegionId>::max())) {
            throw std::length_error("Too many regions for the union-find region-id type.");
        }
        if (regionByPixel.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::length_error("Too many pixels for the residual-tree node-id type.");
        }

        nodeParent_.reserve(numRegions);
        altitude_.reserve(numRegions);
        nextOpenRoot_.reserve(numRegions);
        for (NodeId pixel = 0; pixel < static_cast<NodeId>(regionByPixel.size()); ++pixel) {
            const RegionId region = regionByPixel[static_cast<std::size_t>(pixel)];
            requireRegion(region);
            appendPixel(regions_[static_cast<std::size_t>(region)], pixel);
        }
    }

    /**
     * @brief Emits a residual node for all unresolved content of `region`.
     *
     * Existing open residual roots become children of the new event, and
     * previously unowned pixels become its proper parts.
     *
     * @param region Region identifier used by the operation.
     * @param eventAltitude Altitude or level represented by `eventAltitude`.
     * @return The newly allocated residual node id.
     *
     */
    [[nodiscard]] NodeId emitEvent(RegionId region, T eventAltitude) {
        requireRegion(region);
        if (nodeParent_.size() >= static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::length_error("Too many residual events for the node-id type.");
        }

        RegionLists& lists = regions_[static_cast<std::size_t>(region)];
        if (lists.freshPixelHead == InvalidNode && lists.openRootHead == InvalidNode) {
            throw std::runtime_error("A union-find residual event cannot have empty subtree support.");
        }
        subtreeSupportRecorder_.recordSupportedNode();
        const NodeId event = static_cast<NodeId>(nodeParent_.size());
        nodeParent_.push_back(InvalidNode);
        altitude_.push_back(eventAltitude);
        nextOpenRoot_.push_back(InvalidNode);

        for (NodeId root = lists.openRootHead; root != InvalidNode;) {
            const NodeId next = nextOpenRoot_[static_cast<std::size_t>(root)];
            if (nodeParent_[static_cast<std::size_t>(root)] != InvalidNode) {
                throw std::runtime_error("A union-find residual node received more than one parent.");
            }
            nodeParent_[static_cast<std::size_t>(root)] = event;
            root = next;
        }
        for (NodeId pixel = lists.freshPixelHead; pixel != InvalidNode;) {
            const NodeId next = nextPixel_[static_cast<std::size_t>(pixel)];
            if (properPartOwner_[static_cast<std::size_t>(pixel)] != InvalidNode) {
                throw std::runtime_error("A proper part received more than one direct owner.");
            }
            properPartOwner_[static_cast<std::size_t>(pixel)] = event;
            pixel = next;
        }

        lists.freshPixelHead = InvalidNode;
        lists.freshPixelTail = InvalidNode;
        lists.openRootHead = event;
        lists.openRootTail = event;
        return event;
    }

    /**
     * @brief Returns the number of non-root residual events emitted so far.
     *
     * @return Number of emitted residual events, excluding the reserved root.
     */
    [[nodiscard]] std::size_t numEvents() const noexcept { return nodeParent_.size() - 1; }

    /**
     * @brief Transfers unresolved lists described by a replayed event range.
     * @param survivorRegion Dense root retained by the original contraction.
     * @param absorbedRegions Former live roots absorbed by that contraction.
     */
    void consume(RegionId survivorRegion, std::span<const RegionId> absorbedRegions) {
        requireRegion(survivorRegion);
        RegionLists& survivor = regions_[static_cast<std::size_t>(survivorRegion)];
        for (RegionId absorbed : absorbedRegions) {
            requireRegion(absorbed);
            if (absorbed == survivorRegion) {
                throw std::invalid_argument("A contraction cannot absorb its survivor.");
            }
            RegionLists& source = regions_[static_cast<std::size_t>(absorbed)];
            concatenatePixels(survivor, source);
            concatenateOpenRoots(survivor, source);
            source = RegionLists{};
        }
    }

    /**
     * @brief Attaches remaining content to node 0 and releases output buffers.
     * @param terminalRegion Sole live quotient-partition root.
     * @param terminalAltitude Valuation assigned to root node 0.
     * @return Complete parent, proper-part owner, and altitude buffers.
     */
    [[nodiscard]] Output finalize(RegionId terminalRegion, T terminalAltitude) {
        requireRegion(terminalRegion);
        RegionLists& terminal = regions_[static_cast<std::size_t>(terminalRegion)];
        if (terminal.freshPixelHead == InvalidNode && terminal.openRootHead == InvalidNode) {
            throw std::runtime_error("The union-find residual root cannot have empty subtree support.");
        }
        subtreeSupportRecorder_.recordSupportedNode();
        altitude_[0] = terminalAltitude;
        nodeParent_[0] = 0;

        for (NodeId root = terminal.openRootHead; root != InvalidNode;) {
            const NodeId next = nextOpenRoot_[static_cast<std::size_t>(root)];
            if (nodeParent_[static_cast<std::size_t>(root)] != InvalidNode) {
                throw std::runtime_error("A terminal union-find residual root already has a parent.");
            }
            nodeParent_[static_cast<std::size_t>(root)] = 0;
            root = next;
        }
        for (NodeId pixel = terminal.freshPixelHead; pixel != InvalidNode;) {
            const NodeId next = nextPixel_[static_cast<std::size_t>(pixel)];
            if (properPartOwner_[static_cast<std::size_t>(pixel)] != InvalidNode) {
                throw std::runtime_error("A terminal proper part already has an owner.");
            }
            properPartOwner_[static_cast<std::size_t>(pixel)] = 0;
            pixel = next;
        }

        for (NodeId parent : nodeParent_) {
            if (parent == InvalidNode) {
                throw std::runtime_error("Union-find SDRT assembly left an unparented residual node.");
            }
        }
        for (NodeId owner : properPartOwner_) {
            if (owner == InvalidNode) {
                throw std::runtime_error("Union-find SDRT assembly left an unowned proper part.");
            }
        }

        auto subtreeSupportProof = std::move(subtreeSupportRecorder_).finish(nodeParent_.size());
        return Output{std::move(nodeParent_), std::move(properPartOwner_), std::move(altitude_), std::move(subtreeSupportProof)};
    }

  private:
    /** @brief Intrusive-list endpoints owned by one dense region slot. */
    struct RegionLists {
        /// First pixel whose proper-part owner has not been emitted.
        NodeId freshPixelHead = InvalidNode;
        /// Last pixel in the fresh-pixel list.
        NodeId freshPixelTail = InvalidNode;
        /// First emitted residual root still missing a parent.
        NodeId openRootHead = InvalidNode;
        /// Last residual root in the open-root list.
        NodeId openRootTail = InvalidNode;
    };

    std::vector<RegionLists> regions_;                                        ///< List endpoints indexed by initial region id.
    std::vector<NodeId> properPartOwner_;                                     ///< Output owner indexed by pixel id.
    std::vector<NodeId> nextPixel_;                                           ///< Intrusive next link for fresh pixels.
    std::vector<NodeId> nodeParent_;                                          ///< Output parent indexed by node id.
    std::vector<T> altitude_;                                                 ///< Output valuation indexed by node id.
    std::vector<NodeId> nextOpenRoot_;                                        ///< Intrusive next link for open roots.
    mmcfilters::detail::NativeSubtreeSupportRecorder subtreeSupportRecorder_; ///< O(1)-per-node support evidence.

    /**
     * @brief Rejects a region outside the assembler's dense domain.
     *
     * @param region Region identifier used by the operation.
     */
    void requireRegion(RegionId region) const {
        if (region < 0 || region >= static_cast<RegionId>(regions_.size())) {
            throw std::out_of_range("Union-find residual assembler region id is outside its dense domain.");
        }
    }

    /**
     * @brief Appends one pixel to a region's fresh-pixel list.
     *
     * @param lists Per-node lists used by the operation.
     * @param pixel Pixel identifier used by the operation.
     */
    void appendPixel(RegionLists& lists, NodeId pixel) {
        if (lists.freshPixelHead == InvalidNode) {
            lists.freshPixelHead = pixel;
        } else {
            nextPixel_[static_cast<std::size_t>(lists.freshPixelTail)] = pixel;
        }
        lists.freshPixelTail = pixel;
    }

    /**
     * @brief Concatenates fresh-pixel lists in constant time.
     *
     * @param destination Destination represented by `destination`.
     * @param source Source object or value.
     */
    void concatenatePixels(RegionLists& destination, const RegionLists& source) {
        if (source.freshPixelHead == InvalidNode) {
            return;
        }
        if (destination.freshPixelHead == InvalidNode) {
            destination.freshPixelHead = source.freshPixelHead;
        } else {
            nextPixel_[static_cast<std::size_t>(destination.freshPixelTail)] = source.freshPixelHead;
        }
        destination.freshPixelTail = source.freshPixelTail;
    }

    /**
     * @brief Concatenates unresolved residual-root lists in constant time.
     *
     * @param destination Destination represented by `destination`.
     * @param source Source object or value.
     */
    void concatenateOpenRoots(RegionLists& destination, const RegionLists& source) {
        if (source.openRootHead == InvalidNode) {
            return;
        }
        if (destination.openRootHead == InvalidNode) {
            destination.openRootHead = source.openRootHead;
        } else {
            nextOpenRoot_[static_cast<std::size_t>(destination.openRootTail)] = source.openRootHead;
        }
        destination.openRootTail = source.openRootTail;
    }
};

} // namespace mmcfilters::sdrt::detail
