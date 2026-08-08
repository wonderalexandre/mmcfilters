#pragma once

/**
 * @file MinMaxResidualTreeBuilder.hpp
 * @brief Unrestricted or saturated residual hierarchy from synchronized component trees.
 */

#include "SdrtTiePolicy.hpp"
#include "detail/UnionFindResidualTreeAssembler.hpp"
#include "../adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "../detail/NativeHierarchyValidationDetail.hpp"
#include "../../utils/GenerationStampSet.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt {

/**
 * @brief Lowest-common-ancestor strategy used by the synchronized trees.
 */
enum class SaturatedMinMaxLcaPolicy {
    /** Climb current parent paths, with worst-case cost proportional to height. */
    ParentClimb,

    /**
     * Rebuild an Euler/RMQ snapshot periodically and use it outside subtrees
     * touched since the snapshot, falling back to exact parent climbs inside
     * them.
     */
    BlockedSnapshot,

    /** Maintain the rooted forest with link-cut access operations. */
    LinkCut
};

/**
 * @brief Exact complement-connectivity fallback used after an inconclusive
 *        dual-tree certificate.
 */
enum class SaturatedMinMaxFallbackPolicy {
    /** Explore from one boundary pixel until every other boundary pixel is found. */
    SingleSourceDepthFirst,

    /**
     * Grow all boundary components simultaneously and stop when their labels
     * merge. This is exact and is the default optimized policy.
     */
    BoundaryMultiSource
};

/**
 * @brief Strategy used to obtain the current proper-part boundary.
 */
enum class SaturatedMinMaxBoundaryPolicy {
    /** Rebuild the support and its boundary for every certificate. */
    RecomputeFromSupport,

    /**
     * Cache proper parts and boundary incidences, merge complete prunes
     * small-to-large, and rebuild only nodes affected by a splitting dual
     * update.
     */
    IncrementalSmallToLarge
};

/**
 * @brief Eligibility rule applied to current regional extrema.
 *
 * Both modes use the same synchronized max-tree/min-tree implementation.
 * `SaturatedOnly` retains the exterior-complement connectivity test, whereas
 * `AllRegionalExtrema` accepts every current non-root leaf and therefore
 * builds the unrestricted residual tree.
 */
enum class MinMaxResidualEligibilityPolicy { SaturatedOnly, AllRegionalExtrema };

/**
 * @brief Minimal diagnostics for the synchronized min-tree/max-tree build.
 *
 * The counters retained here participate in consistency checks or public
 * regression tests. Detailed profiling belongs in external benchmarks so the
 * production construction remains free of observational overhead.
 */
struct MinMaxResidualTreeBuildStatistics {
    std::size_t residualEvents = 0;                  ///< Emitted non-root residual nodes.
    std::size_t rejectedExtrema = 0;                 ///< Candidates rejected by saturation.
    std::size_t complementTraversalCertificates = 0; ///< Exact fallback certificates.
};

/**
 * @brief Builds a residual tree from two synchronized mutable component trees.
 *
 * Both component trees use the same symmetric adjacency. By default, a
 * non-exterior leaf is eligible only when its current proper-part support has
 * connected complement. `AllRegionalExtrema` disables that restriction and
 * yields the unrestricted residual tree with the same agenda, updater and
 * output assembly. Accepted max-tree leaves are pruned while the min-tree is
 * updated incrementally, and conversely for accepted min-tree leaves.
 *
 * Saturation is certified without materializing a region-adjacency graph. The
 * fast certificate collects the opposite-tree owners touching the candidate
 * boundary and computes their lowest common ancestor. For a max-tree leaf at
 * level `a`, an opposite min-tree ancestor at a level strictly below `a`
 * provides a connected support excluding the leaf; the dual condition applies
 * to min-tree leaves. When this sufficient proof is inconclusive, the default
 * exact fallback grows simultaneously from every external-boundary pixel and
 * stops when their union-find labels merge. A single-source reference policy
 * remains available. Both are exact because every component of the complement
 * of a proper connected set touches its external boundary. The default LCA
 * policy uses altitude-aligned parent climbs. Exact blocked Euler/RMQ and
 * link-cut policies are available for controlled performance studies; none of
 * these policies changes the emitted hierarchy.
 *
 * A rejected leaf is removed from the agenda until a dual-tree edit reports
 * that node as affected. This is exact because contracting a connected region
 * disjoint from a fixed support cannot change the number of connected
 * components in that support's complement.
 *
 * @tparam T Finite image, component-tree, and residual altitude type.
 * @ingroup sdrt
 */
template <AltitudeValue T> class MinMaxResidualTreeBuilder {
  private:
    /** @brief Defines the `altitude_t` alias used by the component. */
    using altitude_t = T;
    /** @brief Defines the `image_ptr_t` alias used by the component. */
    using image_ptr_t = ImagePtr<altitude_t>;
    /** @brief Defines the `tree_t` alias used by the component. */
    using tree_t = WeightedMorphologicalTree<altitude_t>;
    /** @brief Defines the `Adjustment` alias used by the component. */
    using Adjustment = adjust::DualMinMaxTreeIncrementalFilterLeaf<altitude_t>;

    /** @brief Enumerates the supported polarity values. */
    enum class Polarity : std::uint8_t { Max, Min };

    /** @brief Describes one residual-tree candidate and its deterministic ordering keys. */
    struct Candidate {
        /** @brief Stores the area. */
        int area = 0;
        /** @brief Stores the stable spatial key. */
        NodeId stableSpatialKey = InvalidNode;
        /** @brief Stores the polarity. */
        Polarity polarity = Polarity::Max;
        /** @brief Stores the node identifier. */
        NodeId nodeId = InvalidNode;
    };

    /** @brief Caches the derived properties of one residual-tree candidate. */
    struct CandidateDescriptor {
        /** @brief Stores the area. */
        int area = 0;
        /** @brief Stores the stable spatial key. */
        NodeId stableSpatialKey = InvalidNode;
        /** @brief Indicates whether the candidate support contains the infinity pixel. */
        bool containsInfinity = false;
    };

    /** @brief Orders residual-tree candidates according to the configured tie policy. */
    struct CandidateLess {
        /** @brief Stores the policy. */
        SdrtTiePolicy policy = SdrtTiePolicy::ContrastInvariantSpatial;

        /**
         * @brief Returns the deterministic ordering rank of a component-tree polarity.
         *
         * @param polarity Component-tree polarity to process.
         * @return `0` for max-tree candidates and `1` for min-tree candidates.
         */
        [[nodiscard]] static int polarityRank(Polarity polarity) noexcept { return polarity == Polarity::Max ? 0 : 1; }

        /**
         * @brief Compares two candidates according to the configured deterministic order.
         *
         * @param lhs Left-hand value of the comparison.
         * @param rhs Right-hand value of the comparison.
         * @return `true` when the left-hand candidate precedes the right-hand candidate.
         */
        [[nodiscard]] bool operator()(const Candidate& lhs, const Candidate& rhs) const noexcept {
            if (lhs.area != rhs.area) {
                return lhs.area < rhs.area;
            }
            if (policy == SdrtTiePolicy::MaxBeforeMinThenSpatial && lhs.polarity != rhs.polarity) {
                return polarityRank(lhs.polarity) < polarityRank(rhs.polarity);
            }
            if (lhs.stableSpatialKey != rhs.stableSpatialKey) {
                return lhs.stableSpatialKey < rhs.stableSpatialKey;
            }
            if (lhs.polarity != rhs.polarity) {
                return polarityRank(lhs.polarity) < polarityRank(rhs.polarity);
            }
            return lhs.nodeId < rhs.nodeId;
        }
    };

    /** @brief Stores the active agenda entry associated with one component-tree node. */
    struct AgendaSlot {
        /** @brief Indicates whether the slot currently owns an agenda candidate. */
        bool active = false;
        /** @brief Stores the candidate. */
        Candidate candidate;
    };

    /** @brief Maintains the deterministic agenda of eligible component-tree leaves. */
    class SaturatedLeafAgenda {
      private:
        /** @brief Stores the policy. */
        SdrtTiePolicy policy_;
        /** @brief Stores the infinity pixel. */
        NodeId infinityPixel_ = InvalidNode;
        /** @brief Indicates whether the infinity pixel must be excluded from candidates. */
        bool excludeInfinity_ = true;
        /** @brief Stores the candidates. */
        std::set<Candidate, CandidateLess> candidates_;
        /** @brief Stores the max slots. */
        std::vector<AgendaSlot> maxSlots_;
        /** @brief Stores the min slots. */
        std::vector<AgendaSlot> minSlots_;

        /**
         * @brief Returns the agenda-slot array for one component-tree polarity.
         *
         * @param polarity Component-tree polarity to process.
         * @return A reference to the requested object.
         */
        [[nodiscard]] std::vector<AgendaSlot>& slots(Polarity polarity) { return polarity == Polarity::Max ? maxSlots_ : minSlots_; }

        /**
         * @brief Returns the agenda-slot array for one component-tree polarity.
         *
         * @param polarity Component-tree polarity to process.
         * @return A reference to the requested object.
         */
        [[nodiscard]] const std::vector<AgendaSlot>& slots(Polarity polarity) const { return polarity == Polarity::Max ? maxSlots_ : minSlots_; }

        /**
         * @brief Checks whether a node identifier indexes an agenda-slot array.
         *
         * @param side Agenda-slot array to validate.
         * @param nodeId Identifier of the node processed by the operation.
         * @return `true` when the node identifier is within the slot array; otherwise `false`.
         */
        [[nodiscard]] static bool validSlot(const std::vector<AgendaSlot>& side, NodeId nodeId) noexcept {
            return nodeId >= 0 && nodeId < static_cast<NodeId>(side.size());
        }

        /**
         * @brief Removes one node from the active candidate agenda.
         *
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         */
        void remove(Polarity polarity, NodeId nodeId) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                return;
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            if (!slot.active) {
                return;
            }
            if (candidates_.erase(slot.candidate) != 1) {
                throw std::runtime_error("Min/max residual agenda violates its slot/set invariant.");
            }
            slot = AgendaSlot{};
        }

        /**
         * @brief Indexes all current regional extrema of one component tree.
         *
         * @param tree Tree processed by the operation.
         * @param polarity Component-tree polarity to process.
         */
        void indexTree(const tree_t& tree, Polarity polarity) {
            for (NodeId nodeId : tree.topology().getAliveNodeIds()) {
                update(tree, polarity, nodeId);
            }
        }

      public:
        /**
         * @brief Constructs a `SaturatedLeafAgenda` instance.
         *
         * @param policy Policy controlling the operation.
         * @param infinityPixel Row-major pixel used as the exterior seed.
         * @param excludeInfinity Whether the infinity pixel must be excluded from the agenda.
         */
        SaturatedLeafAgenda(SdrtTiePolicy policy, NodeId infinityPixel, bool excludeInfinity)
            : policy_(policy), infinityPixel_(infinityPixel), excludeInfinity_(excludeInfinity), candidates_(CandidateLess{policy}) {}

        /**
         * @brief Initializes the data structure from the current component trees.
         *
         * @param maxTree Max-tree consumed by the operation.
         * @param minTree Min-tree consumed by the operation.
         */
        void initialize(const tree_t& maxTree, const tree_t& minTree) {
            candidates_ = std::set<Candidate, CandidateLess>(CandidateLess{policy_});
            maxSlots_.assign(static_cast<std::size_t>(maxTree.topology().getNumInternalNodeSlots()), AgendaSlot{});
            minSlots_.assign(static_cast<std::size_t>(minTree.topology().getNumInternalNodeSlots()), AgendaSlot{});
            indexTree(maxTree, Polarity::Max);
            indexTree(minTree, Polarity::Min);
        }

        /**
         * @brief Refreshes the active agenda entry of one component-tree node.
         *
         * @param tree Tree processed by the operation.
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         */
        void update(const tree_t& tree, Polarity polarity, NodeId nodeId) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                return;
            }
            remove(polarity, nodeId);

            const MorphologicalTree& topology = tree.topology();
            if (!topology.isNode(nodeId) || !topology.isAlive(nodeId) || nodeId == topology.getRoot() || !topology.isLeaf(nodeId)) {
                return;
            }
            const NodeId parent = topology.getNodeParent(nodeId);
            if (parent == InvalidNode || parent == nodeId) {
                throw std::runtime_error("Min/max residual leaf violates the parent invariant.");
            }
            const int area = topology.getNumProperParts(nodeId);
            if (area <= 0) {
                return;
            }
            NodeId spatial = std::numeric_limits<NodeId>::max();
            bool containsInfinity = false;
            for (NodeId pixel : topology.getProperParts(nodeId)) {
                spatial = std::min(spatial, pixel);
                containsInfinity = containsInfinity || pixel == infinityPixel_;
            }
            if (spatial == std::numeric_limits<NodeId>::max()) {
                throw std::runtime_error("Min/max residual candidate has an empty proper part.");
            }
            if (excludeInfinity_ && containsInfinity) {
                return;
            }

            Candidate candidate{area, spatial, polarity, nodeId};
            const auto [_, inserted] = candidates_.insert(candidate);
            if (!inserted) {
                throw std::runtime_error("Min/max residual agenda received a duplicate candidate.");
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            slot.active = true;
            slot.candidate = candidate;
        }

        /**
         * @brief Updates an agenda entry from a precomputed candidate descriptor.
         *
         * @param tree Tree processed by the operation.
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         * @param descriptor Descriptor associated with the candidate.
         */
        void updateKnown(const tree_t& tree, Polarity polarity, NodeId nodeId, const CandidateDescriptor& descriptor) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                return;
            }
            remove(polarity, nodeId);

            const MorphologicalTree& topology = tree.topology();
            if (!topology.isNode(nodeId) || !topology.isAlive(nodeId) || nodeId == topology.getRoot() || !topology.isLeaf(nodeId)) {
                return;
            }
            if (descriptor.area <= 0 || descriptor.stableSpatialKey == InvalidNode || topology.getNumProperParts(nodeId) != descriptor.area) {
                throw std::runtime_error("Min/max residual flat-zone metadata violates the agenda invariant.");
            }
            if (excludeInfinity_ && descriptor.containsInfinity) {
                return;
            }

            Candidate candidate{descriptor.area, descriptor.stableSpatialKey, polarity, nodeId};
            const auto [_, inserted] = candidates_.insert(candidate);
            if (!inserted) {
                throw std::runtime_error("Min/max residual agenda received a duplicate metadata candidate.");
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            slot.active = true;
            slot.candidate = candidate;
        }

        /**
         * @brief Removes a rejected node from the active candidate agenda.
         *
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         */
        void reject(Polarity polarity, NodeId nodeId) { remove(polarity, nodeId); }

        /**
         * @brief Selects the first candidate in deterministic agenda order.
         *
         * @return The next candidate, or `std::nullopt` when the agenda is empty.
         */
        [[nodiscard]] std::optional<Candidate> select() const {
            if (candidates_.empty()) {
                return std::nullopt;
            }
            return *candidates_.begin();
        }

        /**
         * @brief Checks whether a node currently has an active agenda entry.
         *
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         * @return `true` when the node has an active agenda entry; otherwise `false`.
         */
        [[nodiscard]] bool contains(Polarity polarity, NodeId nodeId) const {
            const auto& side = slots(polarity);
            return validSlot(side, nodeId) && side[static_cast<std::size_t>(nodeId)].active;
        }
    };

    /** @brief Stores one immutable link in a persistent residual-event stack. */
    struct PersistentStackNode {
        /** @brief Stores the event identifier. */
        int eventId = -1;
        /** @brief Stores the prev. */
        int prev = -1;
    };

    /** @brief Stores the persistent residual-event chains associated with image pixels. */
    struct PersistentChains {
        /** @brief Stores the stack head by pixel. */
        std::vector<int> stackHeadByPixel;
        /** @brief Stores the nodes. */
        std::vector<PersistentStackNode> nodes;
        /** @brief Stores the event valuations. */
        std::vector<altitude_t> eventValuations;
    };

    /**
     * @brief Exact batched dynamic LCA index.
     *
     * At the beginning of each block, an Euler tour and a sparse RMQ table are
     * built in `O(m log m)`. Net parent changes invalidate snapshot subtrees by
     * preorder intervals. A query whose two endpoints remain outside those
     * intervals is answered in `O(1)` on the snapshot; all other queries use the
     * exact current parent paths. Rebuilding occurs after
     * `Theta(sqrt(m))` net topology-changing steps.
     */
    class BlockedDynamicLcaIndex {
      private:
        /** @brief Stores one frame of an iterative tree traversal. */
        struct TraversalFrame {
            /** @brief Stores the node identifier. */
            NodeId nodeId = InvalidNode;
            /** @brief Stores the next child. */
            NodeId nextChild = InvalidNode;
            /** @brief Stores the depth. */
            int depth = 0;
        };

        /** @brief References the tree used by the component. */
        const tree_t* tree_ = nullptr;
        /** @brief Indicates whether the blocked snapshot index is enabled. */
        bool enabled_ = false;
        /** @brief Indicates whether the blocked snapshot requires a rebuild. */
        bool rebuildPending_ = true;
        /** @brief Stores the mutations since rebuild. */
        std::size_t mutationsSinceRebuild_ = 0;
        /** @brief Stores the block size. */
        std::size_t blockSize_ = 1;
        /** @brief Stores the snapshot alive nodes. */
        std::size_t snapshotAliveNodes_ = 0;
        /** @brief Stores the euler. */
        std::vector<NodeId> euler_;
        /** @brief Stores the euler depth. */
        std::vector<int> eulerDepth_;
        /** @brief Stores the first occurrence. */
        std::vector<int> firstOccurrence_;
        /** @brief Stores the preorder entry. */
        std::vector<int> preorderEntry_;
        /** @brief Stores the preorder exit. */
        std::vector<int> preorderExit_;
        /** @brief Stores the rmq log2. */
        std::vector<int> rmqLog2_;
        /** @brief Stores the rmq sparse table. */
        std::vector<int> rmqSparseTable_;
        /** @brief Stores the rmq stride. */
        std::size_t rmqStride_ = 0;
        /** @brief Stores the rmq levels. */
        int rmqLevels_ = 0;
        /** @brief Stores the dirty preorder. */
        std::vector<std::uint8_t> dirtyPreorder_;
        /** @brief Stores the next unmarked. */
        std::vector<int> nextUnmarked_;
        /** @brief Stores the dirty root marks. */
        GenerationStampSet dirtyRootMarks_;
        /** @brief Stores the mirrored parent. */
        std::vector<NodeId> mirroredParent_;
        /** @brief Stores the mirrored alive. */
        std::vector<std::uint8_t> mirroredAlive_;

        /**
         * @brief Selects the shallower of two Euler-tour positions.
         *
         * @param lhs Left-hand value of the comparison.
         * @param rhs Right-hand value of the comparison.
         * @return The Euler-tour index whose occurrence has minimum depth.
         */
        [[nodiscard]] int betterEulerIndex(int lhs, int rhs) const noexcept {
            if (lhs < 0) {
                return rhs;
            }
            if (rhs < 0) {
                return lhs;
            }
            return eulerDepth_[static_cast<std::size_t>(lhs)] <= eulerDepth_[static_cast<std::size_t>(rhs)] ? lhs : rhs;
        }

        /**
         * @brief Appends one node occurrence to the Euler-tour snapshot.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @param depth Tree depth of the appended occurrence.
         */
        void appendEuler(NodeId nodeId, int depth) {
            euler_.push_back(nodeId);
            eulerDepth_.push_back(depth);
        }

        /**
         * @brief Checks whether a node occurrence belongs to a dirty snapshot interval.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @return `true` when the node belongs to a dirty snapshot interval; otherwise `false`.
         */
        [[nodiscard]] bool snapshotNodeIsDirty(NodeId nodeId) const {
            if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= preorderEntry_.size()) {
                return true;
            }
            const int entry = preorderEntry_[static_cast<std::size_t>(nodeId)];
            if (entry < 0) {
                return true;
            }
            return dirtyPreorder_[static_cast<std::size_t>(entry)] != 0;
        }

        /**
         * @brief Finds the next Euler-tour position not covered by a dirty interval.
         *
         * @param position First Euler-tour position to inspect.
         * @return The matching node identifier, or the operation-specific sentinel when absent.
         */
        [[nodiscard]] int findNextUnmarked(int position) {
            int root = position;
            while (nextUnmarked_[static_cast<std::size_t>(root)] != root) {
                root = nextUnmarked_[static_cast<std::size_t>(root)];
            }
            while (nextUnmarked_[static_cast<std::size_t>(position)] != position) {
                const int next = nextUnmarked_[static_cast<std::size_t>(position)];
                nextUnmarked_[static_cast<std::size_t>(position)] = root;
                position = next;
            }
            return root;
        }

        /**
         * @brief Marks an Euler-tour interval as invalidated by topology changes.
         *
         * @param entry Inclusive start of the dirty preorder interval.
         * @param exit Inclusive end of the dirty preorder interval.
         * @return Number of positions newly marked dirty.
         */
        [[nodiscard]] std::size_t insertDirtyInterval(int entry, int exit) {
            std::size_t inserted = 0;
            int position = findNextUnmarked(entry);
            while (position <= exit) {
                dirtyPreorder_[static_cast<std::size_t>(position)] = 1;
                ++inserted;
                nextUnmarked_[static_cast<std::size_t>(position)] = findNextUnmarked(position + 1);
                position = nextUnmarked_[static_cast<std::size_t>(position)];
            }
            return inserted;
        }

        /**
         * @brief Returns the shallowest Euler-tour occurrence in an interval.
         *
         * @param left Inclusive left Euler-tour index.
         * @param right Inclusive right Euler-tour index.
         * @return Euler-tour index of the shallowest occurrence in the interval.
         */
        [[nodiscard]] int rangeMinimum(int left, int right) const {
            const int length = right - left + 1;
            const int level = rmqLog2_[static_cast<std::size_t>(length)];
            const int span = 1 << level;
            const std::size_t row = static_cast<std::size_t>(level) * rmqStride_;
            const int resultLeft = rmqSparseTable_[row + static_cast<std::size_t>(left)];
            const int resultRight = rmqSparseTable_[row + static_cast<std::size_t>(right - span + 1)];
            return betterEulerIndex(resultLeft, resultRight);
        }

        /** @brief Rebuilds the complete Euler-tour and range-minimum snapshot. */
        void rebuild() {
            if (!enabled_ || tree_ == nullptr) {
                return;
            }
            const MorphologicalTree& topology = tree_->topology();
            const std::size_t numSlots = static_cast<std::size_t>(topology.getNumInternalNodeSlots());
            firstOccurrence_.assign(numSlots, -1);
            preorderEntry_.assign(numSlots, -1);
            preorderExit_.assign(numSlots, -1);
            euler_.clear();
            eulerDepth_.clear();
            euler_.reserve(static_cast<std::size_t>(std::max(1, 2 * topology.getNumNodes() - 1)));
            eulerDepth_.reserve(euler_.capacity());

            std::vector<TraversalFrame> stack;
            stack.reserve(static_cast<std::size_t>(std::max(1, topology.getNumNodes())));
            int preorder = 0;
            const auto pushNode = [&](NodeId nodeId, int depth) {
                const std::size_t index = static_cast<std::size_t>(nodeId);
                preorderEntry_[index] = preorder++;
                firstOccurrence_[index] = static_cast<int>(euler_.size());
                appendEuler(nodeId, depth);
                stack.push_back(TraversalFrame{nodeId, topology.getFirstChild(nodeId), depth});
            };
            pushNode(topology.getRoot(), 0);
            while (!stack.empty()) {
                TraversalFrame& frame = stack.back();
                if (frame.nextChild != InvalidNode) {
                    const NodeId child = frame.nextChild;
                    frame.nextChild = topology.getNextSibling(child);
                    pushNode(child, frame.depth + 1);
                    continue;
                }
                preorderExit_[static_cast<std::size_t>(frame.nodeId)] = preorder - 1;
                stack.pop_back();
                if (!stack.empty()) {
                    appendEuler(stack.back().nodeId, stack.back().depth);
                }
            }

            const std::size_t eulerSize = euler_.size();
            rmqLog2_.assign(eulerSize + 1, 0);
            for (std::size_t index = 2; index <= eulerSize; ++index) {
                rmqLog2_[index] = rmqLog2_[index / 2] + 1;
            }
            rmqLevels_ = eulerSize == 0 ? 0 : rmqLog2_[eulerSize] + 1;
            rmqStride_ = eulerSize;
            rmqSparseTable_.assign(static_cast<std::size_t>(rmqLevels_) * rmqStride_, -1);
            for (std::size_t index = 0; index < eulerSize; ++index) {
                rmqSparseTable_[index] = static_cast<int>(index);
            }
            for (int level = 1; level < rmqLevels_; ++level) {
                const std::size_t span = std::size_t{1} << level;
                const std::size_t half = span >> 1U;
                const std::size_t row = static_cast<std::size_t>(level) * rmqStride_;
                const std::size_t previousRow = static_cast<std::size_t>(level - 1) * rmqStride_;
                for (std::size_t index = 0; index + span <= eulerSize; ++index) {
                    rmqSparseTable_[row + index] = betterEulerIndex(rmqSparseTable_[previousRow + index], rmqSparseTable_[previousRow + index + half]);
                }
            }

            snapshotAliveNodes_ = static_cast<std::size_t>(topology.getNumNodes());
            blockSize_ = std::clamp<std::size_t>(
                (static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(std::max<std::size_t>(snapshotAliveNodes_, 1))))) + 1) / 2, 8, 256);
            dirtyPreorder_.assign(snapshotAliveNodes_, std::uint8_t{0});
            nextUnmarked_.resize(snapshotAliveNodes_ + 1);
            std::iota(nextUnmarked_.begin(), nextUnmarked_.end(), 0);
            dirtyRootMarks_.resize(numSlots);
            mutationsSinceRebuild_ = 0;
            rebuildPending_ = false;
        }

      public:
        /**
         * @brief Binds the dynamic index to the current tree topology.
         *
         * @param tree Tree processed by the operation.
         * @param policy Policy controlling the operation.
         */
        void bind(const tree_t& tree, SaturatedMinMaxLcaPolicy policy) {
            tree_ = &tree;
            enabled_ = policy == SaturatedMinMaxLcaPolicy::BlockedSnapshot;
            rebuildPending_ = enabled_;
            mutationsSinceRebuild_ = 0;
            snapshotAliveNodes_ = 0;
            euler_.clear();
            eulerDepth_.clear();
            firstOccurrence_.clear();
            preorderEntry_.clear();
            preorderExit_.clear();
            rmqLog2_.clear();
            rmqSparseTable_.clear();
            rmqStride_ = 0;
            rmqLevels_ = 0;
            dirtyPreorder_.clear();
            nextUnmarked_.clear();
            dirtyRootMarks_ = GenerationStampSet{};
            mirroredParent_.clear();
            mirroredAlive_.clear();
            if (enabled_) {
                const MorphologicalTree& topology = tree.topology();
                const std::size_t numSlots = static_cast<std::size_t>(topology.getNumInternalNodeSlots());
                mirroredParent_.assign(numSlots, InvalidNode);
                mirroredAlive_.assign(numSlots, std::uint8_t{0});
                for (NodeId nodeId : topology.getAliveNodeIds()) {
                    const std::size_t index = static_cast<std::size_t>(nodeId);
                    mirroredAlive_[index] = 1;
                    const NodeId parent = topology.getNodeParent(nodeId);
                    if (parent != nodeId) {
                        mirroredParent_[index] = parent;
                    }
                }
            }
        }

        /**
         * @brief Finds the lowest common ancestor of two current tree nodes.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @return The matching node identifier, or the operation-specific sentinel when absent.
         */
        [[nodiscard]] std::optional<NodeId> findLowestCommonAncestor(NodeId first, NodeId second) {
            if (!enabled_) {
                return std::nullopt;
            }
            if (rebuildPending_) {
                rebuild();
            }
            if (snapshotNodeIsDirty(first) || snapshotNodeIsDirty(second)) {
                return std::nullopt;
            }
            const int firstPosition = firstOccurrence_[static_cast<std::size_t>(first)];
            const int secondPosition = firstOccurrence_[static_cast<std::size_t>(second)];
            if (firstPosition < 0 || secondPosition < 0) {
                return std::nullopt;
            }
            const int left = std::min(firstPosition, secondPosition);
            const int right = std::max(firstPosition, secondPosition);
            const int eulerPosition = rangeMinimum(left, right);
            if (eulerPosition < 0) {
                return std::nullopt;
            }
            const NodeId result = euler_[static_cast<std::size_t>(eulerPosition)];
            if (!tree_->topology().isAlive(result)) {
                return std::nullopt;
            }
            return result;
        }

        /**
         * @brief Updates the dynamic index after component-tree topology changes.
         *
         * @param changedRoots Roots whose topology changed.
         */
        void noteMutation(std::span<const NodeId> changedRoots) {
            if (!enabled_ || changedRoots.empty()) {
                return;
            }
            if (rebuildPending_ || snapshotAliveNodes_ == 0) {
                rebuildPending_ = true;
            }
            const MorphologicalTree& topology = tree_->topology();
            bool topologyChanged = false;
            for (NodeId nodeId : changedRoots) {
                if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= preorderEntry_.size()) {
                    rebuildPending_ = true;
                    continue;
                }
                const std::size_t nodeIndex = static_cast<std::size_t>(nodeId);
                const bool currentAlive = topology.isAlive(nodeId);
                NodeId currentParent = InvalidNode;
                if (currentAlive) {
                    const NodeId parent = topology.getNodeParent(nodeId);
                    if (parent != nodeId) {
                        currentParent = parent;
                    }
                }
                if ((mirroredAlive_[nodeIndex] != 0) == currentAlive && mirroredParent_[nodeIndex] == currentParent) {
                    continue;
                }
                mirroredAlive_[nodeIndex] = currentAlive ? std::uint8_t{1} : std::uint8_t{0};
                mirroredParent_[nodeIndex] = currentParent;
                topologyChanged = true;
                if (rebuildPending_ || snapshotAliveNodes_ == 0) {
                    continue;
                }
                const int entry = preorderEntry_[nodeIndex];
                const int exit = preorderExit_[nodeIndex];
                if (entry < 0 || exit < entry) {
                    rebuildPending_ = true;
                    continue;
                }
                if (dirtyRootMarks_.isMarked(nodeIndex)) {
                    continue;
                }
                dirtyRootMarks_.mark(nodeIndex);
                static_cast<void>(insertDirtyInterval(entry, exit));
            }
            if (topologyChanged) {
                ++mutationsSinceRebuild_;
            }
            if (mutationsSinceRebuild_ >= blockSize_) {
                rebuildPending_ = true;
            }
        }
    };

    /**
     * @brief Exact dynamic rooted-tree LCA index based on link-cut trees.
     *
     * The component-tree root orientation is kept fixed between updates; no
     * evert operation is required. `access(u); access(v)` returns the LCA, and
     * changed parent edges are mirrored after each committed dual-tree edit.
     */
    class LinkCutDynamicLcaIndex {
      private:
        /** @brief Stores one node of the auxiliary link-cut forest. */
        struct LinkCutNode {
            /** @brief Stores the left. */
            NodeId left = InvalidNode;
            /** @brief Stores the right. */
            NodeId right = InvalidNode;
            /** @brief Stores the auxiliary parent. */
            NodeId auxiliaryParent = InvalidNode;
        };

        /** @brief References the tree used by the component. */
        const tree_t* tree_ = nullptr;
        /** @brief Indicates whether the link-cut index is enabled. */
        bool enabled_ = false;
        /** @brief Stores the nodes. */
        std::vector<LinkCutNode> nodes_;
        /** @brief Stores the represented parent. */
        std::vector<NodeId> representedParent_;
        /** @brief Stores the represented alive. */
        std::vector<std::uint8_t> representedAlive_;
        /** @brief Stores the changed scratch. */
        std::vector<NodeId> changedScratch_;

        /**
         * @brief Checks whether a link-cut node is an auxiliary-tree root.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @return `true` when the node is an auxiliary-tree root; otherwise `false`.
         */
        [[nodiscard]] bool isAuxiliaryRoot(NodeId nodeId) const noexcept {
            const NodeId parent = nodes_[static_cast<std::size_t>(nodeId)].auxiliaryParent;
            if (parent == InvalidNode) {
                return true;
            }
            const LinkCutNode& parentNode = nodes_[static_cast<std::size_t>(parent)];
            return parentNode.left != nodeId && parentNode.right != nodeId;
        }

        /**
         * @brief Rotates one node in its auxiliary link-cut tree.
         *
         * @param nodeId Identifier of the node processed by the operation.
         */
        void rotate(NodeId nodeId) {
            LinkCutNode& node = nodes_[static_cast<std::size_t>(nodeId)];
            const NodeId parentId = node.auxiliaryParent;
            LinkCutNode& parent = nodes_[static_cast<std::size_t>(parentId)];
            const NodeId grandParentId = parent.auxiliaryParent;
            const bool isLeft = parent.left == nodeId;
            const NodeId middle = isLeft ? node.right : node.left;

            if (!isAuxiliaryRoot(parentId)) {
                LinkCutNode& grandParent = nodes_[static_cast<std::size_t>(grandParentId)];
                if (grandParent.left == parentId) {
                    grandParent.left = nodeId;
                } else {
                    grandParent.right = nodeId;
                }
            }
            node.auxiliaryParent = grandParentId;
            if (isLeft) {
                node.right = parentId;
                parent.left = middle;
            } else {
                node.left = parentId;
                parent.right = middle;
            }
            parent.auxiliaryParent = nodeId;
            if (middle != InvalidNode) {
                nodes_[static_cast<std::size_t>(middle)].auxiliaryParent = parentId;
            }
        }

        /**
         * @brief Splays one node to the root of its auxiliary link-cut tree.
         *
         * @param nodeId Identifier of the node processed by the operation.
         */
        void splay(NodeId nodeId) {
            while (!isAuxiliaryRoot(nodeId)) {
                const NodeId parent = nodes_[static_cast<std::size_t>(nodeId)].auxiliaryParent;
                if (!isAuxiliaryRoot(parent)) {
                    const NodeId grandParent = nodes_[static_cast<std::size_t>(parent)].auxiliaryParent;
                    const bool nodeIsLeft = nodes_[static_cast<std::size_t>(parent)].left == nodeId;
                    const bool parentIsLeft = nodes_[static_cast<std::size_t>(grandParent)].left == parent;
                    if (nodeIsLeft == parentIsLeft) {
                        rotate(parent);
                    } else {
                        rotate(nodeId);
                    }
                }
                rotate(nodeId);
            }
        }

        /**
         * @brief Exposes the represented-tree path ending at a link-cut node.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @return Last preferred-path root visited while exposing `nodeId`.
         */
        [[nodiscard]] NodeId access(NodeId nodeId) {
            NodeId last = InvalidNode;
            for (NodeId cursor = nodeId; cursor != InvalidNode;) {
                splay(cursor);
                LinkCutNode& cursorNode = nodes_[static_cast<std::size_t>(cursor)];
                const NodeId pathParent = cursorNode.auxiliaryParent;
                cursorNode.right = last;
                if (last != InvalidNode) {
                    nodes_[static_cast<std::size_t>(last)].auxiliaryParent = cursor;
                }
                last = cursor;
                cursor = pathParent;
            }
            splay(nodeId);
            return last;
        }

        /**
         * @brief Cuts a node from its represented-tree parent in the link-cut forest.
         *
         * @param nodeId Identifier of the node processed by the operation.
         */
        void cutRepresentedParent(NodeId nodeId) {
            static_cast<void>(access(nodeId));
            LinkCutNode& node = nodes_[static_cast<std::size_t>(nodeId)];
            if (node.left != InvalidNode) {
                nodes_[static_cast<std::size_t>(node.left)].auxiliaryParent = InvalidNode;
                node.left = InvalidNode;
            }
        }

        /**
         * @brief Links a node to its represented-tree parent in the link-cut forest.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @param parentId Identifier of the parent node.
         */
        void linkRepresentedParent(NodeId nodeId, NodeId parentId) {
            static_cast<void>(access(nodeId));
            nodes_[static_cast<std::size_t>(nodeId)].auxiliaryParent = parentId;
        }

      public:
        /**
         * @brief Binds the dynamic index to the current tree topology.
         *
         * @param tree Tree processed by the operation.
         * @param policy Policy controlling the operation.
         */
        void bind(const tree_t& tree, SaturatedMinMaxLcaPolicy policy) {
            tree_ = &tree;
            enabled_ = policy == SaturatedMinMaxLcaPolicy::LinkCut;
            nodes_.clear();
            representedParent_.clear();
            representedAlive_.clear();
            changedScratch_.clear();
            if (!enabled_) {
                return;
            }
            const MorphologicalTree& topology = tree.topology();
            const std::size_t numSlots = static_cast<std::size_t>(topology.getNumInternalNodeSlots());
            nodes_.assign(numSlots, LinkCutNode{});
            representedParent_.assign(numSlots, InvalidNode);
            representedAlive_.assign(numSlots, std::uint8_t{0});
            for (NodeId nodeId : topology.getAliveNodeIds()) {
                const std::size_t index = static_cast<std::size_t>(nodeId);
                representedAlive_[index] = 1;
                const NodeId parent = topology.getNodeParent(nodeId);
                if (parent != nodeId) {
                    representedParent_[index] = parent;
                    nodes_[index].auxiliaryParent = parent;
                }
            }
            changedScratch_.reserve(64);
        }

        /**
         * @brief Finds the lowest common ancestor of two current tree nodes.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @return The matching node identifier, or the operation-specific sentinel when absent.
         */
        [[nodiscard]] std::optional<NodeId> findLowestCommonAncestor(NodeId first, NodeId second) {
            if (!enabled_) {
                return std::nullopt;
            }
            if (first == second) {
                return first;
            }
            static_cast<void>(access(first));
            const NodeId result = access(second);
            if (result == InvalidNode || !tree_->topology().isAlive(result)) {
                throw std::runtime_error("Min/max residual link-cut LCA returned an invalid node.");
            }
            return result;
        }

        /**
         * @brief Updates the dynamic index after component-tree topology changes.
         *
         * @param changedRoots Roots whose topology changed.
         */
        void noteMutation(std::span<const NodeId> changedRoots) {
            if (!enabled_ || changedRoots.empty()) {
                return;
            }
            const MorphologicalTree& topology = tree_->topology();
            changedScratch_.clear();
            for (NodeId nodeId : changedRoots) {
                if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= nodes_.size()) {
                    continue;
                }
                const std::size_t index = static_cast<std::size_t>(nodeId);
                const bool currentAlive = topology.isAlive(nodeId);
                NodeId currentParent = InvalidNode;
                if (currentAlive) {
                    const NodeId parent = topology.getNodeParent(nodeId);
                    if (parent != nodeId) {
                        currentParent = parent;
                    }
                }
                if ((representedAlive_[index] != 0) != currentAlive || representedParent_[index] != currentParent) {
                    changedScratch_.push_back(nodeId);
                }
            }

            for (NodeId nodeId : changedScratch_) {
                const std::size_t index = static_cast<std::size_t>(nodeId);
                if (representedAlive_[index] != 0 && representedParent_[index] != InvalidNode) {
                    cutRepresentedParent(nodeId);
                }
            }
            for (NodeId nodeId : changedScratch_) {
                const std::size_t index = static_cast<std::size_t>(nodeId);
                const bool currentAlive = topology.isAlive(nodeId);
                NodeId currentParent = InvalidNode;
                if (currentAlive) {
                    const NodeId parent = topology.getNodeParent(nodeId);
                    if (parent != nodeId) {
                        currentParent = parent;
                        linkRepresentedParent(nodeId, currentParent);
                    }
                }
                representedAlive_[index] = currentAlive ? std::uint8_t{1} : std::uint8_t{0};
                representedParent_[index] = currentParent;
            }
        }
    };

    /** @brief Dispatches the configured exact LCA backend. */
    class DynamicLcaIndex {
      private:
        /** @brief Stores the policy. */
        SaturatedMinMaxLcaPolicy policy_ = SaturatedMinMaxLcaPolicy::ParentClimb;
        /** @brief Stores the blocked. */
        BlockedDynamicLcaIndex blocked_;
        /** @brief Stores the link cut. */
        LinkCutDynamicLcaIndex linkCut_;

      public:
        /**
         * @brief Binds the dynamic index to the current tree topology.
         *
         * @param tree Tree processed by the operation.
         * @param policy Policy controlling the operation.
         */
        void bind(const tree_t& tree, SaturatedMinMaxLcaPolicy policy) {
            policy_ = policy;
            blocked_.bind(tree, policy);
            linkCut_.bind(tree, policy);
        }

        /**
         * @brief Finds the lowest common ancestor of two current tree nodes.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @return The matching node identifier, or the operation-specific sentinel when absent.
         */
        [[nodiscard]] std::optional<NodeId> findLowestCommonAncestor(NodeId first, NodeId second) {
            if (policy_ == SaturatedMinMaxLcaPolicy::BlockedSnapshot) {
                return blocked_.findLowestCommonAncestor(first, second);
            }
            if (policy_ == SaturatedMinMaxLcaPolicy::LinkCut) {
                return linkCut_.findLowestCommonAncestor(first, second);
            }
            return std::nullopt;
        }

        /**
         * @brief Updates the dynamic index after component-tree topology changes.
         *
         * @param changedRoots Roots whose topology changed.
         */
        void noteMutation(std::span<const NodeId> changedRoots) {
            blocked_.noteMutation(changedRoots);
            linkCut_.noteMutation(changedRoots);
        }
    };

    /** @brief Recomputes one proper part and its external incidences on demand. */
    class RecomputedProperPartBoundary {
      public:
        /** @brief Scratch storage valid until the next acquisition. */
        struct Bag {
            std::vector<NodeId> supportPixels;             ///< Current direct proper part.
            std::vector<NodeId> externalPixelsByIncidence; ///< Outside incidences.
        };

      private:
        /** @brief Stores the bag. */
        Bag bag_;

      public:
        /**
         * @brief Returns the proper-part boundary required for the current candidate.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @param tree Tree processed by the operation.
         * @param adjacency Adjacency relation used by the operation.
         * @return A reference to the requested object.
         */
        [[nodiscard]] const Bag& acquire(NodeId nodeId, const tree_t& tree, const RegularGridAdjacency2D& adjacency) {
            const MorphologicalTree& topology = tree.topology();
            bag_.supportPixels.clear();
            bag_.externalPixelsByIncidence.clear();
            for (NodeId pixel : topology.getProperParts(nodeId)) {
                bag_.supportPixels.push_back(pixel);
                for (NodeId neighbor : adjacency.getNeighborIndices(pixel)) {
                    if (topology.getProperPartOwner(neighbor) != nodeId) {
                        bag_.externalPixelsByIncidence.push_back(neighbor);
                    }
                }
            }
            return bag_;
        }
    };

    /**
     * @brief Union-only boundary index for the flat zones of the current image.
     *
     * Elementary levelings merge flat zones and never split them. This index
     * therefore remains independent of the internal rewiring of the mutable
     * max-tree and min-tree. Boundary bags keep one outside-pixel entry per
     * original crossing incidence and discard entries that become internal
     * lazily after a union.
     */
    class IncrementalFlatZoneBoundaryCache {
      public:
        /** @brief Non-owning view of a current merged flat zone. */
        struct View {
            NodeId root = InvalidNode;                         ///< Dense flat-zone root.
            std::span<const NodeId> supportPixels;             ///< Current zone support.
            std::span<const NodeId> externalPixelsByIncidence; ///< Outside incidences.
            NodeId stableSpatialKey = InvalidNode;             ///< Tie-breaking key.
            bool containsInfinity = false;                     ///< Whether the zone owns the exterior seed.
        };

      private:
        /** @brief Stores one incremental flat-zone boundary set. */
        struct Zone {
            /** @brief Stores the level. */
            altitude_t level{};
            /** @brief Stores the members. */
            std::vector<NodeId> members;
            /** @brief Stores the external pixels by incidence. */
            std::vector<NodeId> externalPixelsByIncidence;
            /** @brief Stores the stable spatial key. */
            NodeId stableSpatialKey = InvalidNode;
            /** @brief Indicates whether the flat zone contains the infinity pixel. */
            bool containsInfinity = false;
        };

        /** @brief Indicates whether incremental flat-zone caching is enabled. */
        bool enabled_ = true;
        /** @brief Stores the parent. */
        std::vector<NodeId> parent_;
        /** @brief Stores the size. */
        std::vector<std::size_t> size_;
        /** @brief Stores the zones. */
        std::vector<Zone> zones_;

        /**
         * @brief Finds the current disjoint-set representative of a pixel.
         *
         * @param pixel Row-major pixel identifier.
         * @return The matching node identifier, or the operation-specific sentinel when absent.
         */
        [[nodiscard]] NodeId findRoot(NodeId pixel) {
            NodeId root = pixel;
            while (parent_[static_cast<std::size_t>(root)] != root) {
                root = parent_[static_cast<std::size_t>(root)];
            }
            while (pixel != root) {
                const NodeId next = parent_[static_cast<std::size_t>(pixel)];
                parent_[static_cast<std::size_t>(pixel)] = root;
                pixel = next;
            }
            return root;
        }

        /**
         * @brief Unites two equal-altitude pixels during cache initialization.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @return Representative of the initialized disjoint set.
         */
        NodeId uniteInitial(NodeId first, NodeId second) {
            first = findRoot(first);
            second = findRoot(second);
            if (first == second) {
                return first;
            }
            if (size_[static_cast<std::size_t>(first)] < size_[static_cast<std::size_t>(second)]) {
                std::swap(first, second);
            }
            parent_[static_cast<std::size_t>(second)] = first;
            size_[static_cast<std::size_t>(first)] += size_[static_cast<std::size_t>(second)];
            return first;
        }

        /**
         * @brief Moves boundary entries from the smaller set into the larger set.
         *
         * @param target Boundary list that receives all entries.
         * @param source Boundary list consumed and cleared by the merge.
         */
        static void mergeSmallToLarge(std::vector<NodeId>& target, std::vector<NodeId>& source) {
            if (target.size() < source.size()) {
                target.swap(source);
            }
            target.insert(target.end(), source.begin(), source.end());
            source.clear();
        }

        /**
         * @brief Unites two current flat zones and preserves their boundary cache.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @param level Altitude level processed by the operation.
         * @return Representative of the merged current flat zone.
         */
        NodeId uniteCurrent(NodeId first, NodeId second, altitude_t level) {
            first = findRoot(first);
            second = findRoot(second);
            if (first == second) {
                zones_[static_cast<std::size_t>(first)].level = level;
                return first;
            }
            if (zones_[static_cast<std::size_t>(first)].members.size() < zones_[static_cast<std::size_t>(second)].members.size()) {
                std::swap(first, second);
            }
            Zone& target = zones_[static_cast<std::size_t>(first)];
            Zone& source = zones_[static_cast<std::size_t>(second)];
            mergeSmallToLarge(target.members, source.members);
            mergeSmallToLarge(target.externalPixelsByIncidence, source.externalPixelsByIncidence);
            target.stableSpatialKey = std::min(target.stableSpatialKey, source.stableSpatialKey);
            target.containsInfinity = target.containsInfinity || source.containsInfinity;
            target.level = level;
            parent_[static_cast<std::size_t>(second)] = first;
            size_[static_cast<std::size_t>(first)] += size_[static_cast<std::size_t>(second)];
            return first;
        }

      public:
        /**
         * @brief Constructs a `IncrementalFlatZoneBoundaryCache` instance.
         *
         * @param enabled Whether the optional strategy is enabled.
         */
        explicit IncrementalFlatZoneBoundaryCache(bool enabled = true) : enabled_(enabled) {}

        /**
         * @brief Returns the current flat-zone representative of every pixel.
         *
         * @return Region representative indexed by row-major pixel.
         */
        [[nodiscard]] std::vector<detail::RegionId> currentRegionByPixel() {
            std::vector<detail::RegionId> regionByPixel(parent_.size(), detail::InvalidRegion);
            for (NodeId pixel = 0; pixel < static_cast<NodeId>(parent_.size()); ++pixel) {
                regionByPixel[static_cast<std::size_t>(pixel)] = static_cast<detail::RegionId>(findRoot(pixel));
            }
            return regionByPixel;
        }

        /**
         * @brief Returns the current flat-zone root of a pixel.
         *
         * @param pixel Row-major pixel identifier.
         * @return Current disjoint-set representative of `pixel`.
         */
        [[nodiscard]] NodeId currentRoot(NodeId pixel) {
            if (!enabled_ || pixel < 0 || pixel >= static_cast<NodeId>(parent_.size())) {
                throw std::out_of_range("Flat-zone root query lies outside the active domain.");
            }
            return findRoot(pixel);
        }

        /**
         * @brief Initializes the data structure from the current component trees.
         *
         * @param image Image processed by the operation.
         * @param adjacency Adjacency relation used by the operation.
         * @param infinityPixel Row-major pixel used as the exterior seed.
         */
        void initialize(const image_ptr_t& image, const RegularGridAdjacency2D& adjacency, NodeId infinityPixel) {
            parent_.clear();
            size_.clear();
            zones_.clear();
            if (!enabled_) {
                return;
            }
            const std::size_t numPixels = static_cast<std::size_t>(image->getSize());
            parent_.resize(numPixels);
            size_.assign(numPixels, 1);
            zones_.resize(numPixels);
            std::iota(parent_.begin(), parent_.end(), NodeId{0});

            for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
                for (NodeId neighbor : adjacency.getForwardNeighborIndices(pixel)) {
                    if ((*image)[pixel] == (*image)[neighbor]) {
                        uniteInitial(pixel, neighbor);
                    }
                }
            }
            for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
                const NodeId root = findRoot(pixel);
                Zone& zone = zones_[static_cast<std::size_t>(root)];
                zone.level = (*image)[pixel];
                zone.members.push_back(pixel);
                if (zone.stableSpatialKey == InvalidNode) {
                    zone.stableSpatialKey = pixel;
                } else {
                    zone.stableSpatialKey = std::min(zone.stableSpatialKey, pixel);
                }
                zone.containsInfinity = zone.containsInfinity || pixel == infinityPixel;
            }
            for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
                const NodeId pixelRoot = findRoot(pixel);
                for (NodeId neighbor : adjacency.getForwardNeighborIndices(pixel)) {
                    const NodeId neighborRoot = findRoot(neighbor);
                    if (pixelRoot == neighborRoot) {
                        continue;
                    }
                    zones_[static_cast<std::size_t>(pixelRoot)].externalPixelsByIncidence.push_back(neighbor);
                    zones_[static_cast<std::size_t>(neighborRoot)].externalPixelsByIncidence.push_back(pixel);
                }
            }
        }

        /**
         * @brief Returns the proper-part boundary required for the current candidate.
         *
         * @param candidateNode Component-tree node whose boundary is requested.
         * @param primal Component tree that owns `candidateNode`.
         * @return View of the cached flat-zone members and external boundary incidences.
         */
        [[nodiscard]] View acquire(NodeId candidateNode, const tree_t& primal) {
            if (!enabled_) {
                throw std::logic_error("The incremental flat-zone boundary cache is disabled.");
            }
            const MorphologicalTree& topology = primal.topology();
            const auto properParts = topology.getProperParts(candidateNode);
            const auto first = properParts.begin();
            if (first == properParts.end()) {
                throw std::runtime_error("A min/max residual candidate has an empty proper part.");
            }
            const NodeId root = findRoot(*first);
            Zone& zone = zones_[static_cast<std::size_t>(root)];
            if (zone.members.size() != static_cast<std::size_t>(topology.getNumProperParts(candidateNode))) {
                throw std::runtime_error("The incremental flat-zone boundary cache diverged from the component-tree proper part.");
            }
#ifndef NDEBUG
            for (NodeId pixel : zone.members) {
                if (topology.getProperPartOwner(pixel) != candidateNode) {
                    throw std::runtime_error("The incremental flat-zone boundary cache has a foreign proper part.");
                }
            }
#endif
            auto output = zone.externalPixelsByIncidence.begin();
            for (NodeId neighbor : zone.externalPixelsByIncidence) {
                if (findRoot(neighbor) == root) {
                    continue;
                }
                *output++ = neighbor;
            }
            zone.externalPixelsByIncidence.erase(output, zone.externalPixelsByIncidence.end());
            return View{root, zone.members, zone.externalPixelsByIncidence, zone.stableSpatialKey, zone.containsInfinity};
        }

        /**
         * @brief Computes the cached descriptor of a current component-tree node.
         *
         * @param tree Tree processed by the operation.
         * @param nodeId Identifier of the node processed by the operation.
         * @return Descriptor containing area, stable spatial key, and infinity membership.
         */
        [[nodiscard]] CandidateDescriptor describeNode(const tree_t& tree, NodeId nodeId) {
            const MorphologicalTree& topology = tree.topology();
            const auto properParts = topology.getProperParts(nodeId);
            const auto first = properParts.begin();
            if (first == properParts.end()) {
                return CandidateDescriptor{};
            }
            const NodeId root = findRoot(*first);
            const Zone& zone = zones_[static_cast<std::size_t>(root)];
            if (zone.members.size() != static_cast<std::size_t>(topology.getNumProperParts(nodeId))) {
                throw std::runtime_error("Flat-zone metadata diverged from a component-tree proper part.");
            }
            return CandidateDescriptor{static_cast<int>(zone.members.size()), zone.stableSpatialKey, zone.containsInfinity};
        }

        /**
         * @brief Collects the current flat-zone roots that must be merged.
         *
         * @param selectedRoot Current flat-zone root selected for contraction.
         * @param targetLevel Target altitude level of the merge.
         * @param boundaryPixelsByIncidence External boundary pixels, repeated by incidence.
         * @param marks Generation-stamped set used to deduplicate collected roots.
         * @param mergeRoots Roots collected or updated by the merge.
         */
        void collectMergeRoots(NodeId selectedRoot, altitude_t targetLevel, std::span<const NodeId> boundaryPixelsByIncidence, GenerationStampSet& marks,
                               std::vector<NodeId>& mergeRoots) {
            marks.resetAll();
            mergeRoots.clear();
            selectedRoot = findRoot(selectedRoot);
            for (NodeId pixel : boundaryPixelsByIncidence) {
                const NodeId root = findRoot(pixel);
                if (root == selectedRoot || zones_[static_cast<std::size_t>(root)].level != targetLevel || marks.isMarked(static_cast<std::size_t>(root))) {
                    continue;
                }
                marks.mark(static_cast<std::size_t>(root));
                mergeRoots.push_back(root);
            }
            if (mergeRoots.empty()) {
                throw std::runtime_error("An elementary leveling found no adjacent flat zone at its merging level.");
            }
        }

        /**
         * @brief Merges current flat zones at a common altitude level.
         *
         * @param selectedRoot Current flat-zone root selected for contraction.
         * @param mergeRoots Roots collected or updated by the merge.
         * @param targetLevel Target altitude level of the merge.
         * @return Representative of the flat zone after all requested merges.
         */
        [[nodiscard]] NodeId mergeAtLevel(NodeId selectedRoot, std::span<const NodeId> mergeRoots, altitude_t targetLevel) {
            NodeId root = findRoot(selectedRoot);
            for (NodeId neighborRoot : mergeRoots) {
                root = uniteCurrent(root, neighborRoot, targetLevel);
            }
            zones_[static_cast<std::size_t>(root)].level = targetLevel;
            return root;
        }
    };

    /** @brief Owns reusable buffers for exact saturation certification. */
    struct CertificationScratch {
        /** @brief Stores the support marks. */
        GenerationStampSet supportMarks;
        /** @brief Stores the visited marks. */
        GenerationStampSet visitedMarks;
        /** @brief Stores the boundary pixel marks. */
        GenerationStampSet boundaryPixelMarks;
        /** @brief Stores the boundary owner marks. */
        GenerationStampSet boundaryOwnerMarks;
        /** @brief Stores the support owner marks. */
        GenerationStampSet supportOwnerMarks;
        /** @brief Stores the boundary owners. */
        std::vector<NodeId> boundaryOwners;
        /** @brief Stores the boundary pixels. */
        std::vector<NodeId> boundaryPixels;
        /** @brief Stores the support pixels. */
        std::span<const NodeId> supportPixels;
        /** @brief Stores the support owners. */
        std::vector<NodeId> supportOwners;
        /** @brief Stores the flat zone root marks. */
        GenerationStampSet flatZoneRootMarks;
        /** @brief Stores the flat zone merge roots. */
        std::vector<NodeId> flatZoneMergeRoots;
        /** @brief Stores the flat zone root. */
        NodeId flatZoneRoot = InvalidNode;
        /** @brief Stores the frontier. */
        std::vector<NodeId> frontier;
        /** @brief Stores the complement label by pixel. */
        std::vector<int> complementLabelByPixel;
        /** @brief Stores the complement label parent. */
        std::vector<int> complementLabelParent;
        /** @brief Stores the complement label rank. */
        std::vector<std::uint8_t> complementLabelRank;
        /** @brief Stores the dual extremal owner. */
        NodeId dualExtremalOwner = InvalidNode;
        /** @brief Stores the whole support owner. */
        NodeId wholeSupportOwner = InvalidNode;

        /**
         * @brief Constructs a `CertificationScratch` instance.
         *
         * @param numPixels Num pixels used by the operation.
         * @param maxNodeSlots Max node slots used by the operation.
         */
        CertificationScratch(std::size_t numPixels, std::size_t maxNodeSlots)
            : supportMarks(numPixels), visitedMarks(numPixels), boundaryPixelMarks(numPixels), boundaryOwnerMarks(maxNodeSlots),
              supportOwnerMarks(maxNodeSlots), flatZoneRootMarks(numPixels), complementLabelByPixel(numPixels, -1) {
            frontier.reserve(numPixels);
            boundaryOwners.reserve(32);
            boundaryPixels.reserve(64);
            supportOwners.reserve(8);
            flatZoneMergeRoots.reserve(8);
        }
    };

    /** @brief Owns mutable state for one synchronized min/max residual-tree build. */
    struct ConstructionState {
        /** @brief Stores the max tree. */
        std::unique_ptr<tree_t> maxTree;
        /** @brief Stores the min tree. */
        std::unique_ptr<tree_t> minTree;
        /** @brief Stores the adjustment. */
        std::unique_ptr<Adjustment> adjustment;
        /** @brief Stores the max lowest-common-ancestor. */
        DynamicLcaIndex maxLca;
        /** @brief Stores the min lowest-common-ancestor. */
        DynamicLcaIndex minLca;
        /** @brief Stores the agenda. */
        SaturatedLeafAgenda agenda;
        /** @brief Stores the chains. */
        PersistentChains chains;
        /** @brief Stores the certification. */
        CertificationScratch certification;
        /** @brief Stores the max boundary. */
        RecomputedProperPartBoundary maxBoundary;
        /** @brief Stores the min boundary. */
        RecomputedProperPartBoundary minBoundary;
        /** @brief Stores the flat zones. */
        IncrementalFlatZoneBoundaryCache flatZones;
        /** @brief Stores the residual assembler. */
        std::optional<detail::UnionFindResidualTreeAssembler<altitude_t>> residualAssembler;
        /** @brief Stores the absorbed flat zones. */
        std::vector<detail::RegionId> absorbedFlatZones;
        /** @brief Indicates whether incremental boundary maintenance is enabled. */
        bool incrementalBoundary = true;
        /** @brief Stores the cached new head by prev head. */
        std::vector<int> cachedNewHeadByPrevHead;
        /** @brief Stores the cached new head stamp by prev head. */
        std::vector<std::uint32_t> cachedNewHeadStampByPrevHead;
        /** @brief Stores the current event cache stamp. */
        std::uint32_t currentEventCacheStamp = 1;
        /** @brief Indicates whether candidates require saturation certification. */
        bool requiresSaturationCertification = true;

        /**
         * @brief Constructs a `ConstructionState` instance.
         *
         * @param tiePolicy Policy used to resolve deterministic ties.
         * @param infinityPixel Row-major pixel used as the exterior seed.
         * @param numPixels Num pixels used by the operation.
         * @param maxNodeSlots Max node slots used by the operation.
         * @param incrementalBoundary Incremental boundary used by the operation.
         * @param requiresSaturationCertification Requires saturation certification used by the operation.
         */
        ConstructionState(SdrtTiePolicy tiePolicy, NodeId infinityPixel, std::size_t numPixels, std::size_t maxNodeSlots, bool incrementalBoundary,
                          bool requiresSaturationCertification)
            : agenda(tiePolicy, infinityPixel, requiresSaturationCertification), certification(numPixels, maxNodeSlots), flatZones(incrementalBoundary),
              incrementalBoundary(incrementalBoundary), requiresSaturationCertification(requiresSaturationCertification) {}
    };

    /** @brief Stores the adjacency. */
    RegularGridAdjacency2D adjacency_;
    /** @brief Stores the infinity pixel. */
    NodeId infinityPixel_ = InvalidNode;
    /** @brief Stores the tie policy. */
    SdrtTiePolicy tiePolicy_ = SdrtTiePolicy::ContrastInvariantSpatial;
    /** @brief Stores the lowest-common-ancestor policy. */
    SaturatedMinMaxLcaPolicy lcaPolicy_ = SaturatedMinMaxLcaPolicy::ParentClimb;
    /** @brief Stores the fallback policy. */
    SaturatedMinMaxFallbackPolicy fallbackPolicy_ = SaturatedMinMaxFallbackPolicy::BoundaryMultiSource;
    /** @brief Stores the boundary policy. */
    SaturatedMinMaxBoundaryPolicy boundaryPolicy_ = SaturatedMinMaxBoundaryPolicy::IncrementalSmallToLarge;
    /** @brief Stores the eligibility policy. */
    MinMaxResidualEligibilityPolicy eligibilityPolicy_ = MinMaxResidualEligibilityPolicy::SaturatedOnly;
    /** @brief Indicates whether a successfully built result is available. */
    bool built_ = false;
    /** @brief Stores the rows. */
    int rows_ = 0;
    /** @brief Stores the cols. */
    int cols_ = 0;
    /** @brief Stores the root. */
    NodeId root_ = InvalidNode;
    /** @brief Stores the node parent. */
    std::vector<NodeId> nodeParent_;
    /** @brief Stores the proper part owner. */
    std::vector<NodeId> properPartOwner_;
    /** @brief Stores the altitude. */
    std::vector<altitude_t> altitude_;
    /** @brief Stores the topology proof. */
    mmcfilters::detail::NativeTopologyProof topologyProof_;
    /** @brief Stores the statistics. */
    MinMaxResidualTreeBuildStatistics statistics_;

    /** @brief Clears all state associated with the previous build result. */
    void clearResult() {
        built_ = false;
        rows_ = 0;
        cols_ = 0;
        root_ = InvalidNode;
        nodeParent_.clear();
        properPartOwner_.clear();
        altitude_.clear();
        topologyProof_ = mmcfilters::detail::NativeTopologyProof{};
        statistics_ = MinMaxResidualTreeBuildStatistics{};
    }

    /**
     * @brief Validates the inputs required by the construction.
     *
     * @param image Image processed by the operation.
     */
    void requireInputs(const image_ptr_t& image) const {
        if (!image || image->getNumRows() <= 0 || image->getNumCols() <= 0 || image->getSize() <= 0) {
            throw std::invalid_argument("MinMaxResidualTreeBuilder requires a non-null, non-empty image.");
        }
        if (adjacency_.getNumRows() != image->getNumRows() || adjacency_.getNumCols() != image->getNumCols()) {
            throw std::invalid_argument("MinMaxResidualTreeBuilder adjacency domain differs from the image.");
        }
        if (infinityPixel_ < 0 || infinityPixel_ >= image->getSize()) {
            throw std::invalid_argument("MinMaxResidualTreeBuilder infinity pixel lies outside the image.");
        }

        std::vector<std::uint8_t> visited(static_cast<std::size_t>(image->getSize()), std::uint8_t{0});
        std::vector<NodeId> frontier;
        frontier.reserve(static_cast<std::size_t>(image->getSize()));
        frontier.push_back(infinityPixel_);
        visited[static_cast<std::size_t>(infinityPixel_)] = 1;
        std::size_t visitedPixels = 0;
        while (!frontier.empty()) {
            const NodeId pixel = frontier.back();
            frontier.pop_back();
            ++visitedPixels;
            for (NodeId neighbor : adjacency_.getNeighborIndices(pixel)) {
                auto& mark = visited[static_cast<std::size_t>(neighbor)];
                if (mark == 0) {
                    mark = 1;
                    frontier.push_back(neighbor);
                }
            }
        }
        if (visitedPixels != static_cast<std::size_t>(image->getSize())) {
            throw std::invalid_argument("MinMaxResidualTreeBuilder requires a connected adjacency domain.");
        }
    }

    /**
     * @brief Validates one caller-supplied component-tree seed.
     *
     * @param tree Tree processed by the operation.
     * @param image Image processed by the operation.
     * @param kind Kind used by the operation.
     * @param label Label used by the operation.
     */
    static void requireSeedTree(const tree_t& tree, const image_ptr_t& image, MorphologicalTreeKind kind, const char* label) {
        const MorphologicalTree& topology = tree.topology();
        if (topology.getDescriptiveKind() != kind) {
            throw std::invalid_argument(std::string("Min/max residual ") + label + " seed has an unexpected tree type.");
        }
        if (topology.getNumRowsOfGridDomain2D() != image->getNumRows() || topology.getNumColsOfGridDomain2D() != image->getNumCols() ||
            topology.getNumTotalProperParts() != image->getSize()) {
            throw std::invalid_argument(std::string("Min/max residual ") + label + " seed domain differs from the image.");
        }
        tree.validateAltitudeBufferShape();
    }

    /**
     * @brief Initializes mutable state for one synchronized construction.
     *
     * @param image Image processed by the operation.
     * @param minTree Min-tree consumed by the operation.
     * @param maxTree Max-tree consumed by the operation.
     * @return Fully initialized state owning the synchronized component-tree seeds.
     */
    [[nodiscard]] ConstructionState initializeConstruction(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree) const {
        requireSeedTree(maxTree, image, MorphologicalTreeKind::MAX_TREE, "max-tree");
        requireSeedTree(minTree, image, MorphologicalTreeKind::MIN_TREE, "min-tree");
        const std::size_t maxNodeSlots =
            static_cast<std::size_t>(std::max(maxTree.topology().getNumInternalNodeSlots(), minTree.topology().getNumInternalNodeSlots()));
        ConstructionState state(tiePolicy_, infinityPixel_, static_cast<std::size_t>(image->getSize()), maxNodeSlots,
                                boundaryPolicy_ == SaturatedMinMaxBoundaryPolicy::IncrementalSmallToLarge,
                                eligibilityPolicy_ == MinMaxResidualEligibilityPolicy::SaturatedOnly);
        state.maxTree = std::make_unique<tree_t>(std::move(maxTree));
        state.minTree = std::make_unique<tree_t>(std::move(minTree));
        state.adjustment = std::make_unique<Adjustment>(state.minTree.get(), state.maxTree.get(), adjacency_);
        if (state.requiresSaturationCertification) {
            state.minLca.bind(*state.minTree, lcaPolicy_);
            state.maxLca.bind(*state.maxTree, lcaPolicy_);
        }
        if (state.incrementalBoundary) {
            state.flatZones.initialize(image, adjacency_, infinityPixel_);
            state.residualAssembler.emplace(static_cast<std::size_t>(image->getSize()), state.flatZones.currentRegionByPixel());
        }
        if (!state.incrementalBoundary) {
            state.chains.stackHeadByPixel.assign(static_cast<std::size_t>(image->getSize()), -1);
            const std::size_t cacheSize = std::max<std::size_t>(static_cast<std::size_t>(image->getSize()) + 1, 2);
            state.cachedNewHeadByPrevHead.assign(cacheSize, -1);
            state.cachedNewHeadStampByPrevHead.assign(cacheSize, 0);
        }
        state.agenda.initialize(*state.maxTree, *state.minTree);
        return state;
    }

    /**
     * @brief Checks whether a node is alive and remains a regional extremum.
     *
     * @param tree Tree processed by the operation.
     * @param nodeId Identifier of the node processed by the operation.
     * @return `true` when the node is a current regional extremum; otherwise `false`.
     */
    [[nodiscard]] static bool isCurrentCandidateNode(const tree_t& tree, NodeId nodeId) {
        const MorphologicalTree& topology = tree.topology();
        return nodeId >= 0 && topology.isNode(nodeId) && topology.isAlive(nodeId) && nodeId != topology.getRoot() && topology.isLeaf(nodeId) &&
               topology.getNodeParent(nodeId) != InvalidNode && topology.getNodeParent(nodeId) != nodeId && topology.getNumProperParts(nodeId) > 0;
    }

    /**
     * @brief Returns the stable spatial ordering key of a current node.
     *
     * @param topology Component-tree topology that owns `nodeId`.
     * @param nodeId Identifier of the node processed by the operation.
     * @return Smallest row-major proper-part identifier owned by the node.
     */
    [[nodiscard]] static NodeId currentSpatialKey(const MorphologicalTree& topology, NodeId nodeId) {
        NodeId key = std::numeric_limits<NodeId>::max();
        for (NodeId pixel : topology.getProperParts(nodeId)) {
            key = std::min(key, pixel);
        }
        return key == std::numeric_limits<NodeId>::max() ? InvalidNode : key;
    }

    /**
     * @brief Checks whether a cached candidate still describes the current node.
     *
     * @param state Mutable construction state updated by the operation.
     * @param tree Tree processed by the operation.
     * @param candidate Candidate processed by the operation.
     * @return `true` when the cached candidate is current; otherwise `false`.
     */
    [[nodiscard]] static bool isCandidateCurrent(ConstructionState& state, const tree_t& tree, const Candidate& candidate) {
        if (!isCurrentCandidateNode(tree, candidate.nodeId)) {
            return false;
        }
        const MorphologicalTree& topology = tree.topology();
        if (state.incrementalBoundary) {
            const CandidateDescriptor descriptor = state.flatZones.describeNode(tree, candidate.nodeId);
            return descriptor.area == candidate.area && descriptor.stableSpatialKey == candidate.stableSpatialKey;
        }
        return topology.getNumProperParts(candidate.nodeId) == candidate.area && currentSpatialKey(topology, candidate.nodeId) == candidate.stableSpatialKey;
    }

    /**
     * @brief Refreshes one candidate after synchronized tree mutations.
     *
     * @param state Mutable construction state updated by the operation.
     * @param tree Tree processed by the operation.
     * @param polarity Component-tree polarity to process.
     * @param nodeId Identifier of the node processed by the operation.
     */
    static void updateAgendaCandidate(ConstructionState& state, const tree_t& tree, Polarity polarity, NodeId nodeId) {
        if (state.incrementalBoundary && isCurrentCandidateNode(tree, nodeId)) {
            state.agenda.updateKnown(tree, polarity, nodeId, state.flatZones.describeNode(tree, nodeId));
            return;
        }
        state.agenda.update(tree, polarity, nodeId);
    }

    /**
     * @brief Ensures that persistent-chain caches cover all current node slots.
     *
     * @param state Mutable construction state updated by the operation.
     * @param requiredSize Minimum storage size required by the operation.
     */
    static void ensureChainCacheCapacity(ConstructionState& state, std::size_t requiredSize) {
        if (state.cachedNewHeadByPrevHead.size() >= requiredSize) {
            return;
        }
        std::size_t size = state.cachedNewHeadByPrevHead.empty() ? 2 : state.cachedNewHeadByPrevHead.size();
        while (size < requiredSize) {
            size *= 2;
        }
        state.cachedNewHeadByPrevHead.resize(size, -1);
        state.cachedNewHeadStampByPrevHead.resize(size, 0);
    }

    /**
     * @brief Prepends one residual event to a persistent pixel chain.
     *
     * @param state Mutable construction state updated by the operation.
     * @param eventId Identifier of the residual event.
     * @param previousHead Previous persistent-chain head.
     * @return Head index of the newly cached persistent chain.
     */
    [[nodiscard]] static int makeNewHead(ConstructionState& state, int eventId, int previousHead) {
        const std::size_t cacheIndex = static_cast<std::size_t>(previousHead + 1);
        ensureChainCacheCapacity(state, cacheIndex + 1);
        if (state.cachedNewHeadStampByPrevHead[cacheIndex] == state.currentEventCacheStamp) {
            return state.cachedNewHeadByPrevHead[cacheIndex];
        }
        const int head = static_cast<int>(state.chains.nodes.size());
        state.chains.nodes.push_back(PersistentStackNode{eventId, previousHead});
        state.cachedNewHeadByPrevHead[cacheIndex] = head;
        state.cachedNewHeadStampByPrevHead[cacheIndex] = state.currentEventCacheStamp;
        return head;
    }

    /**
     * @brief Records one residual event for every pixel in its support.
     *
     * @param state Mutable construction state updated by the operation.
     * @param altitude Altitude associated with the operation.
     * @param support Pixel or node support processed by the operation.
     */
    static void recordEvent(ConstructionState& state, altitude_t altitude, std::span<const NodeId> support) {
        if (support.empty()) {
            throw std::runtime_error("Min/max residual builder cannot record an empty support.");
        }
        ensureChainCacheCapacity(state, state.chains.nodes.size() + 1);
        ++state.currentEventCacheStamp;
        if (state.currentEventCacheStamp == 0) {
            std::fill(state.cachedNewHeadStampByPrevHead.begin(), state.cachedNewHeadStampByPrevHead.end(), 0);
            state.currentEventCacheStamp = 1;
        }
        const int eventId = static_cast<int>(state.chains.eventValuations.size());
        for (NodeId pixel : support) {
            const std::size_t index = static_cast<std::size_t>(pixel);
            const int previous = state.chains.stackHeadByPixel[index];
            state.chains.stackHeadByPixel[index] = makeNewHead(state, eventId, previous);
        }
        state.chains.eventValuations.push_back(altitude);
    }

    /**
     * @brief Computes a lowest common ancestor by altitude-aligned parent climbs.
     *
     * @param tree Tree processed by the operation.
     * @param first First value processed by the operation.
     * @param second Second value processed by the operation.
     * @return The matching node identifier, or the operation-specific sentinel when absent.
     */
    [[nodiscard]] static NodeId parentClimbLowestCommonAncestor(const tree_t& tree, NodeId first, NodeId second) {
        const MorphologicalTree& topology = tree.topology();
        const bool minTree = topology.getDescriptiveKind() == MorphologicalTreeKind::MIN_TREE;
        NodeId lhs = first;
        NodeId rhs = second;
        while (lhs != rhs) {
            const altitude_t lhsAltitude = tree.getAltitude(lhs);
            const altitude_t rhsAltitude = tree.getAltitude(rhs);
            bool climbLhs = false;
            bool climbRhs = false;
            if (lhsAltitude == rhsAltitude) {
                climbLhs = true;
                climbRhs = true;
            } else if ((minTree && lhsAltitude < rhsAltitude) || (!minTree && lhsAltitude > rhsAltitude)) {
                climbLhs = true;
            } else {
                climbRhs = true;
            }
            if (climbLhs) {
                const NodeId parent = topology.getNodeParent(lhs);
                if (parent == lhs) {
                    throw std::runtime_error("Min/max residual altitude-aligned LCA reached one root before convergence.");
                }
                lhs = parent;
            }
            if (climbRhs) {
                const NodeId parent = topology.getNodeParent(rhs);
                if (parent == rhs) {
                    throw std::runtime_error("Min/max residual altitude-aligned LCA reached one root before convergence.");
                }
                rhs = parent;
            }
            if (!topology.isAlive(lhs) || !topology.isAlive(rhs)) {
                throw std::runtime_error("Min/max residual altitude-aligned LCA reached a dead node.");
            }
        }
        return lhs;
    }

    /**
     * @brief Computes the lowest common ancestor using the configured index.
     *
     * @param tree Tree processed by the operation.
     * @param index Mutable index updated by the operation.
     * @param first First value processed by the operation.
     * @param second Second value processed by the operation.
     * @return The matching node identifier, or the operation-specific sentinel when absent.
     */
    [[nodiscard]] static NodeId lowestCommonAncestor(const tree_t& tree, DynamicLcaIndex& index, NodeId first, NodeId second) {
        const auto indexed = index.findLowestCommonAncestor(first, second);
        if (indexed.has_value()) {
            return *indexed;
        }
        return parentClimbLowestCommonAncestor(tree, first, second);
    }

    /**
     * @brief Checks complement connectivity from one reference seed.
     *
     * @param scratch Reusable scratch storage for the operation.
     * @return `true` when every complement boundary pixel is reachable from the exterior seed.
     */
    [[nodiscard]] bool exactSingleSourceComplementTraversal(CertificationScratch& scratch) const {
        scratch.visitedMarks.resetAll();
        scratch.frontier.clear();
        if (scratch.supportMarks.isMarked(static_cast<std::size_t>(infinityPixel_))) {
            return false;
        }
        if (scratch.boundaryPixels.empty()) {
            throw std::runtime_error("Min/max residual complement traversal requires a non-empty external boundary.");
        }
        const NodeId seed = scratch.boundaryPixels.front();
        scratch.frontier.push_back(seed);
        scratch.visitedMarks.mark(static_cast<std::size_t>(seed));
        std::size_t reachedBoundaryPixels = 1;
        while (!scratch.frontier.empty() && reachedBoundaryPixels < scratch.boundaryPixels.size()) {
            const NodeId pixel = scratch.frontier.back();
            scratch.frontier.pop_back();
            for (NodeId neighbor : adjacency_.getNeighborIndices(pixel)) {
                const std::size_t index = static_cast<std::size_t>(neighbor);
                if (scratch.supportMarks.isMarked(index) || scratch.visitedMarks.isMarked(index)) {
                    continue;
                }
                scratch.visitedMarks.mark(index);
                scratch.frontier.push_back(neighbor);
                if (scratch.boundaryPixelMarks.isMarked(index)) {
                    ++reachedBoundaryPixels;
                }
            }
        }
        return reachedBoundaryPixels == scratch.boundaryPixels.size();
    }

    /**
     * @brief Checks complement connectivity from all boundary seeds.
     *
     * @param scratch Reusable scratch storage for the operation.
     * @return `true` when all complement boundary seeds belong to one connected component.
     */
    [[nodiscard]] bool exactMultiSourceComplementTraversal(CertificationScratch& scratch) const {
        scratch.visitedMarks.resetAll();
        scratch.frontier.clear();
        if (scratch.supportMarks.isMarked(static_cast<std::size_t>(infinityPixel_))) {
            return false;
        }
        if (scratch.boundaryPixels.empty()) {
            throw std::runtime_error("Min/max residual complement traversal requires a non-empty external boundary.");
        }

        const std::size_t numLabels = scratch.boundaryPixels.size();
        scratch.complementLabelParent.resize(numLabels);
        std::iota(scratch.complementLabelParent.begin(), scratch.complementLabelParent.end(), 0);
        scratch.complementLabelRank.assign(numLabels, std::uint8_t{0});

        for (std::size_t label = 0; label < numLabels; ++label) {
            const NodeId pixel = scratch.boundaryPixels[label];
            const std::size_t pixelIndex = static_cast<std::size_t>(pixel);
            scratch.visitedMarks.mark(pixelIndex);
            scratch.complementLabelByPixel[pixelIndex] = static_cast<int>(label);
            scratch.frontier.push_back(pixel);
        }

        std::size_t activeLabels = numLabels;
        const auto findLabel = [&scratch](int label) {
            int root = label;
            while (scratch.complementLabelParent[static_cast<std::size_t>(root)] != root) {
                root = scratch.complementLabelParent[static_cast<std::size_t>(root)];
            }
            while (label != root) {
                const int parent = scratch.complementLabelParent[static_cast<std::size_t>(label)];
                scratch.complementLabelParent[static_cast<std::size_t>(label)] = root;
                label = parent;
            }
            return root;
        };
        const auto uniteLabels = [&scratch, &findLabel, &activeLabels](int first, int second) {
            int lhs = findLabel(first);
            int rhs = findLabel(second);
            if (lhs == rhs) {
                return;
            }
            if (scratch.complementLabelRank[static_cast<std::size_t>(lhs)] < scratch.complementLabelRank[static_cast<std::size_t>(rhs)]) {
                std::swap(lhs, rhs);
            }
            scratch.complementLabelParent[static_cast<std::size_t>(rhs)] = lhs;
            if (scratch.complementLabelRank[static_cast<std::size_t>(lhs)] == scratch.complementLabelRank[static_cast<std::size_t>(rhs)]) {
                ++scratch.complementLabelRank[static_cast<std::size_t>(lhs)];
            }
            --activeLabels;
        };

        std::size_t head = 0;
        while (head < scratch.frontier.size() && activeLabels > 1) {
            const NodeId pixel = scratch.frontier[head++];
            const int pixelLabel = scratch.complementLabelByPixel[static_cast<std::size_t>(pixel)];
            for (NodeId neighbor : adjacency_.getNeighborIndices(pixel)) {
                const std::size_t neighborIndex = static_cast<std::size_t>(neighbor);
                if (scratch.supportMarks.isMarked(neighborIndex)) {
                    continue;
                }
                if (!scratch.visitedMarks.isMarked(neighborIndex)) {
                    scratch.visitedMarks.mark(neighborIndex);
                    scratch.complementLabelByPixel[neighborIndex] = pixelLabel;
                    scratch.frontier.push_back(neighbor);
                    continue;
                }
                uniteLabels(pixelLabel, scratch.complementLabelByPixel[neighborIndex]);
                if (activeLabels == 1) {
                    break;
                }
            }
        }
        return activeLabels == 1;
    }

    /**
     * @brief Runs the configured exact complement-connectivity traversal.
     *
     * @param scratch Reusable scratch storage for the operation.
     * @return `true` when the candidate support has connected complement.
     */
    [[nodiscard]] bool exactComplementTraversal(CertificationScratch& scratch) const {
        if (fallbackPolicy_ == SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst) {
            return exactSingleSourceComplementTraversal(scratch);
        }
        return exactMultiSourceComplementTraversal(scratch);
    }

    /**
     * @brief Prepares one candidate and tests its residual-event eligibility.
     *
     * @param state Mutable construction state updated by the operation.
     * @param candidate Candidate processed by the operation.
     * @param statistics Statistics updated by the operation.
     * @return `true` when the candidate satisfies the configured eligibility policy.
     */
    [[nodiscard]] bool prepareAndTestEligibility(ConstructionState& state, const Candidate& candidate, MinMaxResidualTreeBuildStatistics& statistics) const {
        const tree_t& primal = candidate.polarity == Polarity::Max ? *state.maxTree : *state.minTree;
        const tree_t& dual = candidate.polarity == Polarity::Max ? *state.minTree : *state.maxTree;
        const MorphologicalTree& dualTopology = dual.topology();
        auto& scratch = state.certification;

        scratch.boundaryPixelMarks.resetAll();
        scratch.boundaryOwnerMarks.resetAll();
        scratch.supportOwnerMarks.resetAll();
        scratch.boundaryOwners.clear();
        scratch.boundaryPixels.clear();
        scratch.supportPixels = {};
        scratch.supportOwners.clear();
        scratch.dualExtremalOwner = InvalidNode;
        scratch.wholeSupportOwner = InvalidNode;
        altitude_t dualExtremalAltitude = altitude_t{};
        bool supportContainsInfinity = false;
        bool supportOwnerCertifiedFromFlatZone = false;
        const bool dualIsMaxTree = candidate.polarity == Polarity::Min;
        std::span<const NodeId> boundaryPixelsByIncidence;
        if (state.incrementalBoundary) {
            const auto view = state.flatZones.acquire(candidate.nodeId, primal);
            scratch.flatZoneRoot = view.root;
            scratch.supportPixels = view.supportPixels;
            boundaryPixelsByIncidence = view.externalPixelsByIncidence;
            supportContainsInfinity = view.containsInfinity;
            if (!scratch.supportPixels.empty() && (!supportContainsInfinity || !state.requiresSaturationCertification)) {
                const NodeId owner = dualTopology.getProperPartOwner(scratch.supportPixels.front());
                if (owner != InvalidNode && dualTopology.isAlive(owner) &&
                    dualTopology.getNumProperParts(owner) == static_cast<int>(scratch.supportPixels.size())) {
#ifndef NDEBUG
                    for (NodeId pixel : scratch.supportPixels) {
                        assert(dualTopology.getProperPartOwner(pixel) == owner);
                    }
#endif
                    scratch.supportOwners.push_back(owner);
                    scratch.dualExtremalOwner = owner;
                    scratch.wholeSupportOwner = owner;
                    supportOwnerCertifiedFromFlatZone = true;
                }
            }
        } else {
            RecomputedProperPartBoundary& boundaryCache = candidate.polarity == Polarity::Max ? state.maxBoundary : state.minBoundary;
            const auto& boundaryBag = boundaryCache.acquire(candidate.nodeId, primal, adjacency_);
            scratch.flatZoneRoot = InvalidNode;
            scratch.supportPixels = boundaryBag.supportPixels;
            boundaryPixelsByIncidence = boundaryBag.externalPixelsByIncidence;
        }
        if (!supportOwnerCertifiedFromFlatZone) {
            for (NodeId pixel : scratch.supportPixels) {
                supportContainsInfinity = supportContainsInfinity || pixel == infinityPixel_;
                const NodeId owner = dualTopology.getProperPartOwner(pixel);
                if (owner == InvalidNode || !dualTopology.isAlive(owner)) {
                    throw std::runtime_error("Min/max residual support has an invalid opposite-tree owner.");
                }
                if (!scratch.supportOwnerMarks.isMarked(static_cast<std::size_t>(owner))) {
                    scratch.supportOwnerMarks.mark(static_cast<std::size_t>(owner));
                    scratch.supportOwners.push_back(owner);
                    const altitude_t ownerAltitude = dual.getAltitude(owner);
                    if (scratch.dualExtremalOwner == InvalidNode || (dualIsMaxTree && ownerAltitude < dualExtremalAltitude) ||
                        (!dualIsMaxTree && ownerAltitude > dualExtremalAltitude)) {
                        scratch.dualExtremalOwner = owner;
                        dualExtremalAltitude = ownerAltitude;
                    }
                }
            }
        }
        if (!supportOwnerCertifiedFromFlatZone && scratch.supportOwners.size() == 1) {
            const NodeId owner = scratch.supportOwners.front();
            if (dualTopology.getNumProperParts(owner) == static_cast<int>(scratch.supportPixels.size())) {
                scratch.wholeSupportOwner = owner;
            }
        }
        if (state.requiresSaturationCertification && supportContainsInfinity) {
            return false;
        }

        for (NodeId neighbor : boundaryPixelsByIncidence) {
            const std::size_t neighborIndex = static_cast<std::size_t>(neighbor);
            if (!scratch.boundaryPixelMarks.isMarked(neighborIndex)) {
                scratch.boundaryPixelMarks.mark(neighborIndex);
                scratch.boundaryPixels.push_back(neighbor);
            }
            const NodeId owner = dualTopology.getProperPartOwner(neighbor);
            if (owner == InvalidNode || !dualTopology.isAlive(owner)) {
                throw std::runtime_error("Min/max residual boundary has an invalid opposite-tree owner.");
            }
            if (!scratch.boundaryOwnerMarks.isMarked(static_cast<std::size_t>(owner))) {
                scratch.boundaryOwnerMarks.mark(static_cast<std::size_t>(owner));
                scratch.boundaryOwners.push_back(owner);
            }
        }
        if (scratch.boundaryOwners.empty()) {
            throw std::runtime_error("A non-root min/max leaf has no external boundary.");
        }

        if (!state.requiresSaturationCertification) {
            return true;
        }

        NodeId common = scratch.boundaryOwners.front();
        DynamicLcaIndex& dualLca = candidate.polarity == Polarity::Max ? state.minLca : state.maxLca;
        for (std::size_t index = 1; index < scratch.boundaryOwners.size(); ++index) {
            common = lowestCommonAncestor(dual, dualLca, common, scratch.boundaryOwners[index]);
        }
        const altitude_t source = primal.getAltitude(candidate.nodeId);
        const altitude_t connectingLevel = dual.getAltitude(common);
        const bool certifiedByDualTree = candidate.polarity == Polarity::Max ? connectingLevel < source : connectingLevel > source;
        if (certifiedByDualTree) {
            return true;
        }

        bool commonContainsSupport = false;
        for (NodeId owner : scratch.supportOwners) {
            NodeId cursor = owner;
            while (true) {
                if (cursor == common) {
                    commonContainsSupport = true;
                    break;
                }
                const NodeId parent = dualTopology.getNodeParent(cursor);
                if (parent == cursor) {
                    break;
                }
                cursor = parent;
            }
            if (commonContainsSupport) {
                break;
            }
        }
        if (!commonContainsSupport) {
            return true;
        }

        ++statistics.complementTraversalCertificates;
        scratch.supportMarks.resetAll();
        for (NodeId pixel : scratch.supportPixels) {
            scratch.supportMarks.mark(static_cast<std::size_t>(pixel));
        }
        return exactComplementTraversal(scratch);
    }

    /**
     * @brief Commits one eligible extremum to the synchronized construction.
     *
     * @param state Mutable construction state updated by the operation.
     * @param candidate Candidate processed by the operation.
     */
    static void commitCandidate(ConstructionState& state, const Candidate& candidate) {
        tree_t& primal = candidate.polarity == Polarity::Max ? *state.maxTree : *state.minTree;
        tree_t& dual = candidate.polarity == Polarity::Max ? *state.minTree : *state.maxTree;
        const Polarity dualPolarity = candidate.polarity == Polarity::Max ? Polarity::Min : Polarity::Max;
        const NodeId parent = primal.topology().getNodeParent(candidate.nodeId);
        const auto support = std::span<const NodeId>(state.certification.supportPixels);
        const auto boundaryOwners = std::span<const NodeId>(state.certification.boundaryOwners);
        const NodeId dualExtremalOwner = state.certification.dualExtremalOwner;
        if (support.empty() || boundaryOwners.empty() || dualExtremalOwner == InvalidNode) {
            throw std::runtime_error("Min/max residual commit received an incomplete leaf certificate.");
        }
        const altitude_t targetLevel = primal.getAltitude(parent);
        if (state.incrementalBoundary) {
            state.flatZones.collectMergeRoots(state.certification.flatZoneRoot, targetLevel, state.certification.boundaryPixels,
                                              state.certification.flatZoneRootMarks, state.certification.flatZoneMergeRoots);
        }

        if (candidate.polarity == Polarity::Max) {
            state.adjustment->pruneMaxLeafAndUpdateMinTree(candidate.nodeId, support, dualExtremalOwner, boundaryOwners, state.certification.wholeSupportOwner);
        } else {
            state.adjustment->pruneMinLeafAndUpdateMaxTree(candidate.nodeId, support, dualExtremalOwner, boundaryOwners, state.certification.wholeSupportOwner);
        }
        if (state.incrementalBoundary) {
            const NodeId selectedFlatZoneRoot = state.certification.flatZoneRoot;
            const NodeId mergedFlatZoneRoot =
                state.flatZones.mergeAtLevel(state.certification.flatZoneRoot, state.certification.flatZoneMergeRoots, targetLevel);
            if (state.residualAssembler.has_value()) {
                state.absorbedFlatZones.clear();
                if (selectedFlatZoneRoot != mergedFlatZoneRoot) {
                    state.absorbedFlatZones.push_back(static_cast<detail::RegionId>(selectedFlatZoneRoot));
                }
                for (NodeId root : state.certification.flatZoneMergeRoots) {
                    if (root != mergedFlatZoneRoot) {
                        state.absorbedFlatZones.push_back(static_cast<detail::RegionId>(root));
                    }
                }
                state.residualAssembler->consume(static_cast<detail::RegionId>(mergedFlatZoneRoot), state.absorbedFlatZones);
            }
        }
        state.certification.supportPixels = {};

        if (state.requiresSaturationCertification) {
            DynamicLcaIndex& dualLca = candidate.polarity == Polarity::Max ? state.minLca : state.maxLca;
            dualLca.noteMutation(state.adjustment->getLastTopologyChangedNodes());
        }

        for (NodeId nodeId : state.adjustment->getLastCandidateNodes()) {
            if (state.agenda.contains(dualPolarity, nodeId) || isCurrentCandidateNode(dual, nodeId)) {
                updateAgendaCandidate(state, dual, dualPolarity, nodeId);
            }
        }
        updateAgendaCandidate(state, primal, candidate.polarity, candidate.nodeId);
        if (isCurrentCandidateNode(primal, parent)) {
            updateAgendaCandidate(state, primal, candidate.polarity, parent);
        }
    }

    /**
     * @brief Builds persistent residual-event chains for all pixels.
     *
     * @param state Mutable construction state updated by the operation.
     * @param statistics Statistics updated by the operation.
     * @return Common terminal altitude of the synchronized component-tree roots.
     */
    [[nodiscard]] altitude_t buildChains(ConstructionState& state, MinMaxResidualTreeBuildStatistics& statistics) const {
        int previousArea = 0;
        while (true) {
            const auto selected = state.agenda.select();
            if (!selected.has_value()) {
                break;
            }
            const Candidate candidate = *selected;
            const tree_t& primal = candidate.polarity == Polarity::Max ? *state.maxTree : *state.minTree;
            if (!isCandidateCurrent(state, primal, candidate)) {
                updateAgendaCandidate(state, primal, candidate.polarity, candidate.nodeId);
                continue;
            }

            const bool eligible = prepareAndTestEligibility(state, candidate, statistics);
            if (!eligible) {
                state.agenda.reject(candidate.polarity, candidate.nodeId);
                ++statistics.rejectedExtrema;
                continue;
            }
            if (candidate.area < previousArea) {
                throw std::runtime_error("Min/max residual event areas are not nondecreasing.");
            }
            previousArea = candidate.area;

            const altitude_t eventAltitude = primal.getAltitude(candidate.nodeId);
            if (state.residualAssembler.has_value()) {
                static_cast<void>(state.residualAssembler->emitEvent(static_cast<detail::RegionId>(state.certification.flatZoneRoot), eventAltitude));
            } else {
                recordEvent(state, eventAltitude, state.certification.supportPixels);
            }
            commitCandidate(state, candidate);
            ++statistics.residualEvents;
        }

        if (state.maxTree->topology().getNumNodes() != 1 || state.minTree->topology().getNumNodes() != 1) {
            throw std::runtime_error("Min/max residual agenda was exhausted before both component trees became constant.");
        }
        const NodeId maxRoot = state.maxTree->topology().getRoot();
        const NodeId minRoot = state.minTree->topology().getRoot();
        const altitude_t maxAltitude = state.maxTree->getAltitude(maxRoot);
        const altitude_t minAltitude = state.minTree->getAltitude(minRoot);
        if (maxAltitude != minAltitude) {
            throw std::runtime_error("Min/max residual terminal roots have different altitudes.");
        }
        return maxAltitude;
    }

    /**
     * @brief Materializes hierarchy storage from persistent residual-event chains.
     *
     * @param image Image processed by the operation.
     * @param terminalAltitude Altitude assigned to the terminal residual node.
     * @param chains Chains used by the operation.
     */
    void materializeChains(const image_ptr_t& image, altitude_t terminalAltitude, const PersistentChains& chains) {
        rows_ = image->getNumRows();
        cols_ = image->getNumCols();
        root_ = 0;
        nodeParent_.assign(chains.eventValuations.size() + 1, InvalidNode);
        nodeParent_[0] = 0;
        properPartOwner_.assign(chains.stackHeadByPixel.size(), 0);
        altitude_.assign(chains.eventValuations.size() + 1, altitude_t{});
        altitude_[0] = terminalAltitude;
        for (std::size_t event = 0; event < chains.eventValuations.size(); ++event) {
            altitude_[event + 1] = chains.eventValuations[event];
        }

        std::vector<NodeId> parentByEvent(chains.eventValuations.size(), InvalidNode);
        std::vector<int> ownerEventByHead(chains.nodes.size(), -2);
        std::vector<int> headStack;
        headStack.reserve(32);
        for (int head : chains.stackHeadByPixel) {
            if (head < 0) {
                continue;
            }
            if (static_cast<std::size_t>(head) >= chains.nodes.size()) {
                throw std::runtime_error("Min/max residual chain head lies outside the arena.");
            }
            if (ownerEventByHead[static_cast<std::size_t>(head)] != -2) {
                continue;
            }
            int cursor = head;
            while (cursor >= 0 && ownerEventByHead[static_cast<std::size_t>(cursor)] == -2) {
                headStack.push_back(cursor);
                cursor = chains.nodes[static_cast<std::size_t>(cursor)].prev;
            }
            int ownerEvent = cursor >= 0 ? ownerEventByHead[static_cast<std::size_t>(cursor)] : -1;
            while (!headStack.empty()) {
                const int current = headStack.back();
                headStack.pop_back();
                const auto& node = chains.nodes[static_cast<std::size_t>(current)];
                ownerEvent = ownerEvent >= 0 ? ownerEvent : node.eventId;
                ownerEventByHead[static_cast<std::size_t>(current)] = ownerEvent;
                if (node.prev >= 0) {
                    const int childEvent = chains.nodes[static_cast<std::size_t>(node.prev)].eventId;
                    const NodeId parent = static_cast<NodeId>(node.eventId) + 1;
                    auto& stored = parentByEvent[static_cast<std::size_t>(childEvent)];
                    if (stored == InvalidNode) {
                        stored = parent;
                    } else if (stored != parent) {
                        throw std::runtime_error("Min/max residual chains imply inconsistent parents.");
                    }
                }
            }
        }

        for (NodeId pixel = 0; pixel < static_cast<NodeId>(chains.stackHeadByPixel.size()); ++pixel) {
            const int head = chains.stackHeadByPixel[static_cast<std::size_t>(pixel)];
            if (head < 0) {
                properPartOwner_[static_cast<std::size_t>(pixel)] = 0;
                continue;
            }
            const int ownerEvent = ownerEventByHead[static_cast<std::size_t>(head)];
            properPartOwner_[static_cast<std::size_t>(pixel)] = ownerEvent < 0 ? 0 : static_cast<NodeId>(ownerEvent) + 1;
        }
        for (std::size_t event = 0; event < parentByEvent.size(); ++event) {
            nodeParent_[event + 1] = parentByEvent[event] == InvalidNode ? 0 : parentByEvent[event];
        }

        const NativeHierarchyView<altitude_t> hierarchy{
            nodeParent_, properPartOwner_,           altitude_,
            0,           GridDomain2D{rows_, cols_}, makeHierarchySemantics(MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE, adjacency_)};
        topologyProof_ = mmcfilters::detail::NativeHierarchyValidation::validateComplete(hierarchy);
        for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
            const NodeId owner = properPartOwner_[static_cast<std::size_t>(pixel)];
            if (altitude_[static_cast<std::size_t>(owner)] != (*image)[pixel]) {
                throw std::runtime_error("Min/max residual reconstruction differs from the input image.");
            }
        }
    }

    /**
     * @brief Materializes validated hierarchy storage from assembler output.
     *
     * @param image Image processed by the operation.
     * @param output Output storage consumed or populated by the operation.
     */
    void materializeAssembler(const image_ptr_t& image, typename detail::UnionFindResidualTreeAssembler<altitude_t>::Output output) {
        const NativeHierarchyView<altitude_t> hierarchy{output.nodeParent,
                                                        output.properPartOwner,
                                                        output.altitude,
                                                        0,
                                                        GridDomain2D{image->getNumRows(), image->getNumCols()},
                                                        makeHierarchySemantics(MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE, adjacency_)};
        auto topologyProof = mmcfilters::detail::NativeHierarchyValidation::validateWithEstablishedSupport(
            hierarchy, std::move(output.subtreeSupportProof), [&image, &output](NodeId properPart, NodeId owner) {
                if (output.altitude[static_cast<std::size_t>(owner)] != (*image)[properPart]) {
                    throw std::runtime_error("Min/max residual assembler reconstruction differs from the input image.");
                }
            });
        rows_ = image->getNumRows();
        cols_ = image->getNumCols();
        root_ = 0;
        nodeParent_ = std::move(output.nodeParent);
        properPartOwner_ = std::move(output.properPartOwner);
        altitude_ = std::move(output.altitude);
        topologyProof_ = std::move(topologyProof);
    }

    /** @brief Requires a successfully built result to be available. */
    void requireBuilt() const {
        if (!built_) {
            throw std::logic_error("MinMaxResidualTreeBuilder accessors require a completed build.");
        }
    }

  public:
    /**
     * @brief Configures one reusable synchronized min/max construction.
     *
     * @param adjacency Symmetric adjacency used by both component trees and,
     *        in saturated mode, by complement connectivity.
     * @param infinityPixel Row-major exterior seed excluded from residual
     *        events only in saturated mode.
     * @param tiePolicy Equal-area ordering policy.
     * @param lcaPolicy Dynamic LCA strategy.
     * @param fallbackPolicy Exact complement fallback strategy.
     * @param boundaryPolicy Proper-part boundary maintenance strategy.
     * @param eligibilityPolicy Whether to require saturated supports or accept
     *        every current regional extremum.
     */
    explicit MinMaxResidualTreeBuilder(RegularGridAdjacency2D adjacency, NodeId infinityPixel,
                                       SdrtTiePolicy tiePolicy = SdrtTiePolicy::ContrastInvariantSpatial,
                                       SaturatedMinMaxLcaPolicy lcaPolicy = SaturatedMinMaxLcaPolicy::ParentClimb,
                                       SaturatedMinMaxFallbackPolicy fallbackPolicy = SaturatedMinMaxFallbackPolicy::BoundaryMultiSource,
                                       SaturatedMinMaxBoundaryPolicy boundaryPolicy = SaturatedMinMaxBoundaryPolicy::IncrementalSmallToLarge,
                                       MinMaxResidualEligibilityPolicy eligibilityPolicy = MinMaxResidualEligibilityPolicy::SaturatedOnly)
        : adjacency_(std::move(adjacency)), infinityPixel_(infinityPixel), tiePolicy_(tiePolicy), lcaPolicy_(lcaPolicy), fallbackPolicy_(fallbackPolicy),
          boundaryPolicy_(boundaryPolicy), eligibilityPolicy_(eligibilityPolicy) {}

    /**
     * @brief Consumes caller-built component-tree seeds.
     *
     * The caller is responsible for constructing both seeds with the same
     * adjacency configured in this builder. This overload validates domains,
     * tree kinds, and altitudes, but component trees do not retain enough
     * metadata to prove equality of arbitrary stencils.
     *
     * @param image Image processed by the operation.
     * @param minTree Min-tree consumed by the operation.
     * @param maxTree Max-tree consumed by the operation.
     */
    void build(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree) {
        clearResult();
        MinMaxResidualTreeBuildStatistics localStatistics;
        requireInputs(image);
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(image, "MinMaxResidualTreeBuilder input image");
        buildFromSeeds(image, std::move(minTree), std::move(maxTree), localStatistics);
    }

    /**
     * @brief Returns the number of image rows from the last build.
     *
     * @return The rows.
     */
    [[nodiscard]] int getRows() const {
        requireBuilt();
        return rows_;
    }

    /**
     * @brief Returns the number of image columns from the last build.
     *
     * @return The cols.
     */
    [[nodiscard]] int getCols() const {
        requireBuilt();
        return cols_;
    }

    /**
     * @brief Returns the residual root identifier from the last build.
     *
     * @return The root.
     */
    [[nodiscard]] NodeId getRoot() const {
        requireBuilt();
        return root_;
    }

    /**
     * @brief Returns the configured row-major exterior seed.
     *
     * @return The infinity pixel.
     */
    [[nodiscard]] NodeId getInfinityPixel() const noexcept { return infinityPixel_; }

    /**
     * @brief Returns the shared adjacency configured for both trees.
     *
     * @return The adjacency.
     */
    [[nodiscard]] const RegularGridAdjacency2D& getAdjacency() const noexcept { return adjacency_; }

    /**
     * @brief Returns the configured deterministic tie policy.
     *
     * @return The tie policy.
     */
    [[nodiscard]] SdrtTiePolicy getTiePolicy() const noexcept { return tiePolicy_; }

    /**
     * @brief Returns the configured saturation-certificate LCA policy.
     *
     * @return The lowest-common-ancestor policy.
     */
    [[nodiscard]] SaturatedMinMaxLcaPolicy getLcaPolicy() const noexcept { return lcaPolicy_; }

    /**
     * @brief Returns the configured exact complement fallback policy.
     *
     * @return The fallback policy.
     */
    [[nodiscard]]
    SaturatedMinMaxFallbackPolicy getFallbackPolicy() const noexcept {
        return fallbackPolicy_;
    }

    /**
     * @brief Returns the configured boundary-maintenance policy.
     *
     * @return The boundary policy.
     */
    [[nodiscard]]
    SaturatedMinMaxBoundaryPolicy getBoundaryPolicy() const noexcept {
        return boundaryPolicy_;
    }

    /**
     * @brief Returns the configured extremum eligibility policy.
     *
     * @return The eligibility policy.
     */
    [[nodiscard]]
    MinMaxResidualEligibilityPolicy getEligibilityPolicy() const noexcept {
        return eligibilityPolicy_;
    }

    /**
     * @brief Returns the residual parent buffer from the last build.
     *
     * @return The node parent.
     */
    [[nodiscard]] std::span<const NodeId> getNodeParent() const {
        requireBuilt();
        return nodeParent_;
    }

    /**
     * @brief Returns the direct proper-part owner of every pixel.
     *
     * @return The proper part owner.
     */
    [[nodiscard]] std::span<const NodeId> getProperPartOwner() const {
        requireBuilt();
        return properPartOwner_;
    }

    /**
     * @brief Returns the residual valuation buffer from the last build.
     *
     * @return The altitude.
     */
    [[nodiscard]] std::span<const altitude_t> getAltitude() const {
        requireBuilt();
        return altitude_;
    }

    /**
     * @brief Returns diagnostics collected by the last successful build.
     *
     * @return The statistics.
     */
    [[nodiscard]] const MinMaxResidualTreeBuildStatistics& getStatistics() const {
        requireBuilt();
        return statistics_;
    }

    /**
     * @brief Transfers the last result into a validated native hierarchy.
     * @param semantics Semantic capabilities assigned to the result.
     * @return Validated hierarchy storage owning all residual buffers.
     */
    [[nodiscard]]
    mmcfilters::detail::ValidatedNativeHierarchy<altitude_t> takeValidatedHierarchy(HierarchySemantics semantics) && {
        requireBuilt();
        auto hierarchy =
            mmcfilters::detail::makeValidatedNativeHierarchy<altitude_t>(std::move(nodeParent_), std::move(properPartOwner_), std::move(altitude_), root_,
                                                                         GridDomain2D{rows_, cols_}, std::move(semantics), std::move(topologyProof_));
        built_ = false;
        rows_ = 0;
        cols_ = 0;
        root_ = InvalidNode;
        return hierarchy;
    }

  private:
    /**
     * @brief Builds the residual tree from synchronized min-tree and max-tree seeds.
     *
     * @param image Image processed by the operation.
     * @param minTree Min-tree consumed by the operation.
     * @param maxTree Max-tree consumed by the operation.
     * @param localStatistics Mutable local statistics updated by the operation.
     */
    void buildFromSeeds(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree, MinMaxResidualTreeBuildStatistics& localStatistics) {
        ConstructionState state = initializeConstruction(image, std::move(minTree), std::move(maxTree));
        const altitude_t terminalAltitude = buildChains(state, localStatistics);
        if (state.residualAssembler.has_value()) {
            if (state.residualAssembler->numEvents() != localStatistics.residualEvents) {
                throw std::runtime_error("Min/max residual incremental assembly changed the residual event count.");
            }
            const NodeId terminalFlatZone = state.flatZones.currentRoot(infinityPixel_);
            materializeAssembler(image, state.residualAssembler->finalize(static_cast<detail::RegionId>(terminalFlatZone), terminalAltitude));
        } else {
            materializeChains(image, terminalAltitude, state.chains);
        }
        statistics_ = localStatistics;
        built_ = true;
    }
};

} // namespace mmcfilters::sdrt
