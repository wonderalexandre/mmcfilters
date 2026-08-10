#pragma once

/**
 * @file MinMaxResidualTreeEngineCommon.hpp
 * @brief Shared implementation of specialized synchronized min/max residual-tree engines.
 */

#include "../ResidualTreeBuildStatistics.hpp"
#include "../ResidualTreePolicies.hpp"
#include "ResidualTreeCandidateAgenda.hpp"
#include "ResidualTreeCandidateContext.hpp"
#include "ResidualTreeBoundaryCache.hpp"
#include "ResidualTreeEventAssembler.hpp"
#include "SaturatedDynamicLca.hpp"
#include "SaturatedResidualEligibility.hpp"
#include "../../adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "../../detail/NativeHierarchyValidationDetail.hpp"
#include "../../../utils/GenerationStampSet.hpp"

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
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Shared synchronized min/max engine for the specialized residual-tree builders.
 *
 * Both component trees use the same symmetric adjacency. The saturated
 * specialization accepts a non-exterior leaf only when its current
 * proper-part support has connected complement; the unrestricted
 * specialization accepts every current non-root regional extremum. Both use
 * the same agenda, updater and output assembly. Accepted max-tree leaves are
 * pruned while the min-tree is updated incrementally, and conversely for
 * accepted min-tree leaves.
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
 * The saturated and unrestricted builders instantiate this implementation
 * detail with different eligibility rules. Because that rule is a template
 * argument, saturation-only branches are discarded from the unrestricted
 * specialization.
 *
 * @tparam T Finite image, component-tree, and residual altitude type.
 * @tparam RequiresSaturationCertification Whether complement connectivity is required.
 * @internal
 */
template <AltitudeValue T, bool RequiresSaturationCertification> class MinMaxResidualTreeEngine {
  private:
    /** @brief Defines the `altitude_t` alias used by the component. */
    using altitude_t = T;
    /** @brief Defines the `image_ptr_t` alias used by the component. */
    using image_ptr_t = ImagePtr<altitude_t>;
    /** @brief Defines the `tree_t` alias used by the component. */
    using tree_t = WeightedMorphologicalTree<altitude_t>;
    /** @brief Defines the `Adjustment` alias used by the component. */
    using Adjustment = adjust::DualMinMaxTreeIncrementalFilterLeaf<altitude_t>;
    /** @brief Candidate-agenda type family for the current altitude type. */
    using AgendaTypes = ResidualTreeAgendaTypes<altitude_t>;
    /** @brief Component-tree polarity handled by a candidate. */
    using Polarity = typename AgendaTypes::Polarity;
    /** @brief Candidate record stored in the deterministic agenda. */
    using Candidate = typename AgendaTypes::Candidate;
    /** @brief Derived support properties cached for a candidate. */
    using CandidateDescriptor = typename AgendaTypes::CandidateDescriptor;
    /** @brief Deterministic agenda of current component-tree leaves. */
    using CandidateAgenda = typename AgendaTypes::CandidateAgenda;

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
    /** @brief Dynamic-LCA type family for saturated certification. */
    using LcaTypes = SaturatedLcaTypes<altitude_t>;
    /** @brief Policy-selected dynamic LCA index. */
    using DynamicLcaIndex = typename LcaTypes::DynamicLcaIndex;
    /** @brief Boundary-storage type family for the current altitude type. */
    using BoundaryTypes = ResidualTreeBoundaryTypes<altitude_t>;
    /** @brief On-demand proper-part boundary implementation. */
    using RecomputedProperPartBoundary = typename BoundaryTypes::RecomputedProperPartBoundary;
    /** @brief Incrementally merged flat-zone boundary implementation. */
    using IncrementalFlatZoneBoundaryCache = typename BoundaryTypes::IncrementalFlatZoneBoundaryCache;

    /** @brief Shared candidate support and boundary scratch. */
    using CandidateContextScratch = ResidualTreeCandidateContext;
    /** @brief Saturated-only complement-connectivity state. */
    using SaturationCertificationState = SaturatedResidualEligibilityState<altitude_t>;
    /** @brief Empty eligibility state for unrestricted construction. */
    using UnrestrictedEligibilityState = UnrestrictedResidualEligibilityState;
    /** @brief Eligibility state selected by the engine specialization. */
    using EligibilityState =
        std::conditional_t<RequiresSaturationCertification, SaturationCertificationState, UnrestrictedEligibilityState>;

    /** @brief Owns mutable state for one synchronized min/max residual-tree build. */
    struct ConstructionState {
        /** @brief Stores the max tree. */
        std::unique_ptr<tree_t> maxTree;
        /** @brief Stores the min tree. */
        std::unique_ptr<tree_t> minTree;
        /** @brief Stores the adjustment. */
        std::unique_ptr<Adjustment> adjustment;
        /** @brief Stores the agenda. */
        CandidateAgenda agenda;
        /** @brief Stores the chains. */
        PersistentChains chains;
        /** @brief Stores candidate data shared by both modes. */
        CandidateContextScratch candidateContext;
        /** @brief Stores mode-specific eligibility state. */
        EligibilityState eligibility;
        /** @brief Stores the max boundary. */
        RecomputedProperPartBoundary maxBoundary;
        /** @brief Stores the min boundary. */
        RecomputedProperPartBoundary minBoundary;
        /** @brief Stores the flat zones. */
        IncrementalFlatZoneBoundaryCache flatZones;
        /** @brief Stores the residual assembler. */
        std::optional<ResidualTreeEventAssembler<altitude_t>> residualAssembler;
        /** @brief Stores the absorbed flat zones. */
        std::vector<RegionId> absorbedFlatZones;
        /** @brief Indicates whether incremental boundary maintenance is enabled. */
        bool incrementalBoundary = true;
        /** @brief Stores the cached new head by prev head. */
        std::vector<int> cachedNewHeadByPrevHead;
        /** @brief Stores the cached new head stamp by prev head. */
        std::vector<std::uint32_t> cachedNewHeadStampByPrevHead;
        /** @brief Stores the current event cache stamp. */
        std::uint32_t currentEventCacheStamp = 1;
        /**
         * @brief Constructs a `ConstructionState` instance.
         *
         * @param tiePolicy Policy used to resolve deterministic ties.
         * @param infinityPixel Row-major pixel used as the exterior seed.
         * @param numPixels Num pixels used by the operation.
         * @param maxNodeSlots Max node slots used by the operation.
         * @param incrementalBoundary Incremental boundary used by the operation.
         * @param adjacency Established symmetric pixel adjacency.
         * @param fallbackPolicy Exact complement traversal policy for saturated construction.
         */
        ConstructionState(SdrtTiePolicy tiePolicy, NodeId infinityPixel, std::size_t numPixels, std::size_t maxNodeSlots, bool incrementalBoundary,
                          const RegularGridAdjacency2D& adjacency, SaturatedMinMaxFallbackPolicy fallbackPolicy)
            : agenda(tiePolicy, infinityPixel, RequiresSaturationCertification), candidateContext(numPixels, maxNodeSlots),
              eligibility(numPixels, adjacency, infinityPixel, fallbackPolicy),
              flatZones(incrementalBoundary), incrementalBoundary(incrementalBoundary) {}
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
    ResidualTreeBoundaryPolicy boundaryPolicy_ = ResidualTreeBoundaryPolicy::IncrementalSmallToLarge;
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
    ResidualTreeBuildStatistics statistics_;

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
        statistics_ = ResidualTreeBuildStatistics{};
    }

    /**
     * @brief Validates the inputs required by the construction.
     *
     * @param image Image processed by the operation.
     */
    void requireInputs(const image_ptr_t& image) const {
        if (!image || image->getNumRows() <= 0 || image->getNumCols() <= 0 || image->getSize() <= 0) {
            throw std::invalid_argument("Residual-tree engine requires a non-null, non-empty image.");
        }
        if (adjacency_.getNumRows() != image->getNumRows() || adjacency_.getNumCols() != image->getNumCols()) {
            throw std::invalid_argument("Residual-tree engine adjacency domain differs from the image.");
        }
        if (infinityPixel_ < 0 || infinityPixel_ >= image->getSize()) {
            throw std::invalid_argument("Residual-tree engine infinity pixel lies outside the image.");
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
            throw std::invalid_argument("Residual-tree engine requires a connected adjacency domain.");
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
                                boundaryPolicy_ == ResidualTreeBoundaryPolicy::IncrementalSmallToLarge, adjacency_, fallbackPolicy_);
        state.maxTree = std::make_unique<tree_t>(std::move(maxTree));
        state.minTree = std::make_unique<tree_t>(std::move(minTree));
        state.adjustment = std::make_unique<Adjustment>(state.minTree.get(), state.maxTree.get(), adjacency_);
        if constexpr (RequiresSaturationCertification) {
            state.eligibility.minLca.bind(*state.minTree, lcaPolicy_);
            state.eligibility.maxLca.bind(*state.maxTree, lcaPolicy_);
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
     * @brief Prepares one candidate and tests its residual-event eligibility.
     *
     * @param state Mutable construction state updated by the operation.
     * @param candidate Candidate processed by the operation.
     * @param statistics Statistics updated by the operation.
     * @return `true` when the candidate satisfies the configured eligibility policy.
     */
    [[nodiscard]] bool prepareAndTestEligibility(ConstructionState& state, const Candidate& candidate, ResidualTreeBuildStatistics& statistics) const {
        const tree_t& primal = candidate.polarity == Polarity::Max ? *state.maxTree : *state.minTree;
        const tree_t& dual = candidate.polarity == Polarity::Max ? *state.minTree : *state.maxTree;
        const MorphologicalTree& dualTopology = dual.topology();
        auto& scratch = state.candidateContext;

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
            if (!scratch.supportPixels.empty() && (!supportContainsInfinity || !RequiresSaturationCertification)) {
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
        if constexpr (RequiresSaturationCertification) {
            if (supportContainsInfinity) {
                return false;
            }
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

        if constexpr (!RequiresSaturationCertification) {
            return true;
        } else {
            NodeId common = scratch.boundaryOwners.front();
            DynamicLcaIndex& dualLca =
                candidate.polarity == Polarity::Max ? state.eligibility.minLca : state.eligibility.maxLca;
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
            auto& saturation = state.eligibility;
            saturation.supportMarks.resetAll();
            for (NodeId pixel : scratch.supportPixels) {
                saturation.supportMarks.mark(static_cast<std::size_t>(pixel));
            }
            return saturation.exactComplementTraversal(scratch);
        }
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
        const auto support = std::span<const NodeId>(state.candidateContext.supportPixels);
        const auto boundaryOwners = std::span<const NodeId>(state.candidateContext.boundaryOwners);
        const NodeId dualExtremalOwner = state.candidateContext.dualExtremalOwner;
        if (support.empty() || boundaryOwners.empty() || dualExtremalOwner == InvalidNode) {
            throw std::runtime_error("Min/max residual commit received an incomplete leaf certificate.");
        }
        const altitude_t targetLevel = primal.getAltitude(parent);
        if (state.incrementalBoundary) {
            state.flatZones.collectMergeRoots(state.candidateContext.flatZoneRoot, targetLevel, state.candidateContext.boundaryPixels,
                                              state.candidateContext.flatZoneRootMarks, state.candidateContext.flatZoneMergeRoots);
        }

        if (candidate.polarity == Polarity::Max) {
            state.adjustment->pruneMaxLeafAndUpdateMinTree(candidate.nodeId, support, dualExtremalOwner, boundaryOwners,
                                                          state.candidateContext.wholeSupportOwner);
        } else {
            state.adjustment->pruneMinLeafAndUpdateMaxTree(candidate.nodeId, support, dualExtremalOwner, boundaryOwners,
                                                          state.candidateContext.wholeSupportOwner);
        }
        if (state.incrementalBoundary) {
            const NodeId selectedFlatZoneRoot = state.candidateContext.flatZoneRoot;
            const NodeId mergedFlatZoneRoot =
                state.flatZones.mergeAtLevel(state.candidateContext.flatZoneRoot, state.candidateContext.flatZoneMergeRoots, targetLevel);
            if (state.residualAssembler.has_value()) {
                state.absorbedFlatZones.clear();
                if (selectedFlatZoneRoot != mergedFlatZoneRoot) {
                    state.absorbedFlatZones.push_back(static_cast<RegionId>(selectedFlatZoneRoot));
                }
                for (NodeId root : state.candidateContext.flatZoneMergeRoots) {
                    if (root != mergedFlatZoneRoot) {
                        state.absorbedFlatZones.push_back(static_cast<RegionId>(root));
                    }
                }
                state.residualAssembler->consume(static_cast<RegionId>(mergedFlatZoneRoot), state.absorbedFlatZones);
            }
        }
        state.candidateContext.supportPixels = {};

        if constexpr (RequiresSaturationCertification) {
            DynamicLcaIndex& dualLca =
                candidate.polarity == Polarity::Max ? state.eligibility.minLca : state.eligibility.maxLca;
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
    [[nodiscard]] altitude_t buildChains(ConstructionState& state, ResidualTreeBuildStatistics& statistics) const {
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
                static_cast<void>(state.residualAssembler->emitEvent(static_cast<RegionId>(state.candidateContext.flatZoneRoot), eventAltitude));
            } else {
                recordEvent(state, eventAltitude, state.candidateContext.supportPixels);
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
    void materializeAssembler(const image_ptr_t& image, typename ResidualTreeEventAssembler<altitude_t>::Output output) {
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
            throw std::logic_error("Residual-tree engine accessors require a completed build.");
        }
    }

  public:
    /**
     * @brief Configures one reusable synchronized min/max construction.
     *
     * @param adjacency Symmetric adjacency used by both component trees and,
     *        in the saturated specialization, by complement connectivity.
     * @param infinityPixel Row-major exterior seed excluded from residual
     *        events only by the saturated specialization.
     * @param options Saturation and ordering policies.
     */
    explicit MinMaxResidualTreeEngine(RegularGridAdjacency2D adjacency, NodeId infinityPixel,
                                      SaturatedResidualTreeOptions options)
        requires(RequiresSaturationCertification)
        : adjacency_(std::move(adjacency)), infinityPixel_(infinityPixel), tiePolicy_(options.tiePolicy), lcaPolicy_(options.lcaPolicy),
          fallbackPolicy_(options.fallbackPolicy), boundaryPolicy_(options.boundaryPolicy) {}

    /**
     * @brief Configures unrestricted construction without saturation-only parameters.
     * @param adjacency Symmetric adjacency used by both component trees.
     * @param options Unrestricted ordering policies.
     */
    explicit MinMaxResidualTreeEngine(RegularGridAdjacency2D adjacency, UnrestrictedResidualTreeOptions options)
        requires(!RequiresSaturationCertification)
        : adjacency_(std::move(adjacency)), infinityPixel_(NodeId{0}), tiePolicy_(options.tiePolicy),
          boundaryPolicy_(ResidualTreeBoundaryPolicy::IncrementalSmallToLarge) {}

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
        ResidualTreeBuildStatistics localStatistics;
        requireInputs(image);
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(image, "Residual-tree engine input image");
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
    ResidualTreeBoundaryPolicy getBoundaryPolicy() const noexcept {
        return boundaryPolicy_;
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
    [[nodiscard]] const ResidualTreeBuildStatistics& getStatistics() const {
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
    void buildFromSeeds(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree, ResidualTreeBuildStatistics& localStatistics) {
        ConstructionState state = initializeConstruction(image, std::move(minTree), std::move(maxTree));
        const altitude_t terminalAltitude = buildChains(state, localStatistics);
        if (state.residualAssembler.has_value()) {
            if (state.residualAssembler->numEvents() != localStatistics.residualEvents) {
                throw std::runtime_error("Min/max residual incremental assembly changed the residual event count.");
            }
            const NodeId terminalFlatZone = state.flatZones.currentRoot(infinityPixel_);
            materializeAssembler(image, state.residualAssembler->finalize(static_cast<RegionId>(terminalFlatZone), terminalAltitude));
        } else {
            materializeChains(image, terminalAltitude, state.chains);
        }
        statistics_ = localStatistics;
        built_ = true;
    }
};

} // namespace mmcfilters::sdrt::detail
