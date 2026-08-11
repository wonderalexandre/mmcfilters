#pragma once

/**
 * @file MinMaxResidualTreeEngine.hpp
 * @brief Synchronized min/max contraction engine for residual-tree construction.
 */

#include "../ResidualTreeBuildStatistics.hpp"
#include "../ResidualTreePolicies.hpp"
#include "FlatZonePartition.hpp"
#include "ResidualTreeCandidateAgenda.hpp"
#include "ResidualTreeCandidateContext.hpp"
#include "ResidualTreeCandidatePreparation.hpp"
#include "ResidualTreeEventAssembler.hpp"
#include "ResidualTreeMaterialization.hpp"
#include "SaturatedResidualEligibility.hpp"
#include "../../TreeAltitudeAlgorithms.hpp"
#include "../../adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Synchronizes min-tree and max-tree contractions and emits residual events.
 *
 * Both specializations share candidate ordering, preparation, contraction, and
 * event assembly. The saturated specialization additionally owns a dedicated
 * eligibility evaluator; saturation-only branches and state are absent from
 * the unrestricted specialization.
 *
 * @tparam T Finite image, component-tree, and residual altitude type.
 * @tparam RequiresSaturationCertification Whether complement connectivity is required.
 * @internal
 */
template <AltitudeValue T, bool RequiresSaturationCertification> class MinMaxResidualTreeEngine {
  private:
    using altitude_t = T;                                                       ///< Image and residual altitude type.
    using image_ptr_t = ImagePtr<altitude_t>;                                   ///< Shared input image type.
    using tree_t = WeightedMorphologicalTree<altitude_t>;                       ///< Mutable component-tree type.
    using Adjustment = adjust::DualMinMaxTreeIncrementalFilterLeaf<altitude_t>; ///< Synchronized tree adjustment.
    using AgendaTypes = ResidualTreeAgendaTypes;                                ///< Candidate agenda type family.
    using Polarity = typename AgendaTypes::Polarity;                            ///< Candidate polarity.
    using Candidate = typename AgendaTypes::Candidate;                          ///< Ordered candidate record.
    using CandidateDescriptor = typename AgendaTypes::CandidateDescriptor;      ///< Partition-supplied metadata.
    using CandidateAgenda = typename AgendaTypes::CandidateAgenda;              ///< Deterministic active agenda.
    using FlatZonePartitionType = FlatZonePartition<altitude_t>;                ///< Current flat-zone partition.
    using SaturatedEvaluator = SaturatedEligibilityEvaluator<altitude_t>;       ///< Saturated-only eligibility evaluator.

    /** @brief Configuration retained only by the saturated specialization. */
    struct SaturatedConfiguration {
        NodeId infinityPixel = InvalidNode; ///< Exterior seed pixel.
        SaturatedMinMaxLcaPolicy lcaPolicy = SaturatedMinMaxLcaPolicy::ParentClimb; ///< Dynamic LCA strategy.
        SaturatedMinMaxFallbackPolicy fallbackPolicy = SaturatedMinMaxFallbackPolicy::BoundaryMultiSource; ///< Exact fallback.
    };

    /** @brief Empty mode configuration for unrestricted construction. */
    struct UnrestrictedConfiguration {};

    /** @brief Empty evaluator storage for unrestricted construction. */
    struct NoSaturationEvaluator {};

    using ModeConfiguration =
        std::conditional_t<RequiresSaturationCertification, SaturatedConfiguration, UnrestrictedConfiguration>; ///< Mode-specific options.
    using SaturationStorage =
        std::conditional_t<RequiresSaturationCertification, std::optional<SaturatedEvaluator>, NoSaturationEvaluator>; ///< Mode-specific evaluator storage.

    /** @brief Owns mutable state for one synchronized construction. */
    struct ConstructionState {
        std::unique_ptr<tree_t> maxTree;             ///< Current max-tree.
        std::unique_ptr<tree_t> minTree;             ///< Current min-tree.
        std::unique_ptr<Adjustment> adjustment;      ///< Synchronized tree updater.
        CandidateAgenda agenda;                      ///< Active deterministic candidate agenda.
        ResidualTreeCandidateContext candidateContext; ///< Reusable preparation scratch.
        SaturationStorage saturatedEligibility;      ///< Saturated-only evaluator storage.
        FlatZonePartitionType flatZonePartition;     ///< Current flat-zone union partition.
        ResidualTreeEventAssembler<altitude_t> residualAssembler; ///< Sequential residual event assembler.
        std::vector<RegionId> absorbedFlatZones;     ///< Flat zones absorbed by the current contraction.

        /**
         * @brief Creates mutable state around an initialized flat-zone partition.
         * @param tiePolicy Deterministic equal-area ordering policy.
         * @param numPixels Number of pixels in the common image domain.
         * @param maxNodeSlots Maximum node-slot capacity across both component trees.
         * @param initializedPartition Initial flat-zone partition transferred into the state.
         * @param initialRegionByPixel Initial flat-zone representative indexed by pixel.
         */
        ConstructionState(SdrtTiePolicy tiePolicy, std::size_t numPixels, std::size_t maxNodeSlots,
                          FlatZonePartitionType&& initializedPartition, const std::vector<RegionId>& initialRegionByPixel)
            : agenda(tiePolicy, RequiresSaturationCertification), candidateContext(numPixels, maxNodeSlots),
              flatZonePartition(std::move(initializedPartition)), residualAssembler(numPixels, initialRegionByPixel) {}
    };

    RegularGridAdjacency2D adjacency_; ///< Shared symmetric adjacency.
    SdrtTiePolicy tiePolicy_ = SdrtTiePolicy::ContrastInvariantSpatial; ///< Deterministic equal-area order.
    ModeConfiguration modeConfiguration_; ///< Mode-specific configuration without unrestricted exterior state.
    bool built_ = false;                    ///< Whether a completed result is available.
    int rows_ = 0;                          ///< Rows of the completed result.
    int cols_ = 0;                          ///< Columns of the completed result.
    NodeId root_ = InvalidNode;             ///< Root of the completed result.
    std::vector<NodeId> nodeParent_;        ///< Residual parent buffer.
    std::vector<NodeId> properPartOwner_;   ///< Direct proper-part owner buffer.
    std::vector<altitude_t> altitude_;      ///< Residual altitude buffer.
    mmcfilters::detail::NativeTopologyProof topologyProof_; ///< Established topology proof.
    ResidualTreeBuildStatistics statistics_; ///< Diagnostics from the completed build.

    /** @brief Clears all state associated with the previous result. */
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
     * @brief Validates the common image/adjacency contract and saturated exterior seed.
     * @param image Image whose domain and connectivity are validated.
     */
    void requireInputs(const image_ptr_t& image) const {
        if (!image || image->getNumRows() <= 0 || image->getNumCols() <= 0 || image->getSize() <= 0) {
            throw std::invalid_argument("Residual-tree engine requires a non-null, non-empty image.");
        }
        if (adjacency_.getNumRows() != image->getNumRows() || adjacency_.getNumCols() != image->getNumCols()) {
            throw std::invalid_argument("Residual-tree engine adjacency domain differs from the image.");
        }
        if constexpr (RequiresSaturationCertification) {
            if (modeConfiguration_.infinityPixel < 0 || modeConfiguration_.infinityPixel >= image->getSize()) {
                throw std::invalid_argument("Residual-tree engine infinity pixel lies outside the image.");
            }
        }

        constexpr NodeId domainTraversalSeed = NodeId{0};
        std::vector<std::uint8_t> visited(static_cast<std::size_t>(image->getSize()), std::uint8_t{0});
        std::vector<NodeId> frontier;
        frontier.reserve(static_cast<std::size_t>(image->getSize()));
        frontier.push_back(domainTraversalSeed);
        visited[static_cast<std::size_t>(domainTraversalSeed)] = 1;
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
     * @param tree Component-tree seed to validate.
     * @param image Image whose domain the seed must match.
     * @param kind Required component-tree polarity.
     * @param label Human-readable seed label used in diagnostics.
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
     * @brief Checks whether a node is alive and remains a regional extremum.
     * @param tree Component tree containing the node.
     * @param nodeId Node identifier to inspect.
     * @return `true` when `nodeId` is a live non-root leaf with a non-empty proper part.
     */
    [[nodiscard]] static bool isCurrentCandidateNode(const tree_t& tree, NodeId nodeId) {
        const MorphologicalTree& topology = tree.topology();
        return nodeId >= 0 && topology.isNode(nodeId) && topology.isAlive(nodeId) && nodeId != topology.getRoot() && topology.isLeaf(nodeId) &&
               topology.getNodeParent(nodeId) != InvalidNode && topology.getNodeParent(nodeId) != nodeId && topology.getNumProperParts(nodeId) > 0;
    }

    /**
     * @brief Refreshes one agenda slot exclusively from current partition metadata.
     * @param state Mutable construction state owning the agenda and partition.
     * @param tree Component tree containing the candidate node.
     * @param polarity Polarity of `tree`.
     * @param nodeId Node whose agenda entry is refreshed.
     */
    static void updateAgendaCandidate(ConstructionState& state, const tree_t& tree, Polarity polarity, NodeId nodeId) {
        if (isCurrentCandidateNode(tree, nodeId)) {
            state.agenda.upsert(polarity, nodeId, state.flatZonePartition.describeNode(tree, nodeId));
            return;
        }
        state.agenda.reject(polarity, nodeId);
    }

    /**
     * @brief Indexes every current regional extremum of one component tree.
     * @param state Mutable construction state owning the agenda.
     * @param tree Component tree whose extrema are indexed.
     * @param polarity Polarity of `tree`.
     */
    static void indexTreeCandidates(ConstructionState& state, const tree_t& tree, Polarity polarity) {
        for (NodeId nodeId : tree.topology().getAliveNodeIds()) {
            updateAgendaCandidate(state, tree, polarity, nodeId);
        }
    }

    /**
     * @brief Initializes all mutable state for one synchronized construction.
     * @param image Image shared by the component-tree seeds.
     * @param minTree Min-tree seed transferred into mutable state.
     * @param maxTree Max-tree seed transferred into mutable state.
     * @return Fully initialized synchronized construction state.
     */
    [[nodiscard]] ConstructionState initializeConstruction(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree) const {
        requireSeedTree(maxTree, image, MorphologicalTreeKind::MAX_TREE, "max-tree");
        requireSeedTree(minTree, image, MorphologicalTreeKind::MIN_TREE, "min-tree");
        const std::size_t maxTreeSlots = static_cast<std::size_t>(maxTree.topology().getNumInternalNodeSlots());
        const std::size_t minTreeSlots = static_cast<std::size_t>(minTree.topology().getNumInternalNodeSlots());
        const std::size_t maxNodeSlots = std::max(maxTreeSlots, minTreeSlots);

        FlatZonePartitionType flatZonePartition;
        if constexpr (RequiresSaturationCertification) {
            flatZonePartition.initialize(image, adjacency_, modeConfiguration_.infinityPixel);
        } else {
            flatZonePartition.initialize(image, adjacency_);
        }
        const auto initialRegionByPixel = flatZonePartition.representativesByPixel();
        ConstructionState state(tiePolicy_, static_cast<std::size_t>(image->getSize()), maxNodeSlots, std::move(flatZonePartition), initialRegionByPixel);
        state.maxTree = std::make_unique<tree_t>(std::move(maxTree));
        state.minTree = std::make_unique<tree_t>(std::move(minTree));
        state.adjustment = std::make_unique<Adjustment>(state.minTree.get(), state.maxTree.get(), adjacency_);
        state.agenda.reset(maxTreeSlots, minTreeSlots);
        if constexpr (RequiresSaturationCertification) {
            state.saturatedEligibility.emplace(static_cast<std::size_t>(image->getSize()), adjacency_, modeConfiguration_.infinityPixel,
                                                modeConfiguration_.lcaPolicy, modeConfiguration_.fallbackPolicy, *state.minTree, *state.maxTree);
        }
        indexTreeCandidates(state, *state.maxTree, Polarity::Max);
        indexTreeCandidates(state, *state.minTree, Polarity::Min);
        return state;
    }

    /**
     * @brief Checks whether a cached candidate still matches current partition metadata.
     * @param state Mutable construction state owning the partition.
     * @param tree Component tree containing the candidate.
     * @param candidate Cached agenda candidate.
     * @return `true` when the candidate remains current and its ordering metadata is unchanged.
     */
    [[nodiscard]] static bool isCandidateCurrent(ConstructionState& state, const tree_t& tree, const Candidate& candidate) {
        if (!isCurrentCandidateNode(tree, candidate.nodeId)) {
            return false;
        }
        const CandidateDescriptor descriptor = state.flatZonePartition.describeNode(tree, candidate.nodeId);
        return descriptor.area == candidate.area && descriptor.stableSpatialKey == candidate.stableSpatialKey;
    }

    /**
     * @brief Prepares one candidate and applies mode-specific eligibility.
     * @param state Mutable construction state and candidate scratch storage.
     * @param candidate Current candidate selected from the agenda.
     * @param statistics Build diagnostics updated by saturated certification.
     * @return `true` when the candidate may be contracted in the selected mode.
     */
    [[nodiscard]] static bool prepareAndTestEligibility(ConstructionState& state, const Candidate& candidate,
                                                        ResidualTreeBuildStatistics& statistics) {
        const tree_t& primal = candidate.polarity == Polarity::Max ? *state.maxTree : *state.minTree;
        const tree_t& dual = candidate.polarity == Polarity::Max ? *state.minTree : *state.maxTree;
        const bool containsExteriorSeed =
            prepareResidualTreeCandidate(state.flatZonePartition, state.candidateContext, candidate.nodeId, primal, dual);
        if constexpr (RequiresSaturationCertification) {
            if (!state.saturatedEligibility.has_value()) {
                throw std::logic_error("Saturated residual construction has no eligibility evaluator.");
            }
            return state.saturatedEligibility->isEligible(primal, dual, candidate.nodeId, candidate.polarity == Polarity::Max,
                                                          containsExteriorSeed, state.candidateContext, statistics);
        }
        return true;
    }

    /**
     * @brief Commits one eligible extremum to both component trees and the residual event stream.
     * @param state Mutable synchronized construction state.
     * @param candidate Prepared and eligible agenda candidate.
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
        state.flatZonePartition.collectAdjacentRepresentativesAtLevel(
            state.candidateContext.flatZoneRepresentative, targetLevel, state.candidateContext.boundaryPixels,
            state.candidateContext.flatZoneRepresentativeMarks, state.candidateContext.flatZoneMergeRepresentatives);

        if (candidate.polarity == Polarity::Max) {
            state.adjustment->pruneMaxLeafAndUpdateMinTree(candidate.nodeId, support, dualExtremalOwner, boundaryOwners,
                                                           state.candidateContext.wholeSupportOwner);
        } else {
            state.adjustment->pruneMinLeafAndUpdateMaxTree(candidate.nodeId, support, dualExtremalOwner, boundaryOwners,
                                                           state.candidateContext.wholeSupportOwner);
        }

        const NodeId selectedRepresentative = state.candidateContext.flatZoneRepresentative;
        const NodeId mergedRepresentative = state.flatZonePartition.mergeFlatZonesAtLevel(
            selectedRepresentative, state.candidateContext.flatZoneMergeRepresentatives, targetLevel);
        state.absorbedFlatZones.clear();
        if (selectedRepresentative != mergedRepresentative) {
            state.absorbedFlatZones.push_back(static_cast<RegionId>(selectedRepresentative));
        }
        for (NodeId representative : state.candidateContext.flatZoneMergeRepresentatives) {
            if (representative != mergedRepresentative) {
                state.absorbedFlatZones.push_back(static_cast<RegionId>(representative));
            }
        }
        state.residualAssembler.consume(static_cast<RegionId>(mergedRepresentative), state.absorbedFlatZones);
        state.candidateContext.supportPixels = {};

        if constexpr (RequiresSaturationCertification) {
            state.saturatedEligibility->noteDualMutation(candidate.polarity == Polarity::Max,
                                                          state.adjustment->getLastTopologyChangedNodes());
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
     * @brief Emits all residual events while synchronizing both component trees.
     * @param state Mutable synchronized construction state.
     * @param statistics Build diagnostics updated for emitted and rejected extrema.
     * @return Common terminal altitude of the contracted min-tree and max-tree roots.
     */
    [[nodiscard]] static altitude_t buildEvents(ConstructionState& state, ResidualTreeBuildStatistics& statistics) {
        int previousArea = 0;
        while (const auto selected = state.agenda.select()) {
            const Candidate candidate = *selected;
            const tree_t& primal = candidate.polarity == Polarity::Max ? *state.maxTree : *state.minTree;
            if (!isCandidateCurrent(state, primal, candidate)) {
                updateAgendaCandidate(state, primal, candidate.polarity, candidate.nodeId);
                continue;
            }
            if (!prepareAndTestEligibility(state, candidate, statistics)) {
                state.agenda.reject(candidate.polarity, candidate.nodeId);
                ++statistics.rejectedExtrema;
                continue;
            }
            if (candidate.area < previousArea) {
                throw std::runtime_error("Min/max residual event areas are not nondecreasing.");
            }
            previousArea = candidate.area;

            const altitude_t eventAltitude = primal.getAltitude(candidate.nodeId);
            static_cast<void>(state.residualAssembler.emitEvent(static_cast<RegionId>(state.candidateContext.flatZoneRepresentative), eventAltitude));
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

    /** @brief Requires a successfully built result to be available. */
    void requireBuilt() const {
        if (!built_) {
            throw std::logic_error("Residual-tree engine accessors require a completed build.");
        }
    }

    /**
     * @brief Builds and materializes the residual tree from synchronized seeds.
     * @param image Image shared by both component-tree seeds.
     * @param minTree Min-tree seed transferred into the construction.
     * @param maxTree Max-tree seed transferred into the construction.
     * @param localStatistics Diagnostics accumulated during this build.
     */
    void buildFromSeeds(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree, ResidualTreeBuildStatistics& localStatistics) {
        ConstructionState state = initializeConstruction(image, std::move(minTree), std::move(maxTree));
        const altitude_t terminalAltitude = buildEvents(state, localStatistics);
        if (state.residualAssembler.numEvents() != localStatistics.residualEvents) {
            throw std::runtime_error("Min/max residual incremental assembly changed the residual event count.");
        }

        constexpr NodeId terminalAnchorPixel = NodeId{0};
        const NodeId terminalRepresentative = state.flatZonePartition.representativeOf(terminalAnchorPixel);
        auto assemblerOutput = std::move(state.residualAssembler).finalize(static_cast<RegionId>(terminalRepresentative), terminalAltitude);
        auto result = materializeResidualTree(image, adjacency_, std::move(assemblerOutput));
        rows_ = result.rows;
        cols_ = result.cols;
        root_ = result.root;
        nodeParent_ = std::move(result.nodeParent);
        properPartOwner_ = std::move(result.properPartOwner);
        altitude_ = std::move(result.altitude);
        topologyProof_ = std::move(result.topologyProof);
        statistics_ = localStatistics;
        built_ = true;
    }

  public:
    /**
     * @brief Configures saturated construction with an explicit exterior seed.
     * @param adjacency Symmetric adjacency over the image domain.
     * @param infinityPixel Exterior seed pixel excluded from residual candidates.
     * @param options Saturated construction policies.
     */
    explicit MinMaxResidualTreeEngine(RegularGridAdjacency2D adjacency, NodeId infinityPixel, SaturatedResidualTreeOptions options)
        requires(RequiresSaturationCertification)
        : adjacency_(std::move(adjacency)), tiePolicy_(options.tiePolicy),
          modeConfiguration_{infinityPixel, options.lcaPolicy, options.fallbackPolicy} {}

    /**
     * @brief Configures unrestricted construction without exterior state.
     * @param adjacency Symmetric adjacency over the image domain.
     * @param options Unrestricted construction policies.
     */
    explicit MinMaxResidualTreeEngine(RegularGridAdjacency2D adjacency, UnrestrictedResidualTreeOptions options)
        requires(!RequiresSaturationCertification)
        : adjacency_(std::move(adjacency)), tiePolicy_(options.tiePolicy), modeConfiguration_{} {}

    /**
     * @brief Consumes caller-built synchronized min-tree and max-tree seeds.
     * @param image Image shared by both component-tree seeds.
     * @param minTree Min-tree seed transferred into the construction.
     * @param maxTree Max-tree seed transferred into the construction.
     */
    void build(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree) {
        clearResult();
        ResidualTreeBuildStatistics localStatistics;
        requireInputs(image);
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(image, "Residual-tree engine input image");
        buildFromSeeds(image, std::move(minTree), std::move(maxTree), localStatistics);
    }

    /**
     * @brief Returns the number of rows in the completed hierarchy domain.
     * @return Number of rows.
     */
    [[nodiscard]] int getRows() const {
        requireBuilt();
        return rows_;
    }

    /**
     * @brief Returns the number of columns in the completed hierarchy domain.
     * @return Number of columns.
     */
    [[nodiscard]] int getCols() const {
        requireBuilt();
        return cols_;
    }

    /**
     * @brief Returns the root of the completed residual hierarchy.
     * @return Residual root node identifier.
     */
    [[nodiscard]] NodeId getRoot() const {
        requireBuilt();
        return root_;
    }

    /**
     * @brief Returns the configured saturated exterior seed.
     * @return Row-major exterior seed pixel.
     */
    [[nodiscard]] NodeId getInfinityPixel() const noexcept
        requires(RequiresSaturationCertification)
    {
        return modeConfiguration_.infinityPixel;
    }

    /**
     * @brief Returns the adjacency shared by the construction.
     * @return Configured symmetric grid adjacency.
     */
    [[nodiscard]] const RegularGridAdjacency2D& getAdjacency() const noexcept { return adjacency_; }

    /**
     * @brief Returns the deterministic equal-area ordering policy.
     * @return Configured tie policy.
     */
    [[nodiscard]] SdrtTiePolicy getTiePolicy() const noexcept { return tiePolicy_; }

    /**
     * @brief Returns the saturated dynamic LCA policy.
     * @return Configured LCA policy.
     */
    [[nodiscard]] SaturatedMinMaxLcaPolicy getLcaPolicy() const noexcept
        requires(RequiresSaturationCertification)
    {
        return modeConfiguration_.lcaPolicy;
    }

    /**
     * @brief Returns the saturated exact-complement fallback policy.
     * @return Configured fallback policy.
     */
    [[nodiscard]] SaturatedMinMaxFallbackPolicy getFallbackPolicy() const noexcept
        requires(RequiresSaturationCertification)
    {
        return modeConfiguration_.fallbackPolicy;
    }

    /**
     * @brief Returns the completed residual parent buffer.
     * @return Parent identifier indexed by residual node.
     */
    [[nodiscard]] std::span<const NodeId> getNodeParent() const {
        requireBuilt();
        return nodeParent_;
    }

    /**
     * @brief Returns the completed direct proper-part owner buffer.
     * @return Residual owner indexed by row-major pixel.
     */
    [[nodiscard]] std::span<const NodeId> getProperPartOwner() const {
        requireBuilt();
        return properPartOwner_;
    }

    /**
     * @brief Returns the completed residual altitude buffer.
     * @return Altitude indexed by residual node.
     */
    [[nodiscard]] std::span<const altitude_t> getAltitude() const {
        requireBuilt();
        return altitude_;
    }

    /**
     * @brief Returns diagnostics from the completed construction.
     * @return Build statistics retained by the engine.
     */
    [[nodiscard]] const ResidualTreeBuildStatistics& getStatistics() const {
        requireBuilt();
        return statistics_;
    }

    /**
     * @brief Transfers the completed result into validated native storage.
     * @param semantics Semantic metadata attached to the resulting hierarchy.
     * @return Validated native hierarchy that owns the completed buffers.
     */
    [[nodiscard]] mmcfilters::detail::ValidatedNativeHierarchy<altitude_t> takeValidatedHierarchy(HierarchySemantics semantics) && {
        requireBuilt();
        auto hierarchy = mmcfilters::detail::makeValidatedNativeHierarchy<altitude_t>(
            std::move(nodeParent_), std::move(properPartOwner_), std::move(altitude_), root_, GridDomain2D{rows_, cols_}, std::move(semantics),
            std::move(topologyProof_));
        built_ = false;
        rows_ = 0;
        cols_ = 0;
        root_ = InvalidNode;
        return hierarchy;
    }
};

} // namespace mmcfilters::sdrt::detail
