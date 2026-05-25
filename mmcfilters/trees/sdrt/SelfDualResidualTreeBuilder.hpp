#pragma once

#include "detail/LeafPrunePriorityAgenda.hpp"
#include "../adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters::sdrt {

/**
 * @brief Stateful builder that materializes a self-dual residual tree as native buffers.
 *
 * The construction starts from caller-provided max-tree and min-tree seeds for
 * an image whose pixel type defines the tree altitude type. At each step it selects the smallest currently prunable leaf
 * from either polarity, records the residual support as a persistent per-pixel
 * chain event, and then prunes that leaf while incrementally updating the dual
 * tree. When both component trees collapse to the same constant root, the
 * recorded chains are materialized into native topology and altitude buffers
 * owned by the builder instance.
 *
 * Public contract:
 *
 * - `build()` consumes the caller-provided weighted min/max trees by move;
 * - the input image domain must match both seed trees;
 * - output buffers become available only after a successful `build()`;
 * - accessors return read-only views into storage owned by the builder;
 * - the builder may be reused, and a new successful `build()` replaces the
 *   previously materialized buffers.
 *
 * Ownership note:
 *
 * SDRT construction intentionally owns mutable `WeightedMorphologicalTree<T>`
 * seeds during a build. `WeightedTreeView` is not a replacement here because
 * each leaf prune updates the component-tree pair before the next residual
 * event is selected.
 *
 * Performance intent:
 *
 * The builder avoids rebuilding the candidate frontier after every pruning
 * step. `LeafPrunePriorityAgenda<T>` is updated only around nodes reported by the
 * dual-tree leaf adjuster, while persistent chains share suffixes between
 * pixels with identical residual histories. This keeps the construction close
 * to the cost of the underlying component-tree edits plus logarithmic candidate
 * agenda operations.
 */
template<AltitudeValue T>
class SelfDualResidualTreeBuilder {
private:
    using altitude_t = T;
    using image_ptr_t = ImagePtr<altitude_t>;
    using tree_t = WeightedMorphologicalTree<altitude_t>;

    // Leaf-only adjuster that keeps the min-tree/max-tree pair synchronized
    // after pruning a single candidate leaf from either polarity.
    using SdrtAdjustment = ::mmcfilters::adjust::DualMinMaxTreeIncrementalFilterLeaf<altitude_t>;

    // Local aliases keep the builder readable while the underlying helper
    // types remain scoped as implementation details.
    using CandidateLeaf = detail::CandidateLeaf;
    using ResidualPolarity = detail::ResidualPolarity;

    /**
     * @brief One node in a persistent per-pixel residual stack.
     *
     * `eventId` identifies the pruning event stored at this stack level.
     * `prev` points to the previous stack node, or `-1` for an empty suffix.
     * Many pixels may point to the same stack node, which is how suffix sharing
     * is preserved without copying full event chains per pixel.
     */
    struct PersistentStackNode {
        int eventId = -1;
        int prev = -1;
    };

    /**
     * @brief Persistent-chain store accumulated during one SDRT build.
     *
     * `stackHeadByPixel[p]` points to the latest residual event affecting
     * pixel/proper-part `p`. `persistentStackNodes` stores shared chain nodes,
     * and `eventValuations` maps event ids to the altitude assigned to the SDRT
     * node created for that event.
     */
    struct PersistentChains {
        std::vector<int> stackHeadByPixel;
        std::vector<PersistentStackNode> persistentStackNodes;
        std::vector<altitude_t> eventValuations;
    };

    double radius_ = 1.5;
    bool built_ = false;
    int rows_ = 0;
    int cols_ = 0;
    NodeId root_ = InvalidNode;
    std::vector<NodeId> nodeParent_;
    std::vector<NodeId> properPartOwner_;
    std::vector<altitude_t> altitude_;

    /**
     * Mutable state for one SDRT build.
     *
     * The two weighted trees are edited throughout construction. The agenda
     * stores only current leaf candidates and is updated from the nodes touched
     * by the last dual-tree adjustment. chains stores the persistent residual
     * history that will later be converted into native MAF arrays.
     *
     * cachedNewHeadByPrevHead and cachedNewHeadStampByPrevHead are per-event
     * memoization arrays. During one residual event, many pixels may share the
     * same previous chain head. The cache guarantees that a shared previous head
     * creates only one new PersistentStackNode for that event, preserving suffix
     * sharing in the persistent chains.
     */
    struct ConstructionState {
        AdjacencyRelation adjacency;
        std::unique_ptr<tree_t> maxTree;
        std::unique_ptr<tree_t> minTree;
        std::unique_ptr<SdrtAdjustment> adjustment;
        detail::LeafPrunePriorityAgenda<altitude_t> agenda;
        PersistentChains chains;
        std::vector<int> cachedNewHeadByPrevHead;
        std::vector<std::uint32_t> cachedNewHeadStampByPrevHead;
        std::uint32_t currentEventCacheStamp = 1;

        /**
         * @brief Creates the per-build adjacency context used by dual updates.
         */
        ConstructionState(int rows, int cols, double radius)
            : adjacency(rows, cols, radius) {}
    };

    // Validates the minimal image-domain assumptions required by both
    // component-tree construction and persistent-chain materialization.
    static void requireImage(const image_ptr_t& image) {
        if (!image) {
            throw std::invalid_argument("SelfDualResidualTreeBuilder requires a non-null image.");
        }
        if (image->getNumRows() <= 0 || image->getNumCols() <= 0 || image->getSize() <= 0) {
            throw std::invalid_argument("SelfDualResidualTreeBuilder requires a non-empty 2D image.");
        }
    }

    // Rechecks a candidate selected from the agenda against the current tree.
    // The agenda is maintained incrementally, so this guard protects the build
    // loop from stale entries left by structural edits that changed leaf status
    // or support area before the candidate reached the front of the queue.
    static bool isCandidateStillPrunable(const tree_t& tree, const CandidateLeaf& candidate) {
        const MorphologicalTree& topology = tree.topology();
        if (candidate.nodeId == InvalidNode ||
            candidate.nodeId == topology.getRoot() ||
            !topology.isNode(candidate.nodeId) ||
            !topology.isAlive(candidate.nodeId) ||
            !topology.isLeaf(candidate.nodeId)) {
            return false;
        }
        const NodeId parentId = topology.getNodeParent(candidate.nodeId);
        if (parentId == InvalidNode || parentId == candidate.nodeId) {
            throw std::runtime_error("SDRT candidate violates the component-tree parent invariant.");
        }
        return topology.getNumProperParts(candidate.nodeId) == candidate.area;
    }

    // Tests whether a node would be admitted by the agenda if synchronized now.
    // This is used after dual-tree updates to distinguish affected nodes that
    // should enter the agenda from affected nodes that only need removal of a
    // previous candidate entry.
    static bool isCurrentCandidateNode(const tree_t& tree, NodeId nodeId) {
        const MorphologicalTree& topology = tree.topology();
        if (nodeId == InvalidNode ||
            nodeId == topology.getRoot() ||
            !topology.isNode(nodeId) ||
            !topology.isAlive(nodeId) ||
            !topology.isLeaf(nodeId)) {
            return false;
        }
        const NodeId parentId = topology.getNodeParent(nodeId);
        return parentId != InvalidNode &&
               parentId != nodeId &&
               topology.getNumProperParts(nodeId) > 0;
    }

    // Ensures the per-event chain-head memoization arrays can address
    // prevHead + 1. The +1 encoding reserves index 0 for the sentinel prevHead
    // value -1, which represents an empty previous chain.
    static void ensureChainCacheCapacity(ConstructionState& state, std::size_t requiredSize) {
        if (state.cachedNewHeadByPrevHead.size() >= requiredSize) {
            return;
        }
        std::size_t newSize = state.cachedNewHeadByPrevHead.empty() ? 2 : state.cachedNewHeadByPrevHead.size();
        while (newSize < requiredSize) {
            newSize *= 2;
        }
        state.cachedNewHeadByPrevHead.resize(newSize, -1);
        state.cachedNewHeadStampByPrevHead.resize(newSize, 0);
    }

    // Returns the persistent stack node that represents "eventId on top of
    // prevHead" for the current residual event. Calls with the same prevHead
    // inside one event return the same node id, so pixels with identical chain
    // history share the newly created suffix node.
    static int makeNewHead(ConstructionState& state, int eventId, int prevHead) {
        const auto cacheIndex = static_cast<std::size_t>(prevHead + 1);
        ensureChainCacheCapacity(state, cacheIndex + 1);
        if (state.cachedNewHeadStampByPrevHead[cacheIndex] == state.currentEventCacheStamp) {
            return state.cachedNewHeadByPrevHead[cacheIndex];
        }

        const int newHead = static_cast<int>(state.chains.persistentStackNodes.size());
        state.chains.persistentStackNodes.push_back(PersistentStackNode{eventId, prevHead});
        state.cachedNewHeadByPrevHead[cacheIndex] = newHead;
        state.cachedNewHeadStampByPrevHead[cacheIndex] = state.currentEventCacheStamp;
        return newHead;
    }

    /**
     * Records one pruning event in the persistent-chain representation.
     *
     * valuation is the altitude of the pruned component. support is the proper
     * part set owned by that component at the moment it is pruned. Every support
     * pixel receives a new chain head for this event, while pixels outside the
     * support keep their previous head untouched.
     *
     * The event stamp resets the per-event cache without clearing the backing
     * arrays on every iteration. If the uint32 stamp wraps, all stamps are
     * cleared and numbering restarts from 1.
     */
    template<class PixelRangeT>
    static void recordEvent(ConstructionState& state, altitude_t valuation, PixelRangeT&& support) {
        ensureChainCacheCapacity(state, state.chains.persistentStackNodes.size() + 1);
        ++state.currentEventCacheStamp;
        if (state.currentEventCacheStamp == 0) {
            std::fill(state.cachedNewHeadStampByPrevHead.begin(), state.cachedNewHeadStampByPrevHead.end(), 0);
            state.currentEventCacheStamp = 1;
        }

        const int eventId = static_cast<int>(state.chains.eventValuations.size());
        int countedSupport = 0;
        for (NodeId pixelId : support) {
            if (pixelId < 0 ||
                pixelId >= static_cast<NodeId>(state.chains.stackHeadByPixel.size())) {
                throw std::out_of_range("SDRT builder received an invalid support pixel id.");
            }
            const auto index = static_cast<std::size_t>(pixelId);
            const int prevHead = state.chains.stackHeadByPixel[index];
            state.chains.stackHeadByPixel[index] = makeNewHead(state, eventId, prevHead);
            ++countedSupport;
        }
        if (countedSupport <= 0) {
            throw std::runtime_error("SDRT builder cannot record an empty residual support.");
        }
        state.chains.eventValuations.push_back(valuation);
    }

    static void requireSeedTree(
        const tree_t& tree,
        const image_ptr_t& image,
        MorphologicalTreeKind expectedTreeType,
        const char* label) {
        const MorphologicalTree& topology = tree.topology();
        if (topology.getTreeType() != expectedTreeType) {
            throw std::invalid_argument(std::string("SDRT ") + label + " seed has an unexpected tree type.");
        }
        if (topology.getNumRowsOfImage() != image->getNumRows() ||
            topology.getNumColsOfImage() != image->getNumCols() ||
            topology.getNumTotalProperParts() != image->getSize()) {
            throw std::invalid_argument(std::string("SDRT ") + label + " seed domain must match the input image.");
        }
        tree.validateAltitudeBufferShape();
    }

    // Takes ownership of the initial max-tree/min-tree pair, then creates the
    // leaf adjuster, persistent-chain store and first complete candidate
    // agenda. After this point the build loop maintains the agenda
    // incrementally.
    static ConstructionState initializeConstruction(
        const image_ptr_t& image,
        double radius,
        tree_t&& minTree,
        tree_t&& maxTree) {
        requireSeedTree(maxTree, image, MorphologicalTreeKind::MAX_TREE, "max-tree");
        requireSeedTree(minTree, image, MorphologicalTreeKind::MIN_TREE, "min-tree");

        ConstructionState state(image->getNumRows(), image->getNumCols(), radius);
        state.maxTree = std::make_unique<tree_t>(std::move(maxTree));
        state.minTree = std::make_unique<tree_t>(std::move(minTree));
        state.adjustment =
            std::make_unique<SdrtAdjustment>(state.minTree.get(), state.maxTree.get(), state.adjacency);

        const auto numPixels = static_cast<std::size_t>(image->getSize());
        state.chains.stackHeadByPixel.assign(numPixels, -1);
        state.chains.persistentStackNodes.clear();
        state.chains.eventValuations.clear();
        state.cachedNewHeadByPrevHead.assign(std::max<std::size_t>(numPixels + 1, 2), -1);
        state.cachedNewHeadStampByPrevHead.assign(std::max<std::size_t>(numPixels + 1, 2), 0);
        state.currentEventCacheStamp = 1;
        state.agenda.initialize(*state.maxTree, *state.minTree);
        return state;
    }

    // Returns the final root altitude after both component trees have collapsed
    // to one constant component. The max and min roots must agree; disagreement
    // indicates the dual update sequence failed to preserve a valid bipolar
    // residual state.
    static altitude_t terminalValueForCompletedState(const ConstructionState& state) {
        if (state.maxTree->topology().getNumNodes() != 1 ||
            state.minTree->topology().getNumNodes() != 1) {
            throw std::runtime_error("SDRT builder exhausted candidates before both component trees became constant.");
        }

        const NodeId maxRoot = state.maxTree->topology().getRoot();
        const NodeId minRoot = state.minTree->topology().getRoot();
        if (maxRoot == InvalidNode || minRoot == InvalidNode ||
            !state.maxTree->topology().isAlive(maxRoot) ||
            !state.minTree->topology().isAlive(minRoot)) {
            throw std::runtime_error("SDRT builder produced an invalid terminal root.");
        }

        const altitude_t maxAltitude = state.maxTree->getAltitude(maxRoot);
        const altitude_t minAltitude = state.minTree->getAltitude(minRoot);
        if (maxAltitude != minAltitude) {
            throw std::runtime_error("SDRT builder produced inconsistent terminal max/min root altitudes.");
        }
        return maxAltitude;
    }

    // Applies the selected pruning step and updates the agenda only around the
    // nodes that may have changed. The pruned tree is called primal here; the
    // opposite tree is updated by the adjustment object and is called dual.
    //
    // The leaf adjuster reports affected candidate nodes on the dual side.
    // Those nodes are revalidated if they were already active in the agenda or
    // if they are currently valid leaf candidates. On the primal side, only the
    // pruned node and its former parent can change candidate status.
    static void commitCandidate(ConstructionState& state, const CandidateLeaf& candidate) {
        tree_t& primalTree = candidate.polarity == ResidualPolarity::Max ? *state.maxTree : *state.minTree;
        tree_t& dualTree = candidate.polarity == ResidualPolarity::Max ? *state.minTree : *state.maxTree;
        const ResidualPolarity dualPolarity = candidate.polarity == ResidualPolarity::Max ? ResidualPolarity::Min : ResidualPolarity::Max;
        const NodeId parentId = primalTree.topology().getNodeParent(candidate.nodeId);

        if (candidate.polarity == ResidualPolarity::Max) {
            state.adjustment->pruneMaxLeafAndUpdateMinTree(candidate.nodeId);
        } else {
            state.adjustment->pruneMinLeafAndUpdateMaxTree(candidate.nodeId);
        }

        const auto& candidateNodes = state.adjustment->getLastCandidateNodes();
        for (NodeId nodeId : candidateNodes) {
            if (state.agenda.contains(dualPolarity, nodeId) || isCurrentCandidateNode(dualTree, nodeId)) {
                state.agenda.update(dualTree, dualPolarity, nodeId);
            }
        }

        state.agenda.update(primalTree, candidate.polarity, candidate.nodeId);
        if (isCurrentCandidateNode(primalTree, parentId)) {
            state.agenda.update(primalTree, candidate.polarity, parentId);
        }
    }

    // Main construction loop. It repeatedly selects the smallest candidate,
    // records its residual event using the candidate tree support, and commits
    // the pruning. Stale agenda entries are repaired in place and skipped.
    //
    // The returned altitude is the common terminal value used for the SDRT root
    // when the persistent chains are materialized.
    static altitude_t buildPersistentChains(ConstructionState& state) {
        while (true) {
            auto candidate = state.agenda.select();
            if (!candidate.has_value()) {
                break;
            }

            const tree_t& supportTree = candidate->polarity == ResidualPolarity::Max ? *state.maxTree : *state.minTree;
            if (!isCandidateStillPrunable(supportTree, *candidate)) {
                state.agenda.update(supportTree, candidate->polarity, candidate->nodeId);
                continue;
            }

            recordEvent(state, supportTree.getAltitude(candidate->nodeId), supportTree.topology().getProperParts(candidate->nodeId));
            commitCandidate(state, *candidate);
        }

        return terminalValueForCompletedState(state);
    }

    /**
     * @brief Ensures accessors are used only after successful materialization.
     */
    void requireBuilt() const {
        if (!built_) {
            throw std::logic_error("SelfDualResidualTreeBuilder accessors require a completed build.");
        }
    }

    /**
     * Converts the persistent per-pixel residual chains recorded during
     * construction into the native MAF buffers owned by this builder.
     *
     * The materialized topology uses one node per pruning event plus a terminal
     * root. Event id E maps to node id E + 1; node id 0 is the terminal root.
     */
    void materializeChains(int rows,
                           int cols,
                           altitude_t terminalValue,
                           const PersistentChains& chains) {
        if (rows <= 0 || cols <= 0) {
            throw std::invalid_argument("SDRT builder requires a non-empty image domain for chain materialization.");
        }
        const auto& heads = chains.stackHeadByPixel;
        const auto& nodes = chains.persistentStackNodes;
        const auto& eventValuations = chains.eventValuations;
        if (heads.size() != static_cast<std::size_t>(rows * cols)) {
            throw std::invalid_argument("SDRT builder pixel-chain domain must match rows * cols.");
        }

        rows_ = rows;
        cols_ = cols;
        root_ = 0;
        nodeParent_.assign(eventValuations.size() + 1, InvalidNode);
        nodeParent_[0] = 0;
        properPartOwner_.assign(heads.size(), root_);
        altitude_.assign(eventValuations.size() + 1, altitude_t{});
        altitude_[0] = terminalValue;
        for (std::size_t eventId = 0; eventId < eventValuations.size(); ++eventId) {
            altitude_[eventId + 1] = eventValuations[eventId];
        }

        std::vector<NodeId> parentByEvent(eventValuations.size(), InvalidNode);
        std::vector<int> ownerEventByHead(nodes.size(), -2);
        std::vector<int> headStack;
        headStack.reserve(32);

        for (int head : heads) {
            if (head < 0) {
                continue;
            }
            if (static_cast<std::size_t>(head) >= nodes.size()) {
                throw std::out_of_range("SDRT builder received an invalid residual-chain head.");
            }
            if (ownerEventByHead[static_cast<std::size_t>(head)] != -2) {
                continue;
            }

            int cursor = head;
            while (cursor >= 0) {
                if (static_cast<std::size_t>(cursor) >= nodes.size()) {
                    throw std::out_of_range("SDRT builder found an invalid residual-chain predecessor.");
                }
                if (ownerEventByHead[static_cast<std::size_t>(cursor)] != -2) {
                    break;
                }
                headStack.push_back(cursor);
                cursor = nodes[static_cast<std::size_t>(cursor)].prev;
            }

            int ownerEvent = (cursor >= 0) ? ownerEventByHead[static_cast<std::size_t>(cursor)] : -1;
            while (!headStack.empty()) {
                const int currentHead = headStack.back();
                headStack.pop_back();
                const auto& currentNode = nodes[static_cast<std::size_t>(currentHead)];
                if (currentNode.eventId < 0 ||
                    currentNode.eventId >= static_cast<int>(eventValuations.size())) {
                    throw std::out_of_range("SDRT builder found an invalid event id.");
                }

                ownerEvent = (ownerEvent >= 0) ? ownerEvent : currentNode.eventId;
                ownerEventByHead[static_cast<std::size_t>(currentHead)] = ownerEvent;

                if (currentNode.prev >= 0) {
                    if (static_cast<std::size_t>(currentNode.prev) >= nodes.size()) {
                        throw std::out_of_range("SDRT builder found an invalid residual-chain predecessor.");
                    }
                    const auto childEvent =
                        static_cast<std::size_t>(nodes[static_cast<std::size_t>(currentNode.prev)].eventId);
                    if (childEvent >= eventValuations.size()) {
                        throw std::out_of_range("SDRT builder found an invalid child event id.");
                    }
                    const NodeId candidateParent = static_cast<NodeId>(currentNode.eventId) + 1;
                    if (parentByEvent[childEvent] == InvalidNode) {
                        parentByEvent[childEvent] = candidateParent;
                    } else if (parentByEvent[childEvent] != candidateParent) {
                        throw std::runtime_error("SDRT builder found inconsistent parent candidates.");
                    }
                }
            }
        }

        for (NodeId pixelId = 0; pixelId < static_cast<NodeId>(heads.size()); ++pixelId) {
            const int head = heads[static_cast<std::size_t>(pixelId)];
            if (head < 0) {
                properPartOwner_[static_cast<std::size_t>(pixelId)] = root_;
                continue;
            }
            if (static_cast<std::size_t>(head) >= ownerEventByHead.size()) {
                throw std::out_of_range("SDRT builder received an invalid residual-chain head.");
            }
            const int ownerEvent = ownerEventByHead[static_cast<std::size_t>(head)];
            if (ownerEvent >= static_cast<int>(eventValuations.size())) {
                throw std::out_of_range("SDRT builder found an invalid owner event id.");
            }
            properPartOwner_[static_cast<std::size_t>(pixelId)] =
                ownerEvent < 0 ? root_ : static_cast<NodeId>(ownerEvent) + 1;
        }

        for (std::size_t eventId = 0; eventId < eventValuations.size(); ++eventId) {
            const NodeId nodeId = static_cast<NodeId>(eventId) + 1;
            const NodeId parentId =
                parentByEvent[eventId] == InvalidNode ? root_ : parentByEvent[eventId];
            nodeParent_[static_cast<std::size_t>(nodeId)] = parentId;
        }

        built_ = true;
    }

public:
    /**
     * @brief Creates an SDRT builder with a fixed adjacency radius.
     *
     * @param radius Radius used by the dual-tree adjustment adjacency.
     */
    explicit SelfDualResidualTreeBuilder(double radius = 1.5)
        : radius_(radius) {}

    /**
     * @brief Consumes min/max seed trees and materializes the SDRT buffers.
     *
     * @param image Non-null, non-empty image whose domain matches both seeds.
     * @param minTree Weighted min-tree seed. Moved into the build state.
     * @param maxTree Weighted max-tree seed. Moved into the build state.
     *
     * @throws std::invalid_argument if the image or seed-tree domains are
     * invalid or do not match.
     * @throws std::runtime_error if candidate selection, dual-tree updates, or
     * chain materialization detect an inconsistent state.
     */
    void build(
        const image_ptr_t& image,
        tree_t&& minTree,
        tree_t&& maxTree) {
        built_ = false;
        rows_ = 0;
        cols_ = 0;
        root_ = InvalidNode;
        nodeParent_.clear();
        properPartOwner_.clear();
        altitude_.clear();

        requireImage(image);
        ConstructionState state = initializeConstruction(
            image,
            radius_,
            std::move(minTree),
            std::move(maxTree));
        const altitude_t terminalValue = buildPersistentChains(state);
        materializeChains(
            image->getNumRows(),
            image->getNumCols(),
            terminalValue,
            state.chains);
    }

    /**
     * @brief Returns the row count of the materialized image domain.
     */
    [[nodiscard]] int getRows() const {
        requireBuilt();
        return rows_;
    }

    /**
     * @brief Returns the column count of the materialized image domain.
     */
    [[nodiscard]] int getCols() const {
        requireBuilt();
        return cols_;
    }

    /**
     * @brief Returns the root node id in the materialized internal-node domain.
     */
    [[nodiscard]] NodeId getRoot() const {
        requireBuilt();
        return root_;
    }

    /**
     * @brief Returns parent links indexed by internal `NodeId`.
     *
     * The root is self-parented. The span remains valid until the builder is
     * destroyed or `build()` is called again.
     */
    [[nodiscard]] std::span<const NodeId> getNodeParent() const {
        requireBuilt();
        return nodeParent_;
    }

    /**
     * @brief Returns the smallest owning component for each proper part.
     *
     * The span is indexed by row-major pixel/proper-part id and remains valid
     * until the builder is destroyed or `build()` is called again.
     */
    [[nodiscard]] std::span<const NodeId> getProperPartOwner() const {
        requireBuilt();
        return properPartOwner_;
    }

    /**
     * @brief Returns node altitudes indexed by internal `NodeId`.
     *
     * The span remains valid until the builder is destroyed or `build()` is
     * called again.
     */
    [[nodiscard]] std::span<const altitude_t> getAltitude() const {
        requireBuilt();
        return altitude_;
    }
};


} // namespace mmcfilters::sdrt
