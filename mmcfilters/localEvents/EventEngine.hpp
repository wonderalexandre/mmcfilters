#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../trees/detail/ProperPartEntryNode.hpp"
#include "../utils/Common.hpp"
#include "../utils/Image.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::local_events {

/**
 * @brief Relative integer offset of one sample in a local image window.
 *
 * Offsets are expressed in image coordinates. `rowOffset` moves down for
 * positive values and `colOffset` moves right for positive values.
 */
struct WindowOffset {
    /// Row displacement from the anchor pixel.
    int rowOffset = 0;

    /// Column displacement from the anchor pixel.
    int colOffset = 0;
};

/**
 * @brief Generic event engine for finite binary local computations.
 *
 * @details
 * This engine implements only the structural part of local event
 * computation: entry-node evaluation, event ordering, local-state transitions,
 * and bottom-up bucket aggregation. Attribute-specific headers own the meaning
 * and type of their buckets.
 *
 * The model is based on a finite window sampled around each anchor pixel. For
 * each valid window sample, `entryNode(...)` returns the first node, walking
 * upward from the anchor owner, whose support also contains that sample. Events
 * are then processed from descendants to ancestors, so the binary state only
 * gains visible samples. The policy receives the initial state and every later
 * transition, and can store either raw state counters or an already projected
 * attribute bucket.
 *
 * Returned buffers always use the dense internal node-id space of
 * `MorphologicalTree`; dead slots keep the default value of the policy bucket.
 * This header is public so new local descriptors can reuse the same
 * policy-based event model.
 */
class EventEngine {
  private:
    /**
     * @brief One sample becoming visible at `node`.
     */
    struct EntryEvent {
        /** @brief Stores the node. */
        NodeId node = InvalidNode;
        /** @brief Stores the mask. */
        uint32_t mask = 0;
    };

  public:
    /**
     * @brief Returns the first ancestor of `anchorPixel` that contains `samplePixel`.
     *
     * @details
     * The entry node is the structural point where a binary local sample becomes
     * visible while the anchor support is expanded from a node to its ancestors:
     *
     * - if the anchor owner already contains the sample owner, the entry is the
     *   anchor owner;
     * - if the sample owner is an ancestor of the anchor owner, the entry is the
     *   sample owner;
     * - otherwise the two owners are incomparable and the event enters at their
     *   lowest common ancestor.
     *
     * Invalid pixels or pixels without an owner return `InvalidNode`.
     *
     * @param tree Morphological tree whose proper-part ownership defines the
     * supports being sampled.
     * @param anchorPixel Row-major image-domain proper part used as the window
     * anchor.
     * @param samplePixel Row-major image-domain proper part being tested.
     * @return Dense internal `NodeId` where the sample enters, or `InvalidNode`
     * when either pixel is invalid or unowned.
     */
    static NodeId entryNode(const MorphologicalTree& tree, int anchorPixel, int samplePixel) {
        return detail::properPartEntryNode(tree, anchorPixel, samplePixel);
    }

    /**
     * @brief Returns the entry node for one translated sample around an anchor.
     *
     * @details
     * This overload converts `offset` to an absolute sample pixel in the image
     * domain and returns `InvalidNode` when the translated sample falls outside
     * the image. Out-of-domain samples therefore do not generate local events;
     * attribute policies that need border complements account for them through
     * their own framed-window logic.
     *
     * @param tree Morphological tree whose image domain defines valid samples.
     * @param anchorPixel Row-major image-domain proper part used as the window
     * anchor.
     * @param offset Relative window sample offset.
     * @return Dense internal `NodeId` where the translated sample enters, or
     * `InvalidNode` if the anchor/sample is outside the image or unowned.
     */
    static NodeId entryNode(const MorphologicalTree& tree, int anchorPixel, WindowOffset offset) {
        const int rows = tree.getNumRowsOfGridDomain2D();
        const int cols = tree.getNumColsOfGridDomain2D();
        if (rows <= 0 || cols <= 0) {
            return InvalidNode;
        }
        if (anchorPixel < 0 || anchorPixel >= rows * cols) {
            return InvalidNode;
        }

        const auto [row, col] = ImageUtils::to2D(anchorPixel, cols);
        const int sampleRow = row + offset.rowOffset;
        const int sampleCol = col + offset.colOffset;
        if (sampleRow < 0 || sampleRow >= rows || sampleCol < 0 || sampleCol >= cols) {
            return InvalidNode;
        }

        return entryNode(tree, anchorPixel, ImageUtils::to1D(sampleRow, sampleCol, cols));
    }

    /**
     * @brief Computes local state-change deltas before subtree aggregation.
     *
     * @details
     * The policy must define:
     *
     * - `using Bucket = ...`;
     * - `applyInitial(Bucket&, uint32_t state)`;
     * - `applyTransition(Bucket&, uint32_t oldState, uint32_t newState)`;
     * - `merge(Bucket& parent, const Bucket& child)`.
     *
     * `state` is a bit mask over the samples in `window`. Bit `i` is set after
     * the event for `window[i]` has entered. Because events are sorted
     * bottom-up, transitions are monotone (`newState` is `oldState` plus one or
     * more bits at the same node). The policy is responsible for interpreting a
     * state and adding/removing the previous state's contribution.
     *
     * Returned buckets are the sparse per-node delta representation:
     * aggregating child buckets into parents materializes the final subtree
     * attribute. Attribute implementations can expose this representation when
     * callers need the local-event model itself rather than only final subtree
     * counts.
     *
     * @tparam Policy Policy type defining `Bucket`, `applyInitial`,
     * `applyTransition`, and `merge`.
     * @param tree Tree whose dense internal node-id space indexes the returned
     * bucket vector.
     * @param window Finite local window. At most 32 samples are supported
     * because the active state is encoded in a `uint32_t` mask.
     * @param policy Stateless or externally owned policy object used to update
     * and merge buckets.
     * @return Dense bucket vector sized by `tree.getNumInternalNodeSlots()`.
     * @throws std::invalid_argument If the tree has no image domain, has no
     * live root, or the window contains more than 32 samples.
     */
    template <class Policy>
    static std::vector<typename Policy::Bucket> computeDeltasWithPolicy(const MorphologicalTree& tree, const std::vector<WindowOffset>& window,
                                                                        const Policy& policy) {
        requireValidInput(tree, window);

        std::vector<typename Policy::Bucket> buckets(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        std::vector<EntryEvent> events;
        events.reserve(window.size());

        const int totalPixels = tree.getNumRowsOfGridDomain2D() * tree.getNumColsOfGridDomain2D();
        for (int anchorPixel = 0; anchorPixel < totalPixels; ++anchorPixel) {
            events.clear();
            for (std::size_t i = 0; i < window.size(); ++i) {
                const NodeId entry = entryNode(tree, anchorPixel, window[i]);
                if (entry != InvalidNode) {
                    events.push_back({entry, uint32_t{1} << i});
                }
            }
            if (events.empty()) {
                continue;
            }

            sortEventsBottomUp(tree, events);

            uint32_t state = 0;
            bool hasPreviousState = false;
            for (std::size_t i = 0; i < events.size();) {
                const NodeId eventNode = events[i].node;
                uint32_t eventMask = 0;
                do {
                    eventMask |= events[i].mask;
                    ++i;
                } while (i < events.size() && events[i].node == eventNode);

                const uint32_t oldState = state;
                state |= eventMask;
                auto& bucket = buckets[static_cast<std::size_t>(eventNode)];
                if (!hasPreviousState) {
                    policy.applyInitial(bucket, state);
                    hasPreviousState = true;
                } else {
                    policy.applyTransition(bucket, oldState, state);
                }
            }
        }

        return buckets;
    }

    /**
     * @brief Computes local counters with a caller-owned bucket policy.
     *
     * @details
     * This is `computeDeltasWithPolicy(...)` followed by bottom-up subtree
     * aggregation with the policy's `merge` operation.
     *
     * Let `P` be the number of image pixels, `W` the window size, `N` the number
     * of internal node slots, and `A` the cost of owner/ancestor/LCA entry
     * evaluation. The structural cost is:
     *
     * `O(P * W * A + P * W log W + N + E)`,
     *
     * where `E` is the number of tree edges. Policy transition and merge costs
     * are multiplied into the corresponding event and aggregation terms.
     *
     * @tparam Policy Policy type defining `Bucket`, `applyInitial`,
     * `applyTransition`, and `merge`.
     * @param tree Tree whose dense internal node-id space indexes the returned
     * bucket vector.
     * @param window Finite local window. At most 32 samples are supported.
     * @param policy Stateless or externally owned policy object used to update
     * and merge buckets.
     * @return Dense final subtree bucket vector sized by
     * `tree.getNumInternalNodeSlots()`.
     * @throws std::invalid_argument If the tree/window shape is invalid.
     */
    template <class Policy>
    static std::vector<typename Policy::Bucket> computeWithPolicy(const MorphologicalTree& tree, const std::vector<WindowOffset>& window,
                                                                  const Policy& policy) {
        std::vector<typename Policy::Bucket> buckets = computeDeltasWithPolicy(tree, window, policy);
        aggregateSubtreeBuckets(tree, buckets, policy);
        return buckets;
    }

  private:
    /**
     * @brief Validates the tree/window shape required by the event engine.
     *
     * @param tree Tree topology used by the operation.
     * @param window Finite local sampling window.
     */
    static void requireValidInput(const MorphologicalTree& tree, const std::vector<WindowOffset>& window) {
        if (tree.getNumRowsOfGridDomain2D() <= 0 || tree.getNumColsOfGridDomain2D() <= 0) {
            throw std::invalid_argument("EventEngine requires a non-empty image domain.");
        }
        if (!tree.isAlive(tree.getRoot())) {
            throw std::invalid_argument("EventEngine requires a live tree root.");
        }
        if (window.size() > 32) {
            throw std::invalid_argument("EventEngine supports at most 32 window samples.");
        }
    }

    /**
     * @brief Orders events so descendants are processed before ancestors.
     *
     * @details
     * Valid entry events generated from one anchor lie on the anchor-owner to
     * root chain. The node-id fallback keeps the comparator strict and
     * deterministic if a caller ever passes a topology with stale or unusual
     * ancestry information.
     *
     * @param tree Tree topology used by the operation.
     * @param events Local events processed by the engine.
     */
    static void sortEventsBottomUp(const MorphologicalTree& tree, std::vector<EntryEvent>& events) {
        std::sort(events.begin(), events.end(), [&](const EntryEvent& lhs, const EntryEvent& rhs) {
            if (lhs.node == rhs.node) {
                return false;
            }
            if (tree.isStrictAncestor(lhs.node, rhs.node)) {
                return false;
            }
            if (tree.isStrictAncestor(rhs.node, lhs.node)) {
                return true;
            }
            return lhs.node < rhs.node;
        });
    }

    /**
     * @brief Aggregates local buckets bottom-up into subtree buckets.
     *
     * @details
     * Event processing records contributions at the node where a local state
     * starts or changes. The final per-node attribute value is a subtree
     * accumulation, so every child bucket is merged into its parent in
     * post-order using the policy-defined merge operation.
     *
     * @param tree Tree topology used by the operation.
     * @param values Values read or written by the operation.
     * @param policy Policy controlling the operation.
     */
    template <class Bucket, class Policy>
    static void aggregateSubtreeBuckets(const MorphologicalTree& tree, std::vector<Bucket>& values, const Policy& policy) {
        std::vector<std::pair<NodeId, bool>> stack;
        stack.emplace_back(tree.getRoot(), false);

        while (!stack.empty()) {
            const auto [node, expanded] = stack.back();
            stack.pop_back();
            if (!expanded) {
                stack.emplace_back(node, true);
                for (NodeId child : tree.getChildren(node)) {
                    stack.emplace_back(child, false);
                }
                continue;
            }

            for (NodeId child : tree.getChildren(node)) {
                policy.merge(values[static_cast<std::size_t>(node)], values[static_cast<std::size_t>(child)]);
            }
        }
    }
};

} // namespace mmcfilters::local_events
