#pragma once

/**
 * @file SaturatedDynamicLca.hpp
 * @brief Dynamic LCA strategies used exclusively by saturated certification.
 */

#include "../ResidualTreePolicies.hpp"
#include "../../ValuedMorphologicalTree.hpp"
#include "../../../utils/GenerationStampSet.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Dynamic lowest-common-ancestor strategies for one saturated residual altitude type. */
template <AltitudeValue T> struct SaturatedLcaTypes {
    /** @brief Valued component-tree type indexed by the LCA strategies. */
    using Tree = ValuedMorphologicalTree<T>;


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
            /** @brief Dense node identifier of the node identifier. */
            NodeId nodeId = InvalidNode;
            /** @brief Dense node identifier of the next child. */
            NodeId nextChild = InvalidNode;
            /** @brief Depth. */
            int depth = 0;
        };

        /** @brief References the tree used by the component. */
        const Tree* tree_ = nullptr;
        /** @brief Indicates whether the blocked snapshot index is enabled. */
        bool enabled_ = false;
        /** @brief Indicates whether the blocked snapshot requires a rebuild. */
        bool rebuildPending_ = true;
        /** @brief Mutations since rebuild. */
        std::size_t mutationsSinceRebuild_ = 0;
        /** @brief Block size. */
        std::size_t blockSize_ = 1;
        /** @brief Snapshot alive nodes. */
        std::size_t snapshotAliveNodes_ = 0;
        /** @brief Dense node identifier of the euler. */
        std::vector<NodeId> euler_;
        /** @brief Euler depth buffer. */
        std::vector<int> eulerDepth_;
        /** @brief First occurrence buffer. */
        std::vector<int> firstOccurrence_;
        /** @brief Preorder entry buffer. */
        std::vector<int> preorderEntry_;
        /** @brief Preorder exit buffer. */
        std::vector<int> preorderExit_;
        /** @brief Rmq log2 buffer. */
        std::vector<int> rmqLog2_;
        /** @brief Rmq sparse table buffer. */
        std::vector<int> rmqSparseTable_;
        /** @brief Rmq stride. */
        std::size_t rmqStride_ = 0;
        /** @brief Rmq levels. */
        int rmqLevels_ = 0;
        /** @brief Dirty preorder buffer. */
        std::vector<std::uint8_t> dirtyPreorder_;
        /** @brief Next unmarked buffer. */
        std::vector<int> nextUnmarked_;
        /** @brief Dirty root marks. */
        GenerationStampSet dirtyRootMarks_;
        /** @brief Dense node identifier of the mirrored parent. */
        std::vector<NodeId> mirroredParent_;
        /** @brief Mirrored alive buffer. */
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
            const std::size_t numSlots = static_cast<std::size_t>(topology.numInternalNodeSlots());
            firstOccurrence_.assign(numSlots, -1);
            preorderEntry_.assign(numSlots, -1);
            preorderExit_.assign(numSlots, -1);
            euler_.clear();
            eulerDepth_.clear();
            euler_.reserve(static_cast<std::size_t>(std::max(1, 2 * topology.numNodes() - 1)));
            eulerDepth_.reserve(euler_.capacity());

            std::vector<TraversalFrame> stack;
            stack.reserve(static_cast<std::size_t>(std::max(1, topology.numNodes())));
            int preorder = 0;
            const auto pushNode = [&](NodeId nodeId, int depth) {
                const std::size_t index = static_cast<std::size_t>(nodeId);
                preorderEntry_[index] = preorder++;
                firstOccurrence_[index] = static_cast<int>(euler_.size());
                appendEuler(nodeId, depth);
                stack.push_back(TraversalFrame{nodeId, topology.getFirstChild(nodeId), depth});
            };
            pushNode(topology.root(), 0);
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

            snapshotAliveNodes_ = static_cast<std::size_t>(topology.numNodes());
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
        void bind(const Tree& tree, SaturatedMinMaxLcaPolicy policy) {
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
                const std::size_t numSlots = static_cast<std::size_t>(topology.numInternalNodeSlots());
                mirroredParent_.assign(numSlots, InvalidNode);
                mirroredAlive_.assign(numSlots, std::uint8_t{0});
                for (NodeId nodeId : topology.aliveNodeIds()) {
                    const std::size_t index = static_cast<std::size_t>(nodeId);
                    mirroredAlive_[index] = 1;
                    const NodeId parent = topology.parent(nodeId);
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
                    const NodeId parent = topology.parent(nodeId);
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
            /** @brief Dense node identifier of the left. */
            NodeId left = InvalidNode;
            /** @brief Dense node identifier of the right. */
            NodeId right = InvalidNode;
            /** @brief Dense node identifier of the auxiliary parent. */
            NodeId auxiliaryParent = InvalidNode;
        };

        /** @brief References the tree used by the component. */
        const Tree* tree_ = nullptr;
        /** @brief Indicates whether the link-cut index is enabled. */
        bool enabled_ = false;
        /** @brief Nodes buffer. */
        std::vector<LinkCutNode> nodes_;
        /** @brief Dense node identifier of the represented parent. */
        std::vector<NodeId> representedParent_;
        /** @brief Represented alive buffer. */
        std::vector<std::uint8_t> representedAlive_;
        /** @brief Dense node identifier of the changed scratch. */
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
        void bind(const Tree& tree, SaturatedMinMaxLcaPolicy policy) {
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
            const std::size_t numSlots = static_cast<std::size_t>(topology.numInternalNodeSlots());
            nodes_.assign(numSlots, LinkCutNode{});
            representedParent_.assign(numSlots, InvalidNode);
            representedAlive_.assign(numSlots, std::uint8_t{0});
            for (NodeId nodeId : topology.aliveNodeIds()) {
                const std::size_t index = static_cast<std::size_t>(nodeId);
                representedAlive_[index] = 1;
                const NodeId parent = topology.parent(nodeId);
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
                    const NodeId parent = topology.parent(nodeId);
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
                    const NodeId parent = topology.parent(nodeId);
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
        /** @brief Policy. */
        SaturatedMinMaxLcaPolicy policy_ = SaturatedMinMaxLcaPolicy::ParentClimb;
        /** @brief Blocked. */
        BlockedDynamicLcaIndex blocked_;
        /** @brief Link cut. */
        LinkCutDynamicLcaIndex linkCut_;

      public:
        /**
         * @brief Binds the dynamic index to the current tree topology.
         *
         * @param tree Tree processed by the operation.
         * @param policy Policy controlling the operation.
         */
        void bind(const Tree& tree, SaturatedMinMaxLcaPolicy policy) {
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

};

} // namespace mmcfilters::sdrt::detail
