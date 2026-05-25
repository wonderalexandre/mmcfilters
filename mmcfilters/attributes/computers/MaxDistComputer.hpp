#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "../../utils/Common.hpp"
#include "../../contours/ContoursComputedIncrementally.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"

#include "detail/maxdist/EdtDIFT.hpp"

namespace mmcfilters::attributes::computers {
    /**
     * @brief Computes the MAX_DIST attribute for every node of a morphological tree.
     *
     * MAX_DIST is evaluated by sweeping tree nodes in typed altitude order and
     * maintaining an incremental squared Euclidean distance transform over the
     * support accumulated at the current level. Local contour additions and
     * removals are shared with `ContoursComputedIncrementally`, avoiding a
     * second dense level-image contour pass. The resulting attribute value is
     * the largest squared distance reached from the node contour.
     *
     * The output vector is indexed by internal node id and sized with
     * `tree_.getNumInternalNodeSlots()`, matching the attribute-computation
     * conventions used elsewhere in mmcfilters.
     */
    class MaxDistComputer
    {
    public:
      /**
       * @brief Creates a MAX_DIST computer bound to one tree.
       *
       * `altitude` is the optional uint8 altitude view used by getAttributes().
       * Type-generic callers should prefer the span overload.
       *
       * @param tree Topology whose node supports define the MAX_DIST domains.
       * @param altitude Optional canonical uint8 altitude view indexed by dense
       * internal node id.
       */
      MaxDistComputer(const MorphologicalTree& tree, OptionalAltitudeSpan<std::uint8_t> altitude = {})
        :tree_{tree}, altitude_{altitude}
      {}

      /**
       * @brief Rejects tree kinds for which the altitude sweep is not defined.
       *
       * @throws std::invalid_argument If the tree kind is not `MAX_TREE` or
       * `MIN_TREE`.
       */
      static void requireSupportedTreeKind(const MorphologicalTree& tree)
      {
        switch (tree.getTreeType()) {
          case MorphologicalTreeKind::MAX_TREE:
          case MorphologicalTreeKind::MIN_TREE:
            return;
          case MorphologicalTreeKind::TREE_OF_SHAPES:
            throw std::invalid_argument(
              "MAX_DIST is currently defined only for MAX_TREE and MIN_TREE; TREE_OF_SHAPES is not supported.");
          case MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE:
            throw std::invalid_argument(
              "MAX_DIST is currently defined only for MAX_TREE and MIN_TREE; SELF_DUAL_RESIDUAL_TREE is not supported.");
        }

        throw std::invalid_argument("MAX_DIST received an unsupported morphological tree kind.");
      }

      /**
       * @brief Computes attributes using the constructor-provided uint8 altitude buffer.
       *
       * @throws std::logic_error when no altitude view was provided to the constructor.
       * @throws std::runtime_error when the stored altitude view has an invalid shape.
       */
      [[nodiscard]] std::vector<float> getAttributes() const
      {
        return getAttributes(TreeAltitudeAlgorithms::requireAltitudeSpan(altitude_));
      }

      /**
       * @brief Computes attributes using an explicit typed altitude span.
       *
       * The altitude span must cover the tree's internal node-id space. Values
       * may be integer or floating point, but floating-point altitudes must be
       * finite because the level ordering and contour tests require a total
       * finite order.
       *
       * @param altitude Dense altitude span indexed by internal node id.
       * @return Dense `MAX_DIST` vector indexed by internal node id.
       * @throws std::invalid_argument If the tree kind, altitude shape, or
       * floating-point altitude values are invalid.
       */
      template<AltitudeValue T>
      [[nodiscard]] std::vector<float> getAttributes(std::span<const T> altitude) const
      {
        requireSupportedTreeKind(tree_);
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree_, altitude);
        return getAttributesBySortedLevels(altitude);
      }

    private:
      using ContourDeltaStore = ::mmcfilters::ContoursComputedIncrementally::LocalContourDeltas;

      /**
       * @brief True when an altitude value can participate in level ordering.
       */
      template<AltitudeValue T>
      static bool isFiniteAltitude(T value) noexcept
      {
        if constexpr (std::is_floating_point_v<T>) {
          return std::isfinite(value);
        }
        return true;
      }

      /**
       * @brief Rejects NaN or infinite floating-point altitude values.
       */
      template<AltitudeValue T>
      static void validateFiniteAltitude(std::span<const T> altitude)
      {
        for (std::size_t index = 0; index < altitude.size(); ++index) {
          if (!isFiniteAltitude(altitude[index])) {
            throw std::invalid_argument(
              "MAX_DIST requires finite altitude values; node " +
              std::to_string(index) + " has a non-finite altitude.");
          }
        }
      }

      static void markPixels(std::span<const int> pixels, std::vector<uint8_t>& marks)
      {
        for (int pixelId : pixels) {
          marks[static_cast<std::size_t>(pixelId)] = 1;
        }
      }

      static void clearPixelMarks(std::span<const int> pixels, std::vector<uint8_t>& marks)
      {
        for (int pixelId : pixels) {
          marks[static_cast<std::size_t>(pixelId)] = 0;
        }
      }

      /**
       * @brief Groups nodes by sorted altitude, independently of the altitude
       * value domain.
       *
       * The vector is ordered by the level sweep required by MAX_DIST: descending
       * for max-trees and ascending for min-trees. Stable sorting preserves
       * post-order inside equal-altitude groups, so children are processed before
       * parents when a tree contains flat parent-child levels.
       */
      template<AltitudeValue T>
      std::vector<NodeId> sortedNodesByAltitude(std::span<const T> altitude) const
      {
        std::vector<NodeId> nodes;
        nodes.reserve(static_cast<std::size_t>(tree_.getNumNodes()));
        for (NodeId nodeId : tree_.getPostOrderNodes()) {
          nodes.push_back(nodeId);
        }

        std::stable_sort(nodes.begin(), nodes.end(), [&](NodeId lhs, NodeId rhs) {
          const T lhsAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, lhs);
          const T rhsAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, rhs);
          return tree_.getTreeType() == MorphologicalTreeKind::MAX_TREE
            ? lhsAltitude > rhsAltitude
            : lhsAltitude < rhsAltitude;
        });

        return nodes;
      }

      /**
       * @brief Tests whether two nodes belong to the same altitude group.
       */
      template<AltitudeValue T>
      bool sameAltitude(std::span<const T> altitude, NodeId lhs, NodeId rhs) const
      {
        return TreeAltitudeAlgorithms::getAltitude(altitude, lhs) ==
               TreeAltitudeAlgorithms::getAltitude(altitude, rhs);
      }

      /**
       * @brief Processes one group of nodes that share the same altitude.
       *
       * Children contours are inherited into the parent contour unless the
       * shared contour-delta store says the pixel must be removed at this node.
       * Proper parts are then added as contour seeds when they are local
       * additions, or opened as interior pixels otherwise. Only after all nodes
       * in the level group have been materialized does EdtDIFT propagate labels;
       * this preserves the simultaneous per-level sweep.
       */
      void processLevel(
        std::span<const NodeId> nodes,
        const ContourDeltaStore& contourDeltas,
        detail::maxdist::EdtDIFT& edtDIFT,
        std::vector<std::vector<int>>& contours,
        std::vector<uint8_t>& removalMark,
        std::vector<uint8_t>& contourAdditionMark,
        std::vector<float>& maxDist) const
      {
        if (nodes.empty()) {
          return;
        }

        std::vector<int> toRemove;
        toRemove.reserve(64);
        for (NodeId nodeId : nodes) {
          std::vector<int>& nodeContour = contours[static_cast<std::size_t>(nodeId)];
          nodeContour.clear();

          const auto removals = contourDeltas.removals(nodeId);
          toRemove.clear();
          toRemove.insert(toRemove.end(), removals.begin(), removals.end());
          markPixels(removals, removalMark);

          const auto additions = contourDeltas.additions(nodeId);
          std::size_t reserveSize = additions.size();
          for (NodeId childNodeId : tree_.getChildren(nodeId)) {
            reserveSize += contours[static_cast<std::size_t>(childNodeId)].size();
          }
          nodeContour.reserve(reserveSize);

          for (NodeId childNodeId : tree_.getChildren(nodeId)) {
            std::vector<int>& childContour = contours[static_cast<std::size_t>(childNodeId)];
            for (int pixelId : childContour) {
              if (!removalMark[static_cast<std::size_t>(pixelId)]) {
                nodeContour.push_back(pixelId);
              }
            }
            std::vector<int>().swap(childContour);
          }
          clearPixelMarks(removals, removalMark);

          if (!toRemove.empty()) {
            edtDIFT.treeRemoval(toRemove);
          }

          markPixels(additions, contourAdditionMark);
          for (int pixelId : tree_.getProperParts(nodeId)) {
            edtDIFT.addPixelToBinaryImage(pixelId);

            if (contourAdditionMark[static_cast<std::size_t>(pixelId)]) {
              nodeContour.push_back(pixelId);
              edtDIFT.seed(pixelId);
            }
            else {
              edtDIFT.open(pixelId);
              edtDIFT.insertNeighborsPQueue(pixelId);
            }
          }
          clearPixelMarks(additions, contourAdditionMark);
        }

        // All nodes at the same altitude must enter the binary image before the
        // distance transform propagates. This preserves the mathematical
        // per-level schedule without requiring fixed 0..255 buckets.
        edtDIFT.run();

        for (NodeId nodeId : nodes) {
          maxDist[nodeId] = edtDIFT.maxBedt(contours[static_cast<std::size_t>(nodeId)]);
        }
      }

      /**
       * @brief Runs the complete altitude sweep and returns the attribute vector.
       */
      template<AltitudeValue T>
      std::vector<float> getAttributesBySortedLevels(std::span<const T> altitude) const
      {
        validateFiniteAltitude(altitude);

        std::vector<float> maxDist(tree_.getNumInternalNodeSlots(), 0.0f);
        detail::maxdist::EdtDIFT edtDIFT(tree_.getNumRowsOfImage(), tree_.getNumColsOfImage());

        const ContourDeltaStore contourDeltas =
          ::mmcfilters::ContoursComputedIncrementally::extractContourDeltas(tree_);
        std::vector<std::vector<int>> contours(static_cast<std::size_t>(tree_.getNumInternalNodeSlots()));
        const std::size_t totalPixels = static_cast<std::size_t>(
          tree_.getNumRowsOfImage() * tree_.getNumColsOfImage());
        std::vector<uint8_t> removalMark(totalPixels, 0);
        std::vector<uint8_t> contourAdditionMark(totalPixels, 0);

        std::vector<NodeId> sortedNodes = sortedNodesByAltitude(altitude);
        std::size_t groupBegin = 0;
        while (groupBegin < sortedNodes.size()) {
          std::size_t groupEnd = groupBegin + 1;
          while (groupEnd < sortedNodes.size() &&
                 sameAltitude(altitude, sortedNodes[groupBegin], sortedNodes[groupEnd])) {
            ++groupEnd;
          }

          processLevel(
            std::span<const NodeId>(
              sortedNodes.data() + static_cast<std::ptrdiff_t>(groupBegin),
              groupEnd - groupBegin),
            contourDeltas,
            edtDIFT,
            contours,
            removalMark,
            contourAdditionMark,
            maxDist);

          groupBegin = groupEnd;
        }

        return maxDist;
      }


    private:
      /**
       * @brief Topology whose node supports define the MAX_DIST domains.
       */
      const MorphologicalTree& tree_;

      /**
       * @brief Optional uint8 altitude view for getAttributes().
       */
      OptionalAltitudeSpan<std::uint8_t> altitude_;
    };
} // namespace mmcfilters::attributes::computers
