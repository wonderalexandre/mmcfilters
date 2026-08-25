#pragma once

#include "BitquadAttributeData.hpp"
#include "../../../localAttributes/FiniteWindowLocalAttributeComputer.hpp"
#include "../../../trees/detail/TreeTraversalDetail.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::attributes::computers::detail {

namespace kernel {

/** @brief Pure local decision for signed bitquad-family contributions. */
struct BitquadFamilyLocalDecision {
    using Value = BitquadFamilyIncrement;

    std::size_t anchorPosition = 0; ///< Row-major position that owns the current framed cell.

    /** @brief Creates a rule for one row-major anchor position. @param position Position in `[0, 3]`. */
    explicit BitquadFamilyLocalDecision(std::size_t position) : anchorPosition(position) {}

    /** @brief Evaluates one owned canonical state. @param visibilityState Four-coordinate state. @return Unit family contribution or zero. */
    [[nodiscard]] Value evaluateLocalDecision(local_attributes::BinaryVisibilityState visibilityState) const {
        Value value;
        const BitquadCode code = static_cast<BitquadCode>(visibilityState.bits() & std::uint32_t{0b1111});
        if (!ownsState(code)) {
            return value;
        }
        const std::optional<BitquadFamily> family = bitquadFamily(code);
        if (!family.has_value()) {
            return value;
        }
        switch (*family) {
        case BitquadFamily::Q1:
            value.q1 = 1;
            break;
        case BitquadFamily::Q2:
            value.q2 = 1;
            break;
        case BitquadFamily::QD:
            value.qd = 1;
            break;
        case BitquadFamily::Q3:
            value.q3 = 1;
            break;
        case BitquadFamily::Q4:
            value.q4 = 1;
            break;
        }
        return value;
    }

  private:
    [[nodiscard]] bool ownsState(BitquadCode code) const {
        const std::uint32_t anchorMask = std::uint32_t{1} << anchorPosition;
        return (code & anchorMask) != 0 && (code & (anchorMask - std::uint32_t{1})) == 0;
    }
};

/** @brief Additive algebra for signed bitquad-family events. */
struct BitquadFamilyEventAlgebra {
    using Value = BitquadFamilyIncrement;

    [[nodiscard]] Value additiveIdentity() const { return {}; }

    void addAssign(Value& target, const Value& source) const {
        target.q1 += source.q1;
        target.q2 += source.q2;
        target.qd += source.qd;
        target.q3 += source.q3;
        target.q4 += source.q4;
    }

    /** @brief Subtracts one family increment. @param target Value to update. @param source Value to subtract. */
    void subtractAssign(Value& target, const Value& source) const {
        target.q1 -= source.q1;
        target.q2 -= source.q2;
        target.qd -= source.qd;
        target.q3 -= source.q3;
        target.q4 -= source.q4;
    }

};

/** @brief Pure local decision for signed nonempty-state histogram contributions. */
struct NonemptyBitquadStateHistogramLocalDecision {
    using Value = NonemptyBitquadStateHistogramIncrement;

    std::size_t anchorPosition = 0; ///< Row-major position that owns the current framed cell.

    /** @brief Creates a rule for one row-major anchor position. @param position Position in `[0, 3]`. */
    explicit NonemptyBitquadStateHistogramLocalDecision(std::size_t position) : anchorPosition(position) {}

    /** @brief Evaluates one owned canonical state. @param visibilityState Four-coordinate state. @return Unit state contribution or zero. */
    [[nodiscard]] Value evaluateLocalDecision(local_attributes::BinaryVisibilityState visibilityState) const {
        Value value;
        const BitquadCode code = static_cast<BitquadCode>(visibilityState.bits() & std::uint32_t{0b1111});
        if (ownsState(code)) {
            value.count(code) = 1;
        }
        return value;
    }

  private:
    [[nodiscard]] bool ownsState(BitquadCode code) const {
        const std::uint32_t anchorMask = std::uint32_t{1} << anchorPosition;
        return (code & anchorMask) != 0 && (code & (anchorMask - std::uint32_t{1})) == 0;
    }
};

/** @brief Additive algebra for signed nonempty-state histogram events. */
struct NonemptyBitquadStateHistogramEventAlgebra {
    using Value = NonemptyBitquadStateHistogramIncrement;

    [[nodiscard]] Value additiveIdentity() const { return {}; }

    void addAssign(Value& target, const Value& source) const {
        for (std::size_t index = 0; index < target.bins.size(); ++index) {
            target.bins[index] += source.bins[index];
        }
    }

    /** @brief Subtracts one histogram increment. @param target Value to update. @param source Value to subtract. */
    void subtractAssign(Value& target, const Value& source) const {
        for (std::size_t index = 0; index < target.bins.size(); ++index) {
            target.bins[index] -= source.bins[index];
        }
    }

};

/** @brief Row-major positions p0=TL, p1=TR, p2=BL, p3=BR. */
inline constexpr std::array<local_attributes::WindowOffset, 4> bitquadCellPositions{{{0, 0}, {0, 1}, {1, 0}, {1, 1}}};

/** @brief Immutable canonical windows indexed by row-major anchor position. */
inline const std::array<local_attributes::ObservationWindow, 4> bitquadObservationWindows{
    local_attributes::ObservationWindow{{{0, 0}, {0, 1}, {1, 0}, {1, 1}}},
    local_attributes::ObservationWindow{{{0, -1}, {0, 0}, {1, -1}, {1, 0}}},
    local_attributes::ObservationWindow{{{-1, 0}, {-1, 1}, {0, 0}, {0, 1}}},
    local_attributes::ObservationWindow{{{-1, -1}, {-1, 0}, {0, -1}, {0, 0}}},
};

/**
 * @brief Returns the immutable canonical observation window around one anchor position.
 * @param anchorPosition Row-major cell position in `[0, 3]`.
 * @return Offsets `p_l - p_anchorPosition` in canonical coordinate order.
 */
[[nodiscard]] inline const local_attributes::ObservationWindow& bitquadObservationWindow(std::size_t anchorPosition) {
    return bitquadObservationWindows[anchorPosition];
}

/**
 * @brief Computes the signed five-family increment at every dense node slot.
 * @param tree Established tree topology and grid domain.
 * @return Dense signed family increments.
 */
inline std::vector<BitquadFamilyIncrement> computeBitquadFamilyIncrements(const MorphologicalTree& tree) {
    std::vector<BitquadFamilyIncrement> increments(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (std::size_t anchorPosition = 0; anchorPosition < bitquadCellPositions.size(); ++anchorPosition) {
        local_attributes::detail::kernel::accumulateLocalAttributeIncrementValues(
            tree, bitquadObservationWindow(anchorPosition), BitquadFamilyLocalDecision{anchorPosition}, BitquadFamilyEventAlgebra{}, increments);
    }
    return increments;
}

/**
 * @brief Aggregates signed family increments into nonnegative family counts.
 * @param tree Established tree topology.
 * @param increments Dense signed family increments.
 * @return Dense aggregated family counts.
 */
inline std::vector<BitquadFamilyCounts> aggregateBitquadFamilyIncrements(const MorphologicalTree& tree, std::span<const BitquadFamilyIncrement> increments) {
    std::vector<BitquadFamilyCounts> counts;
    counts.reserve(increments.size());
    for (const BitquadFamilyIncrement& increment : increments) {
        counts.push_back({increment.q1, increment.q2, increment.qd, increment.q3, increment.q4});
    }
    ::mmcfilters::detail::kernel::traversePostOrder(
        tree, tree.root(), [](NodeId) {},
        [&](NodeId parent, NodeId child) {
            BitquadFamilyCounts& parentCounts = counts[static_cast<std::size_t>(parent)];
            const BitquadFamilyCounts& childCounts = counts[static_cast<std::size_t>(child)];
            parentCounts.q1 += childCounts.q1;
            parentCounts.q2 += childCounts.q2;
            parentCounts.qd += childCounts.qd;
            parentCounts.q3 += childCounts.q3;
            parentCounts.q4 += childCounts.q4;
        },
        [](NodeId) {});
    return counts;
}

/**
 * @brief Computes aggregated counts of the five nonempty bitquad families.
 * @param tree Established tree topology and grid domain.
 * @return Dense aggregated family counts.
 */
inline std::vector<BitquadFamilyCounts> computeBitquadFamilyCounts(const MorphologicalTree& tree) {
    return aggregateBitquadFamilyIncrements(tree, computeBitquadFamilyIncrements(tree));
}

/**
 * @brief Computes signed increments for canonical nonempty state codes 1 through 15.
 * @param tree Established tree topology and grid domain.
 * @return Dense signed 15-bin histogram increments.
 */
inline std::vector<NonemptyBitquadStateHistogramIncrement> computeNonemptyBitquadStateHistogramIncrements(const MorphologicalTree& tree) {
    std::vector<NonemptyBitquadStateHistogramIncrement> increments(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (std::size_t anchorPosition = 0; anchorPosition < bitquadCellPositions.size(); ++anchorPosition) {
        local_attributes::detail::kernel::accumulateLocalAttributeIncrementValues(
            tree, bitquadObservationWindow(anchorPosition), NonemptyBitquadStateHistogramLocalDecision{anchorPosition},
            NonemptyBitquadStateHistogramEventAlgebra{}, increments);
    }
    return increments;
}

/**
 * @brief Aggregates nonempty-state increments without inventing an empty-state increment.
 * @param tree Established tree topology.
 * @param increments Dense signed 15-bin histogram increments.
 * @return Dense aggregated 15-bin histograms.
 */
inline std::vector<NonemptyBitquadStateHistogram>
aggregateNonemptyBitquadStateHistogramIncrements(const MorphologicalTree& tree, std::span<const NonemptyBitquadStateHistogramIncrement> increments) {
    std::vector<NonemptyBitquadStateHistogram> histograms;
    histograms.reserve(increments.size());
    for (const NonemptyBitquadStateHistogramIncrement& increment : increments) {
        histograms.push_back({increment.bins});
    }
    ::mmcfilters::detail::kernel::traversePostOrder(
        tree, tree.root(), [](NodeId) {},
        [&](NodeId parent, NodeId child) {
            auto& parentBins = histograms[static_cast<std::size_t>(parent)].bins;
            const auto& childBins = histograms[static_cast<std::size_t>(child)].bins;
            for (std::size_t index = 0; index < parentBins.size(); ++index) {
                parentBins[index] += childBins[index];
            }
        },
        [](NodeId) {});
    return histograms;
}

/**
 * @brief Materializes full 16-bin histograms by deriving the empty count after aggregation.
 * @param numRows Number of image-domain rows.
 * @param numColumns Number of image-domain columns.
 * @param aliveNodes Dense identifiers of live nodes.
 * @param nonemptyHistograms Dense aggregated 15-bin histograms.
 * @return Dense full 16-bin histograms.
 */
inline std::vector<BitquadStateHistogram> materializeEmptyBitquadCount(int numRows, int numColumns, std::span<const NodeId> aliveNodes,
                                                                       std::span<const NonemptyBitquadStateHistogram> nonemptyHistograms) {
    std::vector<BitquadStateHistogram> histograms(nonemptyHistograms.size());
    const int totalCells = (numRows + 1) * (numColumns + 1);
    for (NodeId node : aliveNodes) {
        const NonemptyBitquadStateHistogram& nonempty = nonemptyHistograms[static_cast<std::size_t>(node)];
        BitquadStateHistogram& full = histograms[static_cast<std::size_t>(node)];
        int numNonemptyCells = 0;
        for (BitquadCode code = 1; code <= 15; ++code) {
            full.count(code) = nonempty.count(code);
            numNonemptyCells += nonempty.count(code);
        }
        full.count(0) = totalCells - numNonemptyCells;
    }
    return histograms;
}

} // namespace kernel

/**
 * @brief Canonical finite-window bitquad computation.
 *
 * @details
 * State coordinates use row-major order: bit 0 top-left, bit 1 top-right,
 * bit 2 bottom-left, and bit 3 bottom-right. The lowest-index visible bit owns
 * each nonempty framed cell. State and family increments are signed node
 * contributions; bottom-up aggregation produces counts. Code zero has no
 * family and its histogram count is materialized only after nonempty counts are
 * aggregated.
 */
class BitquadFiniteWindowComputation {
  public:
    using BitquadCode = ::mmcfilters::attributes::computers::detail::BitquadCode;                       ///< Canonical code type.
    using BitquadState = ::mmcfilters::attributes::computers::detail::BitquadState;                     ///< Spatial state type.
    using BitquadFamily = ::mmcfilters::attributes::computers::detail::BitquadFamily;                   ///< Five-family enum.
    using BitquadFamilyIncrement = ::mmcfilters::attributes::computers::detail::BitquadFamilyIncrement; ///< Signed family increment.
    using BitquadFamilyCounts = ::mmcfilters::attributes::computers::detail::BitquadFamilyCounts;       ///< Aggregated family counts.
    using NonemptyBitquadStateHistogramIncrement =
        ::mmcfilters::attributes::computers::detail::NonemptyBitquadStateHistogramIncrement;                          ///< Signed 15-bin increment.
    using NonemptyBitquadStateHistogram = ::mmcfilters::attributes::computers::detail::NonemptyBitquadStateHistogram; ///< Aggregated 15-bin histogram.
    using BitquadStateHistogram = ::mmcfilters::attributes::computers::detail::BitquadStateHistogram;                 ///< Full 16-bin histogram.
    using WindowOffset = local_attributes::WindowOffset;                                                              ///< Reused row-column offset type.

    /** @brief Computes and stores full state histograms and five-family counts. @param tree Tree topology and grid domain. */
    explicit BitquadFiniteWindowComputation(const MorphologicalTree& tree)
        : bitquadStateHistograms_(computeBitquadStateHistograms(tree)), bitquadFamilyCounts_(computeBitquadFamilyCounts(bitquadStateHistograms_)) {}

    /** @brief Returns cached full state histograms. @return Dense node-ordered histograms. */
    [[nodiscard]] const std::vector<BitquadStateHistogram>& bitquadStateHistograms() const noexcept { return bitquadStateHistograms_; }

    /** @brief Returns cached five-family counts. @return Dense node-ordered family counts. */
    [[nodiscard]] const std::vector<BitquadFamilyCounts>& bitquadFamilyCounts() const noexcept { return bitquadFamilyCounts_; }

    /** @brief Encodes a spatial state in canonical row-major order. @param state Spatial state. @return Canonical code. */
    [[nodiscard]] static constexpr BitquadCode bitquadCode(BitquadState state) noexcept {
        return ::mmcfilters::attributes::computers::detail::bitquadCode(state);
    }

    /** @brief Returns the family of a canonical code. @param code Code in `[0, 15]`. @return No value for zero; otherwise the family. */
    [[nodiscard]] static constexpr std::optional<BitquadFamily> bitquadFamily(BitquadCode code) {
        return ::mmcfilters::attributes::computers::detail::bitquadFamily(code);
    }

    /** @brief Returns the complete code-to-family table. @return Optional family for codes 0 through 15. */
    [[nodiscard]] static constexpr std::array<std::optional<BitquadFamily>, 16> bitquadFamilyTable() {
        return {{std::nullopt, BitquadFamily::Q1, BitquadFamily::Q1, BitquadFamily::Q2, BitquadFamily::Q1, BitquadFamily::Q2, BitquadFamily::QD,
                 BitquadFamily::Q3, BitquadFamily::Q1, BitquadFamily::QD, BitquadFamily::Q2, BitquadFamily::Q3, BitquadFamily::Q2, BitquadFamily::Q3,
                 BitquadFamily::Q3, BitquadFamily::Q4}};
    }

    /** @brief Finds the lowest visible bit index. @param code Canonical code. @return No value for zero; otherwise an index in `[0, 3]`. */
    [[nodiscard]] static std::optional<std::size_t> lowestVisibleBitIndex(BitquadCode code) noexcept {
        if (code == 0) {
            return std::nullopt;
        }
        return static_cast<std::size_t>(std::countr_zero(static_cast<unsigned int>(code)));
    }

    /**
     * @brief Returns the canonical observation window for one row-major anchor position.
     * @param anchorPosition Row-major anchor position in `[0, 3]`.
     * @return Canonically ordered four-coordinate observation window.
     */
    [[nodiscard]] static local_attributes::ObservationWindow bitquadObservationWindow(std::size_t anchorPosition) {
        if (anchorPosition >= kernel::bitquadCellPositions.size()) {
            throw std::out_of_range("Bitquad anchor position must be in [0, 3].");
        }
        return kernel::bitquadObservationWindow(anchorPosition);
    }

    /**
     * @brief Computes signed nonempty-state histogram increments for all dense node slots.
     * @param tree Tree topology and grid domain.
     * @return Dense signed 15-bin histogram increments.
     */
    [[nodiscard]] static std::vector<NonemptyBitquadStateHistogramIncrement> computeNonemptyBitquadStateHistogramIncrements(const MorphologicalTree& tree) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(local_attributes::detail::validateFiniteWindowLocalAttributeInput(tree));
        return kernel::computeNonemptyBitquadStateHistogramIncrements(tree);
    }

    /**
     * @brief Aggregates signed nonempty-state increments into 15-bin histograms.
     * @param tree Tree topology used for bottom-up aggregation.
     * @param increments Dense signed histogram increments.
     * @return Dense aggregated nonempty-state histograms.
     */
    [[nodiscard]] static std::vector<NonemptyBitquadStateHistogram>
    aggregateNonemptyBitquadStateHistogramIncrements(const MorphologicalTree& tree, std::span<const NonemptyBitquadStateHistogramIncrement> increments) {
        requireDenseShape(tree, increments.size(), "Nonempty bitquad state histogram increment aggregation");
        return kernel::aggregateNonemptyBitquadStateHistogramIncrements(tree, increments);
    }

    /**
     * @brief Derives bin zero by framed-grid complement after nonempty aggregation.
     * @param tree Tree topology and grid domain.
     * @param nonemptyHistograms Dense aggregated nonempty-state histograms.
     * @return Dense full 16-bin histograms.
     */
    [[nodiscard]] static std::vector<BitquadStateHistogram> materializeEmptyBitquadCount(const MorphologicalTree& tree,
                                                                                         std::span<const NonemptyBitquadStateHistogram> nonemptyHistograms) {
        requireDenseShape(tree, nonemptyHistograms.size(), "Empty bitquad count materialization");
        std::vector<NodeId> aliveNodes;
        aliveNodes.reserve(static_cast<std::size_t>(tree.numNodes()));
        for (NodeId node : tree.aliveNodeIds()) {
            aliveNodes.push_back(node);
        }
        return kernel::materializeEmptyBitquadCount(tree.numRows(), tree.numColumns(), aliveNodes, nonemptyHistograms);
    }

    /**
     * @brief Computes complete 16-bin canonical histograms for every dense node slot.
     * @param tree Tree topology and grid domain.
     * @return Dense full state histograms.
     */
    [[nodiscard]] static std::vector<BitquadStateHistogram> computeBitquadStateHistograms(const MorphologicalTree& tree) {
        const auto increments = computeNonemptyBitquadStateHistogramIncrements(tree);
        const auto nonemptyHistograms = aggregateNonemptyBitquadStateHistogramIncrements(tree, increments);
        return materializeEmptyBitquadCount(tree, nonemptyHistograms);
    }

    /**
     * @brief Computes signed five-family increments for every dense node slot.
     * @param tree Tree topology and grid domain.
     * @return Dense signed family increments.
     */
    [[nodiscard]] static std::vector<BitquadFamilyIncrement> computeBitquadFamilyIncrements(const MorphologicalTree& tree) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(local_attributes::detail::validateFiniteWindowLocalAttributeInput(tree));
        return kernel::computeBitquadFamilyIncrements(tree);
    }

    /**
     * @brief Aggregates signed five-family increments into family counts.
     * @param tree Tree topology used for bottom-up aggregation.
     * @param increments Dense signed family increments.
     * @return Dense aggregated family counts.
     */
    [[nodiscard]] static std::vector<BitquadFamilyCounts> aggregateBitquadFamilyIncrements(const MorphologicalTree& tree,
                                                                                           std::span<const BitquadFamilyIncrement> increments) {
        requireDenseShape(tree, increments.size(), "Bitquad family increment aggregation");
        return kernel::aggregateBitquadFamilyIncrements(tree, increments);
    }

    /**
     * @brief Computes aggregated counts of the five nonempty families.
     * @param tree Tree topology and grid domain.
     * @return Dense aggregated family counts.
     */
    [[nodiscard]] static std::vector<BitquadFamilyCounts> computeBitquadFamilyCounts(const MorphologicalTree& tree) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(local_attributes::detail::validateFiniteWindowLocalAttributeInput(tree));
        return kernel::computeBitquadFamilyCounts(tree);
    }

    /**
     * @brief Projects one full state histogram into five family counts, ignoring code zero.
     * @param histogram Full canonical state histogram.
     * @return Five nonempty-family counts.
     */
    [[nodiscard]] static BitquadFamilyCounts projectBitquadFamilyCounts(const BitquadStateHistogram& histogram) {
        BitquadFamilyCounts counts;
        for (BitquadCode code = 1; code <= 15; ++code) {
            const int count = histogram.count(code);
            switch (*bitquadFamily(code)) {
            case BitquadFamily::Q1:
                counts.q1 += count;
                break;
            case BitquadFamily::Q2:
                counts.q2 += count;
                break;
            case BitquadFamily::QD:
                counts.qd += count;
                break;
            case BitquadFamily::Q3:
                counts.q3 += count;
                break;
            case BitquadFamily::Q4:
                counts.q4 += count;
                break;
            }
        }
        return counts;
    }

    /**
     * @brief Projects full state histograms into five family counts.
     * @param histograms Full canonical state histograms.
     * @return Corresponding five-family counts.
     */
    [[nodiscard]] static std::vector<BitquadFamilyCounts> computeBitquadFamilyCounts(std::span<const BitquadStateHistogram> histograms) {
        std::vector<BitquadFamilyCounts> counts;
        counts.reserve(histograms.size());
        for (const BitquadStateHistogram& histogram : histograms) {
            counts.push_back(projectBitquadFamilyCounts(histogram));
        }
        return counts;
    }

    /**
     * @brief Projects one signed nonempty-state increment into a signed family increment.
     * @param histogramIncrement Signed 15-bin state increment.
     * @return Signed five-family increment.
     */
    [[nodiscard]] static BitquadFamilyIncrement projectBitquadFamilyIncrement(const NonemptyBitquadStateHistogramIncrement& histogramIncrement) {
        BitquadFamilyIncrement increment;
        for (BitquadCode code = 1; code <= 15; ++code) {
            const int value = histogramIncrement.count(code);
            switch (*bitquadFamily(code)) {
            case BitquadFamily::Q1:
                increment.q1 += value;
                break;
            case BitquadFamily::Q2:
                increment.q2 += value;
                break;
            case BitquadFamily::QD:
                increment.qd += value;
                break;
            case BitquadFamily::Q3:
                increment.q3 += value;
                break;
            case BitquadFamily::Q4:
                increment.q4 += value;
                break;
            }
        }
        return increment;
    }

    /**
     * @brief Projects per-node full histograms onto the smallest-node map.
     * @param tree Tree supplying the smallest-node map.
     * @param nodeHistograms Dense per-node full histograms.
     * @return One full histogram per proper part.
     */
    [[nodiscard]] static std::vector<BitquadStateHistogram> projectBitquadStateHistogramsToProperParts(const MorphologicalTree& tree,
                                                                                                       std::span<const BitquadStateHistogram> nodeHistograms) {
        requireDenseShape(tree, nodeHistograms.size(), "Bitquad state proper-part projection");
        return projectToSmallestNodeMap(tree, nodeHistograms);
    }

    /**
     * @brief Projects per-node family counts onto the smallest-node map.
     * @param tree Tree supplying the smallest-node map.
     * @param nodeCounts Dense per-node family counts.
     * @return One family-count value per proper part.
     */
    [[nodiscard]] static std::vector<BitquadFamilyCounts> projectBitquadFamilyCountsToProperParts(const MorphologicalTree& tree,
                                                                                                  std::span<const BitquadFamilyCounts> nodeCounts) {
        requireDenseShape(tree, nodeCounts.size(), "Bitquad family count proper-part projection");
        return projectToSmallestNodeMap(tree, nodeCounts);
    }

    /**
     * @brief Projects per-node family increments onto the smallest-node map.
     * @param tree Tree supplying the smallest-node map.
     * @param nodeIncrements Dense per-node family increments.
     * @return One family-increment value per proper part.
     */
    [[nodiscard]] static std::vector<BitquadFamilyIncrement>
    projectBitquadFamilyIncrementsToProperParts(const MorphologicalTree& tree, std::span<const BitquadFamilyIncrement> nodeIncrements) {
        requireDenseShape(tree, nodeIncrements.size(), "Bitquad family increment proper-part projection");
        return projectToSmallestNodeMap(tree, nodeIncrements);
    }

  private:
    /**
     * @brief Requires exact dense-node-slot coverage.
     * @param tree Tree defining the dense node-slot space.
     * @param size Number of supplied entries.
     * @param context Diagnostic operation name.
     */
    static void requireDenseShape(const MorphologicalTree& tree, std::size_t size, const char* context) {
        if (size != static_cast<std::size_t>(tree.numInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " requires exactly one entry per internal node slot.");
        }
    }

    /**
     * @brief Copies dense node values through the smallest-node map.
     * @tparam Value Node-value type.
     * @param tree Tree supplying the smallest-node map.
     * @param nodeValues Dense values indexed by node identifier.
     * @return Values indexed by proper-part identifier.
     */
    template <class Value> [[nodiscard]] static std::vector<Value> projectToSmallestNodeMap(const MorphologicalTree& tree, std::span<const Value> nodeValues) {
        std::vector<Value> projected(static_cast<std::size_t>(tree.numPixels()));
        for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
            const NodeId smallestNode = tree.smallestNode(pixel);
            if (smallestNode == InvalidNode || !tree.isAlive(smallestNode)) {
                throw std::runtime_error("Bitquad proper-part projection found a pixel without a live smallest node.");
            }
            projected[static_cast<std::size_t>(pixel)] = nodeValues[static_cast<std::size_t>(smallestNode)];
        }
        return projected;
    }

    std::vector<BitquadStateHistogram> bitquadStateHistograms_; ///< Cached full state histograms.
    std::vector<BitquadFamilyCounts> bitquadFamilyCounts_;      ///< Cached five-family counts.
};

} // namespace mmcfilters::attributes::computers::detail
