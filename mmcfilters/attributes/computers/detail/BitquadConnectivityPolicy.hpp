#pragma once

#include "../../../trees/MorphologicalTree.hpp"
#include "../../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../../trees/detail/MorphologicalTreeConstructionContextQueries.hpp"

#include <optional>
#include <span>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail {

/** @brief Polarity of one non-root shape relative to its parent shape. */
enum class ShapePolarity {
    Lower, ///< The node altitude is strictly lower than its parent altitude.
    Upper  ///< The node altitude is strictly higher than its parent altitude.
};

/** @brief Foreground connectivity selected for bitquad scalar projection. */
enum class BitquadConnectivity { Four, Eight };

/**
 * @brief Explicit connectivity choices used after bitquad-family counting.
 *
 * @details
 * The root has no shape polarity and therefore has its own connectivity entry.
 * Lower and upper entries apply only to non-root shapes. This policy belongs
 * to scalar materialization; it is independent of finite-window computation and
 * bitquad-family aggregation.
 */
class BitquadConnectivityPolicy {
  private:
    BitquadConnectivity rootConnectivity_ = BitquadConnectivity::Eight;       ///< Connectivity for the unpolarized root.
    BitquadConnectivity lowerShapeConnectivity_ = BitquadConnectivity::Eight; ///< Connectivity for lower shapes.
    BitquadConnectivity upperShapeConnectivity_ = BitquadConnectivity::Eight; ///< Connectivity for upper shapes.
    bool shapePolarityRequired_ = false;                                      ///< Whether lower/upper choices differ.

  public:
    /** @brief Constructs explicit root, lower-shape, and upper-shape choices. @param rootConnectivity Root choice. @param lowerShapeConnectivity Lower-shape choice. @param upperShapeConnectivity Upper-shape choice. @param shapePolarityRequired Whether non-root polarity is required. */
    BitquadConnectivityPolicy(BitquadConnectivity rootConnectivity, BitquadConnectivity lowerShapeConnectivity,
                              BitquadConnectivity upperShapeConnectivity, bool shapePolarityRequired) noexcept
        : rootConnectivity_(rootConnectivity), lowerShapeConnectivity_(lowerShapeConnectivity),
          upperShapeConnectivity_(upperShapeConnectivity), shapePolarityRequired_(shapePolarityRequired) {}

    /** @brief Returns whether non-root shape polarity is needed to select connectivity. @return True when lower/upper choices differ. */
    [[nodiscard]] bool requiresShapePolarity() const noexcept { return shapePolarityRequired_; }

    /** @brief Returns the connectivity used by the unpolarized root. @return Root connectivity. */
    [[nodiscard]] BitquadConnectivity rootConnectivity() const noexcept { return rootConnectivity_; }

    /** @brief Returns the connectivity used by lower shapes. @return Lower-shape connectivity. */
    [[nodiscard]] BitquadConnectivity lowerShapeConnectivity() const noexcept { return lowerShapeConnectivity_; }

    /** @brief Returns the connectivity used by upper shapes. @return Upper-shape connectivity. */
    [[nodiscard]] BitquadConnectivity upperShapeConnectivity() const noexcept { return upperShapeConnectivity_; }

    /**
     * @brief Selects connectivity from an optional shape polarity.
     * @param shapePolarity `std::nullopt` for the unpolarized root; otherwise the non-root shape polarity.
     * @return Connectivity used by the scalar bitquad formulas.
     */
    [[nodiscard]] BitquadConnectivity connectivity(std::optional<ShapePolarity> shapePolarity) const noexcept {
        if (!shapePolarity.has_value()) {
            return rootConnectivity_;
        }
        return *shapePolarity == ShapePolarity::Upper ? upperShapeConnectivity_ : lowerShapeConnectivity_;
    }

    /** @brief Tests whether the selected scalar projection uses four-connectivity. @param shapePolarity Optional non-root polarity. @return True for four-connectivity. */
    [[nodiscard]] bool uses4Connectivity(std::optional<ShapePolarity> shapePolarity) const noexcept {
        return connectivity(shapePolarity) == BitquadConnectivity::Four;
    }
};

/** @brief Converts a validated canonical grid adjacency into a bitquad connectivity. */
inline BitquadConnectivity bitquadConnectivity(const RegularGridAdjacency2D& adjacency) {
    if (!adjacency.isCanonical4Or8Connectivity()) {
        throw std::invalid_argument("Bitquad scalar projection requires canonical 4- or 8-connectivity.");
    }
    return adjacency.is4connectivity() ? BitquadConnectivity::Four : BitquadConnectivity::Eight;
}

/**
 * @brief Builds the scalar-projection policy from typed construction context.
 * @param tree Tree whose retained construction context supplies the adjacency choices.
 * @param nodeAltitudesAvailable Whether a node-altitude buffer is available for shape-polarity selection.
 * @return Explicit root/lower/upper bitquad connectivity policy.
 */
inline BitquadConnectivityPolicy makeBitquadConnectivityPolicy(const MorphologicalTree& tree, bool nodeAltitudesAvailable) {
    if (const RegularGridAdjacency2D* adjacency = ::mmcfilters::detail::constructionAdjacency(tree)) {
        const BitquadConnectivity connectivity = bitquadConnectivity(*adjacency);
        return {connectivity, connectivity, connectivity, false};
    }

    const auto projectionAdjacencies = ::mmcfilters::detail::currentBitquadProjectionAdjacencies(tree);
    if (!projectionAdjacencies.has_value()) {
        throw std::invalid_argument("Bitquad scalar attributes require a construction adjacency or topographic convention.");
    }

    const BitquadConnectivity lowerShapeConnectivity = bitquadConnectivity(projectionAdjacencies->minAdjacency);
    const BitquadConnectivity upperShapeConnectivity = bitquadConnectivity(projectionAdjacencies->maxAdjacency);
    const bool shapePolarityRequired = lowerShapeConnectivity != upperShapeConnectivity;
    if (shapePolarityRequired && !nodeAltitudesAvailable) {
        throw std::invalid_argument("Bitquad scalar projection requires node altitudes when lower- and upper-shape connectivity differ.");
    }

    // The root is not an upper or lower shape. Preserve the established
    // eight-connected root projection when the complementary choices differ;
    // equal choices naturally remain uniform at the root.
    const BitquadConnectivity rootConnectivity =
        shapePolarityRequired ? BitquadConnectivity::Eight : lowerShapeConnectivity;
    return {rootConnectivity, lowerShapeConnectivity, upperShapeConnectivity, shapePolarityRequired};
}

/**
 * @brief Derives one non-root shape's polarity from exact node altitudes.
 * @param tree Tree containing the node and its parent relation.
 * @param nodeAltitudes Dense node-altitude buffer.
 * @param node Node whose polarity is requested.
 * @return No value for the root; otherwise `Lower` or `Upper`.
 * @throws std::runtime_error If a non-root node has the same altitude as its parent.
 */
template <AltitudeValue T>
[[nodiscard]] std::optional<ShapePolarity> shapePolarity(const MorphologicalTree& tree, std::span<const T> nodeAltitudes, NodeId node) {
    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, nodeAltitudes);
    if (tree.isRoot(node)) {
        return std::nullopt;
    }

    const NodeId parent = tree.parent(node);
    const T nodeAltitude = nodeAltitudes[static_cast<std::size_t>(node)];
    const T parentAltitude = nodeAltitudes[static_cast<std::size_t>(parent)];
    if (nodeAltitude > parentAltitude) {
        return ShapePolarity::Upper;
    }
    if (nodeAltitude < parentAltitude) {
        return ShapePolarity::Lower;
    }
    throw std::runtime_error("Bitquad scalar projection cannot derive shape polarity from equal node and parent altitudes.");
}

} // namespace mmcfilters::attributes::computers::detail
