#include "support/TestSupport.hpp"

#include "mmcfilters/utils/RegularGridAdjacency2D.hpp"

#include <array>
#include <cmath>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

static_assert(std::ranges::forward_range<RegularGridAdjacency2D::AdjacentIndexRange>);
static_assert(std::ranges::forward_range<RegularGridAdjacency2D::NeighborIndexRange>);
static_assert(std::ranges::forward_range<RegularGridAdjacency2D::ForwardNeighborIndexRange>);
static_assert(std::is_same_v<decltype(std::declval<MorphologicalTree&>().getUniformGridAdjacency2D()), const RegularGridAdjacency2D*>);
static_assert(std::is_same_v<decltype(std::declval<MorphologicalTree&>().getUniformGridAdjacency2D()), const RegularGridAdjacency2D*>);
static_assert(std::is_same_v<decltype(std::declval<MorphologicalTree&>().getDirectionalGridAdjacency2D()), const DirectionalGridAdjacency2D*>);

namespace {

void testRadiusStencilOrderAndBounds() {
    const RegularGridAdjacency2D adjacency(3, 3, 1.5);

    requireVectorEqual(collectNodeIds(adjacency.getAdjacentIndices(1, 1)), std::vector<NodeId>{4, 3, 0, 1, 2, 5, 8, 7, 6},
                       "8-neighborhood order including origin");
    requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(4)), std::vector<NodeId>{3, 0, 1, 2, 5, 8, 7, 6}, "8-neighborhood order");
    requireVectorEqual(collectNodeIds(adjacency.getForwardNeighborIndices(4)), std::vector<NodeId>{5, 8, 7, 6}, "8-neighborhood forward order");
    requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(0)), std::vector<NodeId>{1, 4, 3}, "corner bounds filtering");

    requireThrows<std::out_of_range>([&] { (void)adjacency.getNeighborIndices(-1); }, "negative linear adjacency index");
    requireThrows<std::out_of_range>([&] { (void)adjacency.getAdjacentIndices(3, 0); }, "row outside adjacency domain");
}

void testNestedAndInterleavedTraversal() {
    const RegularGridAdjacency2D adjacency(3, 3, 1.5);
    const auto outer = adjacency.getNeighborIndices(4);
    auto outerIterator = outer.begin();
    requireEqual(*outerIterator, 3, "first outer neighbor");

    requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(0)), std::vector<NodeId>{1, 4, 3}, "nested traversal");
    ++outerIterator;
    requireEqual(*outerIterator, 0, "outer iterator survives nested traversal");

    const auto second = adjacency.getNeighborIndices(8);
    auto secondIterator = second.begin();
    requireEqual(*secondIterator, 7, "independent second traversal");
    ++outerIterator;
    requireEqual(*outerIterator, 1, "first traversal advances independently");
    ++secondIterator;
    requireEqual(*secondIterator, 4, "second traversal advances independently");
}

void testStructuringElementFactories() {
    const std::array<GridOffset2D, 5> cross{{
        {1, 0},
        {0, -1},
        {0, 0},
        {-1, 0},
        {0, 1},
    }};
    const auto custom = RegularGridAdjacency2D::fromStructuringElement(3, 3, std::span<const GridOffset2D>(cross));
    require(custom.getShape() == RegularGridAdjacencyShape::StructuringElement, "custom stencil shape");
    require(custom.is4connectivity(), "custom cross is 4-connectivity");
    require(!custom.is8connectivity(), "custom cross is not 8-connectivity");
    requireVectorEqual(collectNodeIds(custom.getNeighborIndices(4)), std::vector<NodeId>{3, 1, 5, 7}, "custom cross canonical order");
    require(custom.isAdjacent(4, 3), "custom cross contains horizontal offset");
    require(!custom.isAdjacent(4, 0), "custom cross excludes long diagonal offset");

    const auto rectangle = RegularGridAdjacency2D::rectangular(5, 5, 1, 2);
    requireEqual(rectangle.getSize(), 15, "rectangular stencil size");
    requireEqual(static_cast<int>(collectNodeIds(rectangle.getNeighborIndices(12)).size()), 14, "rectangular center neighbor count");
    requireEqual(static_cast<int>(collectNodeIds(rectangle.getForwardNeighborIndices(12)).size()), 7, "rectangular forward half");
    requireNear(rectangle.getRadius(), std::sqrt(5.0), 1e-12, "rectangular bounding radius");

    const auto horizontal = RegularGridAdjacency2D::horizontalLine(1, 5, 2);
    requireVectorEqual(collectNodeIds(horizontal.getNeighborIndices(2)), std::vector<NodeId>{1, 0, 3, 4}, "horizontal line");

    const auto vertical = RegularGridAdjacency2D::verticalLine(5, 1, 2);
    requireVectorEqual(collectNodeIds(vertical.getNeighborIndices(2)), std::vector<NodeId>{1, 0, 3, 4}, "vertical line");

    const auto diagonal = RegularGridAdjacency2D::line(7, 7, 2, 1);
    requireVectorEqual(collectNodeIds(diagonal.getNeighborIndices(24)), std::vector<NodeId>{16, 9, 32, 39}, "oriented digital line");
}

void testStructuringElementValidation() {
    const std::array<GridOffset2D, 2> asymmetric{{{0, 0}, {0, 1}}};
    requireThrowsContaining<std::invalid_argument>([&] { (void)RegularGridAdjacency2D::fromStructuringElement(3, 3, asymmetric); }, "centrally symmetric",
                                                   "asymmetric structuring element");

    const std::array<GridOffset2D, 2> missingOrigin{{{0, -1}, {0, 1}}};
    requireThrowsContaining<std::invalid_argument>([&] { (void)RegularGridAdjacency2D::fromStructuringElement(3, 3, missingOrigin); }, "origin",
                                                   "structuring element without origin");

    const std::array<GridOffset2D, 2> duplicateOrigin{{{0, 0}, {0, 0}}};
    requireThrowsContaining<std::invalid_argument>([&] { (void)RegularGridAdjacency2D::fromStructuringElement(3, 3, duplicateOrigin); }, "duplicate",
                                                   "duplicate structuring-element offset");
}

void testComponentTreeUsesExplicitAdjacency() {
    const auto image = makeImage(1, 5, {3, 1, 4, 1, 5});
    const auto horizontal = RegularGridAdjacency2D::horizontalLine(1, 5, 1);
    const auto tree = MorphologicalTreeFactory::createMaxTree(image, horizontal);
    const auto* stored = tree.topology().getUniformGridAdjacency2D();
    require(stored != nullptr, "component tree stores explicit adjacency");
    requireVectorEqual(collectNodeIds(stored->getNeighborIndices(2)), std::vector<NodeId>{1, 3}, "stored explicit adjacency");
    requireThrowsContaining<std::invalid_argument>([&] { (void)AttributeComputation::computeSingleTopologyAttribute(tree.topology(), BITQUADS_AREA); },
                                                   "canonical 4- or 8-connectivity", "BitQuads rejects unsupported custom adjacency");

    requireThrowsContaining<std::invalid_argument>(
        [&] { (void)MorphologicalTreeFactory::createMinTree(image, RegularGridAdjacency2D::rectangular(2, 5, 1, 1)); }, "must match",
        "component tree rejects adjacency from another domain");
}

void testCanonicalCustomStencilIsEquivalent() {
    const auto image = makeImage(3, 3,
                                 {
                                     3,
                                     2,
                                     2,
                                     3,
                                     4,
                                     1,
                                     0,
                                     4,
                                     1,
                                 });
    const std::array<GridOffset2D, 5> cross{{
        {0, 0},
        {-1, 0},
        {0, -1},
        {0, 1},
        {1, 0},
    }};

    const auto radiusTree = MorphologicalTreeFactory::createMaxTree(image, 1.0);
    const auto customTree = MorphologicalTreeFactory::createMaxTree(image, RegularGridAdjacency2D::fromStructuringElement(3, 3, cross));
    const auto radiusExport = exportHigraHierarchy(radiusTree);
    const auto customExport = exportHigraHierarchy(customTree);
    requireVectorEqual(customExport.first, radiusExport.first, "custom 4-connectivity parent equivalence");
    requireVectorEqual(customExport.second, radiusExport.second, "custom 4-connectivity altitude equivalence");
    for (NodeId properPart = 0; properPart < radiusTree.topology().getNumTotalProperParts(); ++properPart) {
        requireEqual(customTree.topology().getProperPartOwner(properPart), radiusTree.topology().getProperPartOwner(properPart),
                     "custom 4-connectivity owner equivalence");
    }
}

} // namespace

int main() {
    testRadiusStencilOrderAndBounds();
    testNestedAndInterleavedTraversal();
    testStructuringElementFactories();
    testStructuringElementValidation();
    testComponentTreeUsesExplicitAdjacency();
    testCanonicalCustomStencilIsEquivalent();
    return 0;
}
