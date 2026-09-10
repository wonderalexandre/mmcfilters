#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/computers/FilledShapeAttributeComputer.hpp"

#include <array>
#include <numbers>
#include <random>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;
using attributes::computers::FilledShapeAttributeComputer;

namespace {

// Raster oracle: flood the complementary background from a padded image frame.
std::vector<PixelId> filledPixels(const MorphologicalTree& tree, NodeId node, bool foregroundFour) {
    const int rows = tree.numRows() + 2;
    const int columns = tree.numColumns() + 2;
    std::vector<unsigned char> occupied(static_cast<std::size_t>(rows * columns), 0);
    for (NodeId descendant : tree.subtreeNodes(node)) {
        for (PixelId pixel : tree.properPart(descendant)) {
            occupied[static_cast<std::size_t>((pixel / tree.numColumns() + 1) * columns + pixel % tree.numColumns() + 1)] = 1;
        }
    }
    std::vector<int> queue{0};
    occupied[0] = 2;
    for (std::size_t i = 0; i < queue.size(); ++i) {
        const int row = queue[i] / columns;
        const int column = queue[i] % columns;
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if ((dx == 0 && dy == 0) || (!foregroundFour && dx != 0 && dy != 0)) continue;
                const int y = row + dy, x = column + dx;
                if (y < 0 || x < 0 || y >= rows || x >= columns) continue;
                const int next = y * columns + x;
                if (occupied[static_cast<std::size_t>(next)] == 0) {
                    occupied[static_cast<std::size_t>(next)] = 2;
                    queue.push_back(next);
                }
            }
        }
    }
    std::vector<PixelId> result;
    for (int y = 1; y + 1 < rows; ++y) {
        for (int x = 1; x + 1 < columns; ++x) {
            if (occupied[static_cast<std::size_t>(y * columns + x)] != 2) {
                result.push_back((y - 1) * tree.numColumns() + x - 1);
            }
        }
    }
    return result;
}

template <AltitudeValue T>
void verifyRasterOracle(const ValuedMorphologicalTree<T>& valued, bool foregroundFour, bool directional = false) {
    const MorphologicalTree& tree = valued.topology();
    const auto result = AttributeComputation::computeAttributes<double>(valued, {AttributeGroup::FilledShape});
    requireEqual(result.first.NUM_ATTRIBUTES, 12, "filled family column count");
    const auto& names = result.first;
    for (NodeId node : tree.aliveNodeIds()) {
        auto value = [&](Attribute attribute) { return result.second[names.linearIndex(node, attribute)]; };
        bool nodeFour = foregroundFour;
        if (directional && !tree.isRoot(node)) {
            const bool lower = valued.nodeAltitude(node) < valued.nodeAltitude(tree.parent(node));
            nodeFour = lower ? foregroundFour : !foregroundFour;
        }
        const auto pixels = filledPixels(tree, node, nodeFour);
        const double area = static_cast<double>(pixels.size());
        double column = 0, row = 0;
        for (PixelId pixel : pixels) {
            column += pixel % tree.numColumns();
            row += pixel / tree.numColumns();
        }
        column /= area;
        row /= area;
        double mu20 = 0, mu02 = 0, mu11 = 0;
        for (PixelId pixel : pixels) {
            const double x = pixel % tree.numColumns() - column;
            const double y = pixel / tree.numColumns() - row;
            mu20 += x * x;
            mu02 += y * y;
            mu11 += x * y;
        }
        double supportArea = 0, supportColumn = 0, supportRow = 0;
        for (NodeId descendant : tree.subtreeNodes(node)) {
            for (PixelId pixel : tree.properPart(descendant)) {
                ++supportArea;
                supportColumn += pixel % tree.numColumns();
                supportRow += pixel / tree.numColumns();
            }
        }
        const double discriminant = std::sqrt((mu20 - mu02) * (mu20 - mu02) + 4 * mu11 * mu11);
        const double lambda1 = mu20 + mu02 + discriminant;
        const double lambda2 = std::max(0., mu20 + mu02 - discriminant);
        const double eccentricity = lambda1 <= 1e-12 ? 1 : (lambda2 <= 1e-12 ? 1e6 : std::min(lambda1 / lambda2, 1e6));
        const double angle = (std::abs(mu20 - mu02) + std::abs(mu11) <= 1e-12) ? 0 :
            std::abs(0.5 * std::atan2(2 * mu11, mu20 - mu02)) * 180 / std::numbers::pi;
        requireNear(value(FilledArea), area, 1e-10, "raster filled area");
        requireNear(value(FilledCentroidColumn), column, 1e-10, "raster filled column centroid");
        requireNear(value(FilledCentroidRow), row, 1e-10, "raster filled row centroid");
        requireNear(value(FilledLengthMajorAxis), std::sqrt(2 * lambda1 / area), 1e-9, "raster major axis");
        requireNear(value(FilledLengthMinorAxis), std::sqrt(2 * lambda2 / area), 1e-7, "raster minor axis");
        requireNear(value(FilledAxisOrientation), angle, 1e-9, "raster orientation");
        requireNear(value(FilledEccentricity), eccentricity, 1e-7, "raster eccentricity");
        requireNear(value(FilledInertia), (mu20 + mu02) / (area * area), 1e-10, "raster inertia");
        const double compactness = mu20 + mu02 > 1e-12 ? area / (2 * std::numbers::pi * (mu20 + mu02)) : 0;
        const double circularity = lambda1 <= 1e-12 ? 1 : (lambda2 <= 1e-12 ? 0 : lambda2 / lambda1);
        requireNear(value(FilledCompactness), compactness, 1e-10, "raster compactness");
        requireNear(value(FilledCircularity), circularity, 1e-10, "raster circularity");
        requireNear(value(HoleAreaFraction), 1 - supportArea / area, 1e-10, "raster hole fraction");
        requireNear(value(FilledCentroidDisplacementNormalized),
                    std::hypot(column - supportColumn / supportArea, row - supportRow / supportArea) / std::sqrt(area),
                    1e-10, "raster normalized centroid displacement");
    }
}

void verifyGeometryAndRouting() {
    auto ring = ImageUInt8::create(8, 10, 0);
    for (int row = 1; row <= 5; ++row) {
        for (int column = 2; column <= 8; ++column) {
            (*ring)[row * 10 + column] = 2;
        }
    }
    (*ring)[2 * 10 + 3] = 0;
    (*ring)[3 * 10 + 3] = 0;
    (*ring)[4 * 10 + 7] = 0;
    for (double radius : {1., 1.5}) {
        auto valued = makeValuedComponentTree(ring, true, radius);
        verifyRasterOracle(*valued, radius == 1.);
        const auto topology = AttributeComputation::computeTopologyAttributes<double>(valued->topology(), {AttributeGroup::FilledShape});
        const auto typed = AttributeComputation::computeAttributes<double>(*valued, {AttributeGroup::FilledShape});
        require(topology.second == typed.second, "uniform topology and valued routing agree");
        const auto scalar = AttributeComputation::computeSingleAttribute<float>(*valued, FilledArea);
        requireEqual(scalar.first.NUM_ATTRIBUTES, 1, "individual requests expose no extra columns");
        for (NodeId node : valued->topology().aliveNodeIds()) {
            requireNear(static_cast<double>(scalar.second[scalar.first.linearIndex(node, FilledArea)]),
                        typed.second[typed.first.linearIndex(node, FilledArea)], 1e-6, "float area scalar");
        }
    }

    // Filled moment conventions agree with the existing point-mass attributes.
    auto rectangle = makeValuedComponentTree(ImageUInt8::create(3, 5, 1), true);
    const auto filled = AttributeComputation::computeAttributes<double>(*rectangle, {AttributeGroup::FilledShape});
    const std::array<std::pair<Attribute, Attribute>, 8> pairs{{
        {FilledArea, Area}, {FilledLengthMajorAxis, LengthMajorAxis}, {FilledLengthMinorAxis, LengthMinorAxis},
        {FilledAxisOrientation, AxisOrientation}, {FilledEccentricity, Eccentricity}, {FilledInertia, Inertia},
        {FilledCompactness, Compactness}, {FilledCircularity, Circularity}}};
    for (auto [filledAttribute, ordinaryAttribute] : pairs) {
        const auto ordinary = AttributeComputation::computeSingleAttribute<double>(*rectangle, ordinaryAttribute);
        const NodeId root = rectangle->topology().root();
        requireNear(filled.second[filled.first.linearIndex(root, filledAttribute)],
                    ordinary.second[ordinary.first.linearIndex(root, ordinaryAttribute)], 1e-10, "ordinary moment convention");
        const auto scalar = AttributeComputation::computeSingleAttribute<double>(*rectangle, filledAttribute);
        requireNear(scalar.second[scalar.first.linearIndex(root, filledAttribute)],
                    filled.second[filled.first.linearIndex(root, filledAttribute)], 1e-10, "scalar moment selection");
    }
}

void verifyRandomAndDegenerateSupports() {
    for (auto image : {makeImage(1, 1, {4}), makeImage(1, 5, {0, 2, 2, 2, 0}),
                       makeImage(5, 1, {0, 2, 2, 2, 0}), makeImage(3, 3, {1, 0, 2, 0, 2, 0, 0, 0, 0}),
                       makeImage(2, 2, {2, 0, 0, 2})}) {
        for (double radius : {1., 1.5}) {
            for (bool maximum : {false, true}) {
                verifyRasterOracle(*makeValuedComponentTree(image, maximum, radius), radius == 1.);
            }
        }
    }
    std::mt19937 random(711);
    for (int iteration = 0; iteration < 40; ++iteration) {
        auto image = ImageUInt8::create(6, 7);
        for (PixelId pixel = 0; pixel < 42; ++pixel) (*image)[pixel] = static_cast<uint8_t>(random() % 5);
        for (double radius : {1., 1.5}) {
            for (bool maximum : {false, true}) {
                verifyRasterOracle(*makeValuedComponentTree(image, maximum, radius), radius == 1.);
            }
        }
    }
}

void verifyGroupsAndUnitRows() {
    const auto& group = ATTRIBUTE_GROUPS.at(AttributeGroup::FilledShape);
    require(std::equal(group.begin(), group.end(), FilledShapeAttributeComputer::producedAttributes.begin(),
                       FilledShapeAttributeComputer::producedAttributes.end()), "complete filled group");
    for (AttributeGroup other : {AttributeGroup::Shape, AttributeGroup::Moments, AttributeGroup::Boundary}) {
        for (Attribute attribute : group) {
            const auto& members = ATTRIBUTE_GROUPS.at(other);
            require(std::find(members.begin(), members.end(), attribute) == members.end(), "exclusive thematic group");
        }
    }
    auto tree = makeValuedComponentTree(makeComponentTreeFixture(), true);
    const auto names = AttributeNames::fromList(group);
    const std::vector<PixelId> pixels{10, 0, 15};
    std::vector<double> buffer(pixels.size() * group.size(), -1);
    FilledShapeAttributeComputer::computeUnitRows(UnitAttributeComputeContext<double>{tree->topology(), pixels, buffer, names, group});
    for (NodeId row = 0; row < 3; ++row) {
        for (Attribute attribute : group) {
            double expected = 0;
            if (attribute == FilledArea || attribute == FilledEccentricity || attribute == FilledCircularity) expected = 1;
            if (attribute == FilledCentroidRow) expected = pixels[static_cast<std::size_t>(row)] / 4;
            if (attribute == FilledCentroidColumn) expected = pixels[static_cast<std::size_t>(row)] % 4;
            requireEqual(buffer[names.linearIndex(row, attribute)], expected, "unit row in supplied pixel order");
        }
    }
}

void verifyShapeConnectivityAndUnsupportedSupports() {
    for (auto immersion : {TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
        for (bool upper : {false, true}) {
            auto image = ImageUInt8::create(5, 5, uint8_t{1});
            (*image)[6] = (*image)[12] = upper ? 2 : 0;
            auto valued = MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(image, makeTopographicConvention(image, immersion));
            verifyRasterOracle(valued, immersion == TestTopographicImmersion::Min4Max8, true);
            bool rejected = false;
            try {
                static_cast<void>(AttributeComputation::computeSingleTopologyAttribute(valued.topology(), FilledArea));
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            require(rejected, "unequal directional adjacency requires altitudes");
        }
    }
    auto disconnected = makeTreeFromHigraParent({3, 4, 3, 4, 4}, 1, 3, true, 1.);
    requireThrowsContaining<std::logic_error>(
        [&] { static_cast<void>(AttributeComputation::computeSingleTopologyAttribute(*disconnected, FilledArea)); },
        "requires exactly one external boundary", "multiple external boundaries cannot define one filled region");
}

} // namespace

int main() {
    verifyGeometryAndRouting();
    verifyRandomAndDegenerateSupports();
    verifyGroupsAndUnitRows();
    verifyShapeConnectivityAndUnsupportedSupports();
    return 0;
}
