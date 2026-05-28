#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/Attributes.hpp"

#include <algorithm>
#include <cctype>
#include <concepts>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

#ifndef MMCFILTERS_TEST_DATA_DIR
#define MMCFILTERS_TEST_DATA_DIR "dat"
#endif

ImageUInt8Ptr makeConstantImage(int rows, int cols, std::uint8_t value) {
    return ImageUInt8::create(rows, cols, value);
}

ImageUInt8Ptr makeConcentricImage(bool brightCenter) {
    constexpr int rows = 8;
    constexpr int cols = 8;
    auto image = ImageUInt8::create(rows, cols);

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int ring = std::min({row, col, rows - 1 - row, cols - 1 - col});
            const int value = brightCenter ? ring * 32 : (3 - ring) * 32;
            (*image)[ImageUtils::to1D(row, col, cols)] = static_cast<std::uint8_t>(value);
        }
    }

    return image;
}

ImageUInt8Ptr makeLineImage(bool horizontal, bool brightLine) {
    constexpr int rows = 9;
    constexpr int cols = 9;
    const std::uint8_t foreground = brightLine ? 200 : 0;
    const std::uint8_t background = brightLine ? 0 : 200;
    auto image = ImageUInt8::create(rows, cols, background);

    if (horizontal) {
        const int row = rows / 2;
        for (int col = 1; col + 1 < cols; ++col) {
            (*image)[ImageUtils::to1D(row, col, cols)] = foreground;
        }
    } else {
        const int col = cols / 2;
        for (int row = 1; row + 1 < rows; ++row) {
            (*image)[ImageUtils::to1D(row, col, cols)] = foreground;
        }
    }

    return image;
}

ImageUInt8Ptr makeDiagonalImage() {
    constexpr int rows = 9;
    constexpr int cols = 9;
    auto image = ImageUInt8::create(rows, cols, static_cast<std::uint8_t>(0));
    for (int i = 1; i + 1 < rows; ++i) {
        (*image)[ImageUtils::to1D(i, i, cols)] = 200;
    }
    return image;
}

ImageUInt8Ptr makeTwoPlateauImage(bool brightPeaks) {
    constexpr int rows = 11;
    constexpr int cols = 11;
    auto image = ImageUInt8::create(rows, cols, static_cast<std::uint8_t>(brightPeaks ? 0 : 96));

    auto fillBox = [&](int rowBegin, int rowEnd, int colBegin, int colEnd, std::uint8_t value) {
        for (int row = rowBegin; row <= rowEnd; ++row) {
            for (int col = colBegin; col <= colEnd; ++col) {
                (*image)[ImageUtils::to1D(row, col, cols)] = value;
            }
        }
    };

    constexpr std::uint8_t mid = 48;
    const std::uint8_t extreme = brightPeaks ? 96 : 0;

    fillBox(1, 4, 1, 4, mid);
    fillBox(6, 9, 6, 9, mid);
    fillBox(2, 3, 2, 3, extreme);
    fillBox(7, 8, 7, 8, extreme);

    return image;
}

ImageUInt8Ptr makeRingImage(bool brightRing) {
    constexpr int rows = 9;
    constexpr int cols = 9;
    auto image = ImageUInt8::create(rows, cols, static_cast<std::uint8_t>(brightRing ? 0 : 160));
    const std::uint8_t ringValue = brightRing ? 160 : 0;

    for (int i = 1; i + 1 < rows; ++i) {
        (*image)[ImageUtils::to1D(1, i, cols)] = ringValue;
        (*image)[ImageUtils::to1D(rows - 2, i, cols)] = ringValue;
        (*image)[ImageUtils::to1D(i, 1, cols)] = ringValue;
        (*image)[ImageUtils::to1D(i, cols - 2, cols)] = ringValue;
    }

    return image;
}

ImageUInt8Ptr makeCheckerboardImage() {
    constexpr int rows = 8;
    constexpr int cols = 8;
    auto image = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            (*image)[ImageUtils::to1D(row, col, cols)] =
                ((row + col) % 2 == 0) ? static_cast<std::uint8_t>(220) : static_cast<std::uint8_t>(20);
        }
    }
    return image;
}

std::string readPgmToken(std::istream& input) {
    std::string token;
    char ch = '\0';
    while (input.get(ch)) {
        if (std::isspace(static_cast<unsigned char>(ch))) {
            continue;
        }
        if (ch == '#') {
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        token.push_back(ch);
        break;
    }

    while (input.get(ch)) {
        if (std::isspace(static_cast<unsigned char>(ch))) {
            break;
        }
        if (ch == '#') {
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            break;
        }
        token.push_back(ch);
    }

    if (token.empty()) {
        throw std::runtime_error("Unexpected end of PGM file.");
    }
    return token;
}

ImageUInt8Ptr loadPgmAsSampledUInt8(const std::filesystem::path& path, int targetRows, int targetCols) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Failed to open PGM file: " + path.string());
    }

    const std::string magic = readPgmToken(input);
    const int cols = std::stoi(readPgmToken(input));
    const int rows = std::stoi(readPgmToken(input));
    const int maxValue = std::stoi(readPgmToken(input));
    require(cols > 0 && rows > 0 && maxValue > 0, "PGM dimensions and max value must be positive");

    std::vector<int> source(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols), 0);
    if (magic == "P5") {
        for (int index = 0; index < rows * cols; ++index) {
            if (maxValue < 256) {
                char value = '\0';
                input.get(value);
                source[static_cast<std::size_t>(index)] =
                    static_cast<unsigned char>(value);
            } else {
                char high = '\0';
                char low = '\0';
                input.get(high);
                input.get(low);
                source[static_cast<std::size_t>(index)] =
                    (static_cast<unsigned char>(high) << 8) | static_cast<unsigned char>(low);
            }
            if (!input) {
                throw std::runtime_error("Truncated PGM pixel data: " + path.string());
            }
        }
    } else if (magic == "P2") {
        for (int index = 0; index < rows * cols; ++index) {
            source[static_cast<std::size_t>(index)] = std::stoi(readPgmToken(input));
        }
    } else {
        throw std::runtime_error("Unsupported PGM magic in " + path.string() + ": " + magic);
    }

    const int outRows = std::min(rows, targetRows);
    const int outCols = std::min(cols, targetCols);
    auto output = ImageUInt8::create(outRows, outCols);
    for (int row = 0; row < outRows; ++row) {
        const int sourceRow = row * rows / outRows;
        for (int col = 0; col < outCols; ++col) {
            const int sourceCol = col * cols / outCols;
            const int raw = source[static_cast<std::size_t>(ImageUtils::to1D(sourceRow, sourceCol, cols))];
            const int scaled = (raw * 255 + (maxValue / 2)) / maxValue;
            (*output)[ImageUtils::to1D(row, col, outCols)] =
                static_cast<std::uint8_t>(std::clamp(scaled, 0, 255));
        }
    }
    return output;
}

struct NumericFixture {
    std::string label;
    ImageUInt8Ptr image;
    bool isMaxTree;
    double adjacencyRadius;
};

std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> makeWeightedTree(const NumericFixture& fixture) {
    return makeWeightedComponentTree(fixture.image, fixture.isMaxTree, fixture.adjacencyRadius);
}

std::string treeKindLabel(bool isMaxTree) {
    return isMaxTree ? "max-tree" : "min-tree";
}

template <std::floating_point Real>
void requireAllAttributesAreFinite(
    const NumericFixture& fixture,
    const AttributeNames& names,
    std::span<const Real> values,
    int numRows,
    const std::string& outputSpaceLabel,
    std::unordered_map<Attribute, bool>& observedFiniteValues) {
    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::ALL);

    requireEqual(
        names.NUM_ATTRIBUTES,
        static_cast<int>(allAttributes.size()),
        fixture.label + " " + outputSpaceLabel + " ALL group stride");
    requireEqual(
        values.size(),
        static_cast<std::size_t>(numRows) * allAttributes.size(),
        fixture.label + " " + outputSpaceLabel + " ALL group buffer size");

    for (Attribute attribute : allAttributes) {
        require(
            names.contains(attribute),
            fixture.label + " " + outputSpaceLabel + " missing attribute " + AttributeNames::toString(attribute));
    }

    for (int row = 0; row < numRows; ++row) {
        for (Attribute attribute : allAttributes) {
            const Real value = values[static_cast<std::size_t>(names.linearIndex(row, attribute))];
            if (std::isfinite(value)) {
                observedFiniteValues[attribute] = true;
                continue;
            }

            std::ostringstream message;
            message << fixture.label << " " << outputSpaceLabel
                    << " " << treeKindLabel(fixture.isMaxTree)
                    << " radius=" << fixture.adjacencyRadius
                    << " row=" << row
                    << " attribute=" << AttributeNames::toString(attribute)
                    << " produced a non-finite value " << value;
            throw std::runtime_error(message.str());
        }
    }
}

template <std::floating_point Real>
void requireAllAttributesAreFiniteOnInternalNodes(
    const NumericFixture& fixture,
    std::unordered_map<Attribute, bool>& observedFiniteValues) {
    const auto weighted = makeWeightedTree(fixture);
    const MorphologicalTree& tree = weighted->topology();
    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::ALL);

    const auto computed = AttributeComputation::computeAttributes<Real>(*weighted, {AttributeGroup::ALL});
    const AttributeNames& names = computed.attributeNames();
    const std::vector<Real>& values = computed.values();

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : allAttributes) {
            const Real value = values[names.linearIndex(nodeId, attribute)];
            if (std::isfinite(value)) {
                observedFiniteValues[attribute] = true;
                continue;
            }

            std::ostringstream message;
            message << fixture.label << " internal "
                    << treeKindLabel(fixture.isMaxTree)
                    << " radius=" << fixture.adjacencyRadius
                    << " node=" << nodeId
                    << " attribute=" << AttributeNames::toString(attribute)
                    << " produced a non-finite value " << value;
            throw std::runtime_error(message.str());
        }
    }
}

template <std::floating_point Real>
void requireAllAttributesAreFiniteInExportedHigraSpace(
    const NumericFixture& fixture,
    std::unordered_map<Attribute, bool>& observedFiniteValues) {
    const auto weighted = makeWeightedTree(fixture);
    const auto [parent, altitude] = weighted->exportHigraHierarchy();
    auto imported = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(parent),
        std::span<const std::uint8_t>(altitude),
        weighted->topology().getNumRowsOfImage(),
        weighted->topology().getNumColsOfImage(),
        fixture.isMaxTree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE,
        AdjacencyRelation(
            weighted->topology().getNumRowsOfImage(),
            weighted->topology().getNumColsOfImage(),
            fixture.adjacencyRadius));

    const auto computed = AttributeComputation::computeAttributes<Real>(
        imported,
        {AttributeGroup::ALL},
        NodeIdSpace::HIGRA);
    requireAllAttributesAreFinite<Real>(
        fixture,
        computed.attributeNames(),
        computed.values(),
        imported.topology().getNodeIdSpaceSize(NodeIdSpace::HIGRA),
        "exported-higra",
        observedFiniteValues);
}

void addTreeCases(
    std::vector<NumericFixture>& fixtures,
    const std::string& label,
    ImageUInt8Ptr image,
    std::span<const double> radii) {
    for (double radius : radii) {
        fixtures.push_back({label + "-max-r" + std::to_string(radius), image, true, radius});
        fixtures.push_back({label + "-min-r" + std::to_string(radius), image, false, radius});
    }
}

} // namespace

int main() {
    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::ALL);
    std::unordered_map<Attribute, bool> observedFiniteValues;
    for (Attribute attribute : allAttributes) {
        observedFiniteValues[attribute] = false;
    }

    std::vector<NumericFixture> fixtures;
    const std::vector<double> componentTreeRadii{1.0, 1.5};
    addTreeCases(fixtures, "unit", makeConstantImage(1, 1, 7), componentTreeRadii);
    addTreeCases(fixtures, "constant", makeConstantImage(5, 5, 42), componentTreeRadii);
    addTreeCases(fixtures, "concentric-bright", makeConcentricImage(true), componentTreeRadii);
    addTreeCases(fixtures, "concentric-dark", makeConcentricImage(false), componentTreeRadii);
    addTreeCases(fixtures, "horizontal-line-bright", makeLineImage(true, true), componentTreeRadii);
    addTreeCases(fixtures, "vertical-line-dark", makeLineImage(false, false), componentTreeRadii);
    addTreeCases(fixtures, "diagonal-bright", makeDiagonalImage(), componentTreeRadii);
    addTreeCases(fixtures, "two-plateau-bright", makeTwoPlateauImage(true), componentTreeRadii);
    addTreeCases(fixtures, "two-plateau-dark", makeTwoPlateauImage(false), componentTreeRadii);
    addTreeCases(fixtures, "ring-bright", makeRingImage(true), componentTreeRadii);
    addTreeCases(fixtures, "ring-dark", makeRingImage(false), componentTreeRadii);
    addTreeCases(fixtures, "checkerboard", makeCheckerboardImage(), componentTreeRadii);

    const auto dataDir = std::filesystem::path(MMCFILTERS_TEST_DATA_DIR);
    addTreeCases(fixtures, "real-lena", loadPgmAsSampledUInt8(dataDir / "lena.pgm", 48, 48), componentTreeRadii);
    addTreeCases(fixtures, "real-brain", loadPgmAsSampledUInt8(dataDir / "brain2.pgm", 48, 48), componentTreeRadii);
    addTreeCases(fixtures, "real-wrist", loadPgmAsSampledUInt8(dataDir / "wrist.pgm", 48, 48), componentTreeRadii);

    for (const NumericFixture& fixture : fixtures) {
        requireAllAttributesAreFiniteOnInternalNodes<float>(fixture, observedFiniteValues);
        requireAllAttributesAreFiniteInExportedHigraSpace<float>(fixture, observedFiniteValues);
        requireAllAttributesAreFiniteOnInternalNodes<double>(fixture, observedFiniteValues);
        requireAllAttributesAreFiniteInExportedHigraSpace<double>(fixture, observedFiniteValues);
    }

    for (Attribute attribute : allAttributes) {
        require(
            observedFiniteValues.at(attribute),
            "attribute " + AttributeNames::toString(attribute) + " must have at least one finite value across numeric fixtures");
    }

    return 0;
}
