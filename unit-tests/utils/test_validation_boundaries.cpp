#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/utils/Contract.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

ImageUInt8Ptr makeImage() {
    auto image = ImageUInt8::create(5, 5);
    const std::uint8_t values[] = {1, 1, 1, 1, 1, 1, 3, 3, 3, 1, 1, 3, 7, 3, 1, 1, 3, 3, 3, 1, 1, 1, 1, 1, 1};
    for (int pixel = 0; pixel < image->getSize(); ++pixel) {
        (*image)[pixel] = values[pixel];
    }
    return image;
}

long double requireScientificVerticalSlice() {
    const auto image = makeImage();
    const auto weighted = MorphologicalTreeFactory::createMaxTree(image, 1.0);
    const auto maxDist = AttributeComputation::computeSingleAttribute<double>(weighted, MAX_DIST);
    long double checksum = 0;

    requireEqual(maxDist.first.NUM_ATTRIBUTES, 1, "MAX_DIST layout width");
    requireEqual(maxDist.second.size(), static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()), "MAX_DIST dense node domain");
    for (NodeId node : weighted.topology().getAliveNodeIds()) {
        const double value = maxDist.second[static_cast<std::size_t>(node)];
        require(std::isfinite(value), "MAX_DIST values must be finite");
        checksum += static_cast<long double>(value) * static_cast<long double>(node + 1);
    }

    const std::vector<AttributeOrGroup> aggregateRequest{AREA,
                                                         VOLUME,
                                                         RELATIVE_VOLUME,
                                                         LEVEL,
                                                         MEAN_LEVEL,
                                                         VARIANCE_LEVEL,
                                                         GRAY_HEIGHT,
                                                         HEIGHT_NODE,
                                                         NUM_DESCENDANTS_NODE,
                                                         BOX_WIDTH,
                                                         RECTANGULARITY,
                                                         CENTRAL_MOMENT_20,
                                                         HU_MOMENT_1,
                                                         INERTIA,
                                                         CONTOUR_PERIMETER,
                                                         BITQUADS_AREA};
    const auto aggregates = AttributeComputation::computeAttributes<double>(weighted, aggregateRequest);
    requireEqual(aggregates.first.NUM_ATTRIBUTES, 16, "aggregate kernel-chain layout width");
    requireEqual(aggregates.second.size(),
                 static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots() * aggregates.first.NUM_ATTRIBUTES),
                 "aggregate kernel-chain dense node domain");
    constexpr std::array<Attribute, 16> aggregateAttributes{AREA,
                                                            VOLUME,
                                                            RELATIVE_VOLUME,
                                                            LEVEL,
                                                            MEAN_LEVEL,
                                                            VARIANCE_LEVEL,
                                                            GRAY_HEIGHT,
                                                            HEIGHT_NODE,
                                                            NUM_DESCENDANTS_NODE,
                                                            BOX_WIDTH,
                                                            RECTANGULARITY,
                                                            CENTRAL_MOMENT_20,
                                                            HU_MOMENT_1,
                                                            INERTIA,
                                                            CONTOUR_PERIMETER,
                                                            BITQUADS_AREA};
    for (NodeId node : weighted.topology().getAliveNodeIds()) {
        for (std::size_t attributeIndex = 0; attributeIndex < aggregateAttributes.size(); ++attributeIndex) {
            const Attribute attribute = aggregateAttributes[attributeIndex];
            const double value = aggregates.second[static_cast<std::size_t>(aggregates.first.linearIndex(node, attribute))];
            require(std::isfinite(value), "aggregate kernel-chain values must be finite");
            checksum += static_cast<long double>(value) * static_cast<long double>(node + 1) * static_cast<long double>(attributeIndex + 2);
        }
    }

    std::vector<bool> keep(static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()), true);
    AttributeFilters<std::uint8_t> filters(weighted);
    const auto filtered = filters.filteringByDirectRule(keep);
    require(filtered->isEqual(image), "direct filtering with an all-true criterion must reconstruct the source image");
    for (int pixel = 0; pixel < filtered->getSize(); ++pixel) {
        checksum += static_cast<long double>((*filtered)[pixel]) * static_cast<long double>(pixel + 1);
    }
    return checksum;
}

void requireCheckedBoundaryDiagnostics() {
    if constexpr (contract::validationsEnabled) {
        const auto weighted = MorphologicalTreeFactory::createMaxTree(makeImage(), 1.0);
        AttributeFilters<std::uint8_t> filters(weighted);
        std::vector<bool> shortCriterion(1, true);
        requireThrows<std::invalid_argument>([&] { static_cast<void>(filters.filteringByDirectRule(shortCriterion)); },
                                             "checked direct-filter boundary must reject a short criterion");
        requireThrows<std::invalid_argument>([&] { static_cast<void>(weighted.topology().getNodeParent(InvalidNode)); },
                                             "checked primitive tree access must reject an invalid node");
        requireThrows<std::invalid_argument>([&] { static_cast<void>(MorphologicalTreeFactory::createMaxTree(makeImage(), 0.0)); },
                                             "checked construction boundary must reject a disconnected adjacency");
    }
}

} // namespace

int main() {
#if MMCFILTERS_CONTRACT_MODE == MMCFILTERS_CONTRACT_CHECKED
    static_assert(contract::validationsEnabled);
#else
    static_assert(!contract::validationsEnabled);
#endif
    const long double checksum = requireScientificVerticalSlice();
    requireCheckedBoundaryDiagnostics();
    std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << "scientific-slice-checksum=" << checksum << '\n';
    return 0;
}
