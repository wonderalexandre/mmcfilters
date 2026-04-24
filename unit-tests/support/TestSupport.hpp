#pragma once

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "mmcfilters/utils/Common.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"

namespace mmcfilters::unit_tests {

inline void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class T, class U>
inline void requireEqual(const T &actual, const U &expected, const std::string &label) {
    if (!(actual == expected)) {
        std::ostringstream oss;
        oss << label << ": expected `" << expected << "` but got `" << actual << "`";
        throw std::runtime_error(oss.str());
    }
}

template <class T>
inline void requireVectorEqual(const std::vector<T> &actual, const std::vector<T> &expected, const std::string &label) {
    if (actual == expected) {
        return;
    }

    std::ostringstream oss;
    oss << label << ": expected [";
    for (size_t i = 0; i < expected.size(); ++i) {
        if (i != 0) {
            oss << ", ";
        }
        oss << +expected[i];
    }
    oss << "] but got [";
    for (size_t i = 0; i < actual.size(); ++i) {
        if (i != 0) {
            oss << ", ";
        }
        oss << +actual[i];
    }
    oss << "]";
    throw std::runtime_error(oss.str());
}

template <class T>
inline void requireNear(const T& actual, const T& expected, const T& tolerance, const std::string& label) {
    if (std::abs(actual - expected) <= tolerance) {
        return;
    }

    std::ostringstream oss;
    oss << label << ": expected `" << expected << "` +/- `" << tolerance << "` but got `" << actual << "`";
    throw std::runtime_error(oss.str());
}

inline ImageUInt8Ptr makeImage(int rows, int cols, std::initializer_list<uint8_t> values) {
    requireEqual(static_cast<int>(values.size()), rows * cols, "image buffer size");
    auto image = ImageUInt8::create(rows, cols);
    int index = 0;
    for (uint8_t value : values) {
        (*image)[index++] = value;
    }
    return image;
}

inline ImageUInt8Ptr makeComponentTreeFixture() {
    return makeImage(
        4,
        4,
        {
            3, 3, 2, 2,
            3, 4, 4, 2,
            1, 4, 5, 2,
            1, 1, 5, 0,
        }
    );
}

inline std::vector<NodeId> collectNodeIds(auto range) {
    std::vector<NodeId> ids;
    for (NodeId id : range) {
        ids.push_back(id);
    }
    return ids;
}

inline std::vector<NodeId> exportParentArray(const MorphologicalTree& tree) {
    std::vector<NodeId> parent(
        static_cast<size_t>(tree.getNumTotalProperParts() + tree.getNumInternalNodeSlots()),
        InvalidNode
    );

    for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        parent[static_cast<size_t>(properPart)] = tree.getSmallestComponent(properPart);
    }

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        parent[static_cast<size_t>(tree.getNumTotalProperParts() + nodeId)] = tree.getNodeParent(nodeId);
    }

    return parent;
}

inline std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy(const WeightedMorphologicalTree& tree) {
    return tree.exportHigraHierarchy();
}

inline int computeAreaAttribute(MorphologicalTree& tree, NodeId nodeId) {
    auto [attrNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(tree, AREA);
    return static_cast<int>(buffer[attrNames.linearIndex(nodeId, AREA)]);
}

template <class PixelType>
inline void requireImageShape(const std::shared_ptr<Image<PixelType>> &image, int rows, int cols) {
    require(static_cast<bool>(image), "image must not be null");
    requireEqual(image->getNumRows(), rows, "image rows");
    requireEqual(image->getNumCols(), cols, "image cols");
}

template <class PixelType>
inline std::vector<PixelType> collectImageValues(const std::shared_ptr<Image<PixelType>>& image) {
    std::vector<PixelType> values(image->getSize());
    for (int i = 0; i < image->getSize(); ++i) {
        values[static_cast<std::size_t>(i)] = (*image)[i];
    }
    return values;
}

} // namespace mmcfilters::unit_tests
