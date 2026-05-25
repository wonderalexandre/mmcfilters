#pragma once

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/utils/Image.hpp"
#include "mmcfilters/utils/Common.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
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

template <class Exception = std::exception, class Fn>
inline void requireThrows(Fn&& fn, const std::string& label) {
    try {
        fn();
    } catch (const Exception&) {
        return;
    } catch (const std::exception& ex) {
        std::ostringstream oss;
        oss << label << ": expected requested exception type but got `" << ex.what() << "`";
        throw std::runtime_error(oss.str());
    }

    throw std::runtime_error(label + ": expected an exception");
}

template <class Exception = std::exception, class Fn>
inline void requireThrowsContaining(Fn&& fn, const std::string& expectedMessageFragment, const std::string& label) {
    try {
        fn();
    } catch (const Exception& ex) {
        const std::string message = ex.what();
        if (message.find(expectedMessageFragment) != std::string::npos) {
            return;
        }

        std::ostringstream oss;
        oss << label << ": expected exception message containing `" << expectedMessageFragment
            << "` but got `" << message << "`";
        throw std::runtime_error(oss.str());
    } catch (const std::exception& ex) {
        std::ostringstream oss;
        oss << label << ": expected requested exception type but got `" << ex.what() << "`";
        throw std::runtime_error(oss.str());
    }

    throw std::runtime_error(label + ": expected an exception");
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

inline std::shared_ptr<MorphologicalTree> makeComponentTree(ImageUInt8Ptr image, bool isMaxtree, double radius = 1.5) {
    auto weighted = isMaxtree
        ? MorphologicalTreeFactory::createMaxTree(image, radius)
        : MorphologicalTreeFactory::createMinTree(image, radius);
    return std::make_shared<MorphologicalTree>(weighted.topology().clone());
}

inline std::shared_ptr<MorphologicalTree> makeTreeOfShapes(
    ImageUInt8Ptr image,
    ToSInterpolation interpolation = ToSInterpolation::SelfDual) {
    auto weighted = MorphologicalTreeFactory::createTreeOfShapes(image, interpolation);
    return std::make_shared<MorphologicalTree>(weighted.topology().clone());
}

inline std::shared_ptr<MorphologicalTree> makeTreeFromHigraParent(
    const std::vector<NodeId>& parent,
    int rows,
    int cols,
    bool isMaxtree,
    double radius = 1.5) {
    std::vector<std::uint8_t> altitude(parent.size(), std::uint8_t{});
    auto weighted = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(parent),
        std::span<const std::uint8_t>(altitude),
        rows,
        cols,
        isMaxtree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE,
        AdjacencyRelation(rows, cols, radius));
    return std::make_shared<MorphologicalTree>(
        weighted.topology().clone());
}

inline std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> makeWeightedComponentTree(ImageUInt8Ptr image, bool isMaxtree, double radius = 1.5) {
    return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
        isMaxtree
            ? MorphologicalTreeFactory::createMaxTree(image, radius)
            : MorphologicalTreeFactory::createMinTree(image, radius));
}

inline std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> makeWeightedTreeOfShapes(
    ImageUInt8Ptr image,
    ToSInterpolation interpolation = ToSInterpolation::SelfDual) {
    return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
        MorphologicalTreeFactory::createTreeOfShapes(image, interpolation));
}

inline std::vector<NodeId> collectNodeIds(auto range) {
    std::vector<NodeId> ids;
    for (NodeId id : range) {
        ids.push_back(id);
    }
    return ids;
}

inline std::pair<std::vector<NodeId>, std::vector<std::uint8_t>> exportFlatHigraHierarchy(const MorphologicalTree& tree) {
    AltitudeBuffer<std::uint8_t> altitude(static_cast<size_t>(tree.getNumInternalNodeSlots()), std::uint8_t{});
    return TreeAltitudeAlgorithms::exportHigraHierarchy(tree, std::span<const std::uint8_t>(altitude));
}

inline std::pair<std::vector<NodeId>, std::vector<std::uint8_t>> exportHigraHierarchy(const WeightedMorphologicalTree<std::uint8_t>& tree) {
    return tree.exportHigraHierarchy();
}

inline int computeAreaViaAttributeFacade(const MorphologicalTree& tree, NodeId nodeId) {
    auto [attrNames, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, AREA);
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
