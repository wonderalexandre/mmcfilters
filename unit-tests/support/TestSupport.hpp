#pragma once

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/utils/Image.hpp"
#include "mmcfilters/utils/Common.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTree.hpp"

namespace mmcfilters::unit_tests {

inline void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class T, class U> inline void requireEqual(const T& actual, const U& expected, const std::string& label) {
    if (!(actual == expected)) {
        std::ostringstream oss;
        oss << label << ": expected `" << expected << "` but got `" << actual << "`";
        throw std::runtime_error(oss.str());
    }
}

template <class T> inline void requireVectorEqual(const std::vector<T>& actual, const std::vector<T>& expected, const std::string& label) {
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

template <class T> inline void requireNear(const T& actual, const T& expected, const T& tolerance, const std::string& label) {
    if (std::abs(actual - expected) <= tolerance) {
        return;
    }

    std::ostringstream oss;
    oss << label << ": expected `" << expected << "` +/- `" << tolerance << "` but got `" << actual << "`";
    throw std::runtime_error(oss.str());
}

template <class Exception = std::exception, class Fn> inline void requireThrows(Fn&& fn, const std::string& label) {
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
        oss << label << ": expected exception message containing `" << expectedMessageFragment << "` but got `" << message << "`";
        throw std::runtime_error(oss.str());
    } catch (const std::exception& ex) {
        std::ostringstream oss;
        oss << label << ": expected requested exception type but got `" << ex.what() << "`";
        throw std::runtime_error(oss.str());
    }

    throw std::runtime_error(label + ": expected an exception");
}

inline ImageUInt8Ptr makeImage(int rows, int columns, std::initializer_list<uint8_t> values) {
    requireEqual(static_cast<int>(values.size()), rows * columns, "image buffer size");
    auto image = ImageUInt8::create(rows, columns);
    int index = 0;
    for (uint8_t value : values) {
        (*image)[index++] = value;
    }
    return image;
}

inline ImageUInt8Ptr makeComponentTreeFixture() {
    return makeImage(4, 4,
                     {
                         3,
                         3,
                         2,
                         2,
                         3,
                         4,
                         4,
                         2,
                         1,
                         4,
                         5,
                         2,
                         1,
                         1,
                         5,
                         0,
                     });
}

inline std::shared_ptr<MorphologicalTree> makeComponentTree(ImageUInt8Ptr image, bool isMaxtree, double radius = 1.5) {
    auto valuedTree = isMaxtree ? MorphologicalTreeFactory::createMaxTree(image, radius) : MorphologicalTreeFactory::createMinTree(image, radius);
    return std::make_shared<MorphologicalTree>(valuedTree.topology().clone());
}

enum class TestTopographicImmersion { SelfDualSpan, Min4Max8, Min8Max4 };

// Test conventions publish exact doubled units so that a single helper covers
// every immersion, including the self-dual span. Helpers dedicated to the 8-bit
// encoding are defined below.
inline TopographicConvention makeTopographicConvention(int rows, int columns,
                                                        TestTopographicImmersion immersion = TestTopographicImmersion::SelfDualSpan,
                                                        TopographicDomainExtension domainExtension = TopographicDomainExtension::ExteriorRing,
                                                        PixelId infinityPixel = 0,
                                                        TopographicAltitudeEncoding altitudeEncoding = TopographicAltitudeEncoding::ExactDoubled) {
    TreeOfShapesImmersion specification = SelfDualSpanImmersion{};
    if (immersion == TestTopographicImmersion::Min4Max8) {
        specification = ComplementaryGridImmersion{
            ComplementaryAdjacencies{RegularGridAdjacency2D(rows, columns, 1.0), RegularGridAdjacency2D(rows, columns, 1.5)}};
    } else if (immersion == TestTopographicImmersion::Min8Max4) {
        specification = ComplementaryGridImmersion{
            ComplementaryAdjacencies{RegularGridAdjacency2D(rows, columns, 1.5), RegularGridAdjacency2D(rows, columns, 1.0)}};
    }
    return TopographicConvention{std::move(specification), domainExtension, infinityPixel, altitudeEncoding};
}

inline TopographicConvention makeTopographicConvention(const ImageUInt8Ptr& image,
                                                        TestTopographicImmersion immersion = TestTopographicImmersion::SelfDualSpan,
                                                        TopographicDomainExtension domainExtension = TopographicDomainExtension::ExteriorRing,
                                                        PixelId infinityPixel = 0,
                                                        TopographicAltitudeEncoding altitudeEncoding = TopographicAltitudeEncoding::ExactDoubled) {
    return makeTopographicConvention(image->getNumRows(), image->getNumColumns(), immersion, domainExtension, infinityPixel, altitudeEncoding);
}

inline std::shared_ptr<MorphologicalTree> makeTreeOfShapes(ImageUInt8Ptr image,
                                                           TestTopographicImmersion immersion = TestTopographicImmersion::SelfDualSpan) {
    auto valuedTree = MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(image, makeTopographicConvention(image, immersion));
    return std::make_shared<MorphologicalTree>(valuedTree.topology().clone());
}

template <class Altitude>
inline std::vector<Altitude> makeStrictHigraAltitude(const std::vector<NodeId>& parent, int numPixels, bool increasingFromRoot) {
    static_assert(std::is_integral_v<Altitude>);
    require(numPixels > 0 && numPixels < static_cast<NodeId>(parent.size()), "strict Higra altitude requires leaves followed by internal nodes");

    const auto unknown = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> depth(parent.size(), unknown);
    std::vector<NodeId> path;
    std::size_t maxDepth = 0;

    for (NodeId higraNode = numPixels; higraNode < static_cast<NodeId>(parent.size()); ++higraNode) {
        NodeId cursor = higraNode;
        path.clear();
        while (depth[static_cast<std::size_t>(cursor)] == unknown) {
            const NodeId parentNode = parent[static_cast<std::size_t>(cursor)];
            require(parentNode >= numPixels && parentNode < static_cast<NodeId>(parent.size()), "strict Higra altitude requires internal parents");
            if (parentNode == cursor) {
                depth[static_cast<std::size_t>(cursor)] = 0;
                break;
            }
            require(path.size() < parent.size(), "strict Higra altitude requires an acyclic hierarchy");
            path.push_back(cursor);
            cursor = parentNode;
        }

        std::size_t cursorDepth = depth[static_cast<std::size_t>(cursor)];
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            ++cursorDepth;
            depth[static_cast<std::size_t>(*it)] = cursorDepth;
            maxDepth = std::max(maxDepth, cursorDepth);
        }
    }

    require(maxDepth <= static_cast<std::size_t>(std::numeric_limits<Altitude>::max()), "strict Higra altitude exceeds the requested value type");

    std::vector<Altitude> altitude(parent.size(), Altitude{});
    for (NodeId higraNode = numPixels; higraNode < static_cast<NodeId>(parent.size()); ++higraNode) {
        const std::size_t nodeDepth = depth[static_cast<std::size_t>(higraNode)];
        altitude[static_cast<std::size_t>(higraNode)] = static_cast<Altitude>(increasingFromRoot ? nodeDepth : maxDepth - nodeDepth);
    }
    for (PixelId pixel = 0; pixel < numPixels; ++pixel) {
        const NodeId smallestNodeId = parent[static_cast<std::size_t>(pixel)];
        altitude[static_cast<std::size_t>(pixel)] = altitude[static_cast<std::size_t>(smallestNodeId)];
    }
    return altitude;
}

inline std::shared_ptr<MorphologicalTree> makeTreeFromHigraParent(const std::vector<NodeId>& parent, int rows, int columns, bool isMaxtree, double radius = 1.5) {
    const auto altitude = makeStrictHigraAltitude<std::uint8_t>(parent, rows * columns, isMaxtree);
    auto valuedTree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), rows, columns,
                                                                    isMaxtree ? MorphologicalTreeKind::MaxTree : MorphologicalTreeKind::MinTree,
                                                                    RegularGridAdjacency2D(rows, columns, radius));
    return std::make_shared<MorphologicalTree>(valuedTree.topology().clone());
}

inline std::shared_ptr<ValuedMorphologicalTree<std::uint8_t>> makeValuedComponentTree(ImageUInt8Ptr image, bool isMaxtree, double radius = 1.5) {
    return std::make_shared<ValuedMorphologicalTree<std::uint8_t>>(isMaxtree ? MorphologicalTreeFactory::createMaxTree(image, radius)
                                                                               : MorphologicalTreeFactory::createMinTree(image, radius));
}

inline std::shared_ptr<ValuedMorphologicalTree<ToSGrayLevel>> makeValuedTreeOfShapes(ImageUInt8Ptr image,
                                                                                        TestTopographicImmersion immersion = TestTopographicImmersion::SelfDualSpan) {
    return std::make_shared<ValuedMorphologicalTree<ToSGrayLevel>>(
        MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(image, makeTopographicConvention(image, immersion)));
}

// Builds the 8-bit tree of shapes published by the default convention.
inline std::shared_ptr<ValuedMorphologicalTree<std::uint8_t>>
makeValuedTreeOfShapesUInt8(ImageUInt8Ptr image, ComplementaryPairing pairing = ComplementaryPairing::Min4Max8) {
    TopographicConvention convention{CanonicalComplementaryGridImmersion{pairing}, TopographicDomainExtension::None, PixelId{0},
                                     TopographicAltitudeEncoding::UInt8};
    return std::make_shared<ValuedMorphologicalTree<std::uint8_t>>(
        MorphologicalTreeFactory::createTreeOfShapes<std::uint8_t>(image, std::move(convention)));
}

inline std::vector<NodeId> collectNodeIds(auto range) {
    std::vector<NodeId> ids;
    for (NodeId id : range) {
        ids.push_back(id);
    }
    return ids;
}

inline std::vector<PixelId> collectPixelIds(auto range) {
    std::vector<PixelId> ids;
    for (PixelId id : range) {
        ids.push_back(id);
    }
    return ids;
}

inline std::pair<std::vector<NodeId>, std::vector<std::uint8_t>> exportFlatHigraHierarchy(const MorphologicalTree& tree) {
    NodeAltitudeBuffer<std::uint8_t> altitude(static_cast<size_t>(tree.numInternalNodeSlots()), std::uint8_t{});
    return TreeAltitudeAlgorithms::exportHigraHierarchy(tree, std::span<const std::uint8_t>(altitude));
}

inline std::pair<std::vector<NodeId>, std::vector<std::uint8_t>> exportHigraHierarchy(const ValuedMorphologicalTree<std::uint8_t>& tree) {
    return tree.exportHigraHierarchy();
}

inline int computeAreaViaAttributeFacade(const MorphologicalTree& tree, NodeId nodeId) {
    auto [attrNames, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, Area);
    return static_cast<int>(buffer[attrNames.linearIndex(nodeId, Area)]);
}

template <class PixelType> inline void requireImageShape(const std::shared_ptr<Image<PixelType>>& image, int rows, int columns) {
    require(static_cast<bool>(image), "image must not be null");
    requireEqual(image->getNumRows(), rows, "image rows");
    requireEqual(image->getNumColumns(), columns, "image columns");
}

template <class PixelType> inline std::vector<PixelType> collectImageValues(const std::shared_ptr<Image<PixelType>>& image) {
    std::vector<PixelType> values(image->getSize());
    for (int i = 0; i < image->getSize(); ++i) {
        values[static_cast<std::size_t>(i)] = (*image)[i];
    }
    return values;
}

} // namespace mmcfilters::unit_tests
