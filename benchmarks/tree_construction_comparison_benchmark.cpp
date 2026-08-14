#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "stb_image.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>

using namespace mmcfilters;

namespace {

using Clock = std::chrono::steady_clock;
using UInt8Tree = ValuedMorphologicalTree<std::uint8_t>;
using TreeOfShapesTree = ValuedMorphologicalTree<ToSGrayLevel>;
using Tree = std::variant<UInt8Tree, TreeOfShapesTree>;

enum class Algorithm : std::size_t { TreeOfShapesMax4cMin8c, TreeOfShapesSelfDual, MaxTree8c, MinTree8c, UnrestrictedResidualTree8c, SaturatedResidualTree8c };

constexpr std::array<Algorithm, 6> Algorithms{
    Algorithm::TreeOfShapesMax4cMin8c,     Algorithm::TreeOfShapesSelfDual,   Algorithm::MaxTree8c, Algorithm::MinTree8c,
    Algorithm::UnrestrictedResidualTree8c, Algorithm::SaturatedResidualTree8c};

constexpr std::array<std::string_view, 3> Resolutions{"480p", "720p", "1080p"};

struct Options {
    std::filesystem::path dataRoot;
    std::filesystem::path outputPath;
    int repetitions = 5;
};

struct Measurement {
    double milliseconds = 0.0;
    int nodes = 0;
};

[[nodiscard]] Options parseOptions(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        throw std::invalid_argument("usage: tree_construction_comparison_benchmark "
                                    "<ICDAR data root> <output.csv> [repetitions]");
    }
    Options options;
    options.dataRoot = argv[1];
    options.outputPath = argv[2];
    if (argc == 4) {
        options.repetitions = std::stoi(argv[3]);
    }
    if (options.repetitions <= 0) {
        throw std::invalid_argument("repetitions must be positive");
    }
    return options;
}

[[nodiscard]] ImageUInt8Ptr loadImage(const std::filesystem::path& path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    using Buffer = std::unique_ptr<unsigned char, decltype(&stbi_image_free)>;
    Buffer pixels(stbi_load(path.string().c_str(), &width, &height, &channels, 1), &stbi_image_free);
    if (!pixels || width <= 0 || height <= 0) {
        throw std::runtime_error("could not load image: " + path.string());
    }
    auto image = ImageUInt8::create(height, width);
    std::copy_n(pixels.get(), image->getSize(), image->rawData());
    return image;
}

[[nodiscard]] std::string_view algorithmName(Algorithm algorithm) {
    switch (algorithm) {
    case Algorithm::TreeOfShapesMax4cMin8c:
        return "tos_max4c_min8c";
    case Algorithm::TreeOfShapesSelfDual:
        return "tos_self_dual";
    case Algorithm::MaxTree8c:
        return "max_tree_8c";
    case Algorithm::MinTree8c:
        return "min_tree_8c";
    case Algorithm::UnrestrictedResidualTree8c:
        return "residual_unrestricted_8c";
    case Algorithm::SaturatedResidualTree8c:
        return "residual_saturated_8c";
    }
    throw std::invalid_argument("unsupported benchmark algorithm");
}

[[nodiscard]] Tree buildTree(Algorithm algorithm, const ImageUInt8Ptr& image, const RegularGridAdjacency2D& adjacency8c) {
    switch (algorithm) {
    case Algorithm::TreeOfShapesMax4cMin8c:
        return MorphologicalTreeFactory::createTreeOfShapes(
            image, TopographicConvention{ComplementaryGridImmersion{ComplementaryAdjacencies{
                       RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.5),
                       RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.0)}}});
    case Algorithm::TreeOfShapesSelfDual:
        return MorphologicalTreeFactory::createTreeOfShapes(image);
    case Algorithm::MaxTree8c:
        return MorphologicalTreeFactory::createMaxTree(image, adjacency8c);
    case Algorithm::MinTree8c:
        return MorphologicalTreeFactory::createMinTree(image, adjacency8c);
    case Algorithm::UnrestrictedResidualTree8c:
        return MorphologicalTreeFactory::createUnrestrictedResidualTree(
            image, adjacency8c, sdrt::UnrestrictedResidualTreeOptions{});
    case Algorithm::SaturatedResidualTree8c:
        return MorphologicalTreeFactory::createSaturatedResidualTree(
            image, adjacency8c, PixelId{0}, sdrt::SaturatedResidualTreeOptions{});
    }
    throw std::invalid_argument("unsupported benchmark algorithm");
}

[[nodiscard]] bool reconstructsExactly(const Tree& tree, const ImageUInt8Ptr& image) {
    return std::visit(
        [&](const auto& concreteTree) {
            const auto reconstruction = concreteTree.reconstructFromNodeAltitudes();
            if (reconstruction->getNumRows() != image->getNumRows() || reconstruction->getNumColumns() != image->getNumColumns() ||
                reconstruction->getSize() != image->getSize()) {
                return false;
            }
            using Altitude = typename std::remove_cvref_t<decltype(concreteTree)>::AltitudeType;
            for (int index = 0; index < image->getSize(); ++index) {
                const auto expected = [&] {
                    if constexpr (std::is_same_v<Altitude, ToSGrayLevel>) {
                        return static_cast<ToSGrayLevel>(2u * (*image)[index]);
                    } else {
                        return static_cast<Altitude>((*image)[index]);
                    }
                }();
                if ((*reconstruction)[index] != expected) {
                    return false;
                }
            }
            return true;
        },
        tree);
}

[[nodiscard]] Measurement measureConstruction(Algorithm algorithm, const ImageUInt8Ptr& image, const RegularGridAdjacency2D& adjacency8c) {
    const auto start = Clock::now();
    Tree tree = buildTree(algorithm, image, adjacency8c);
    const auto stop = Clock::now();
    const int nodeCount = std::visit([](const auto& concreteTree) { return concreteTree.topology().numNodes(); }, tree);
    return {std::chrono::duration<double, std::milli>(stop - start).count(), nodeCount};
}

[[nodiscard]] std::filesystem::path imagePath(const std::filesystem::path& root, std::string_view resolution, int imageIndex) {
    std::ostringstream filename;
    filename << "test_" << std::setw(3) << std::setfill('0') << imageIndex << ".png";
    return root / ("icdar_" + std::string(resolution)) / filename.str();
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::ofstream output(options.outputPath);
        if (!output) {
            throw std::runtime_error("could not open output: " + options.outputPath.string());
        }
        output << "resolution,image,rows,columns,pixels,repetition,position,"
                  "algorithm,construction_ms,nodes\n"
               << std::fixed << std::setprecision(6);

        for (std::size_t resolutionIndex = 0; resolutionIndex < Resolutions.size(); ++resolutionIndex) {
            const std::string_view resolution = Resolutions[resolutionIndex];
            for (int imageIndex = 0; imageIndex < 10; ++imageIndex) {
                const std::filesystem::path path = imagePath(options.dataRoot, resolution, imageIndex);
                const auto image = loadImage(path);
                const RegularGridAdjacency2D adjacency8c(image->getNumRows(), image->getNumColumns(), 1.5);

                std::array<int, Algorithms.size()> expectedNodes{};
                for (Algorithm algorithm : Algorithms) {
                    Tree warmup = buildTree(algorithm, image, adjacency8c);
                    if (!reconstructsExactly(warmup, image)) {
                        throw std::runtime_error(std::string(algorithmName(algorithm)) + " failed exact reconstruction for " + path.string());
                    }
                    expectedNodes[static_cast<std::size_t>(algorithm)] =
                        std::visit([](const auto& concreteTree) { return concreteTree.topology().numNodes(); }, warmup);
                }

                for (int repetition = 0; repetition < options.repetitions; ++repetition) {
                    const std::size_t orderOffset =
                        (resolutionIndex + static_cast<std::size_t>(imageIndex) + static_cast<std::size_t>(repetition)) % Algorithms.size();
                    for (std::size_t position = 0; position < Algorithms.size(); ++position) {
                        const Algorithm algorithm = Algorithms[(orderOffset + position) % Algorithms.size()];
                        const Measurement measurement = measureConstruction(algorithm, image, adjacency8c);
                        if (measurement.nodes != expectedNodes[static_cast<std::size_t>(algorithm)]) {
                            throw std::runtime_error(std::string(algorithmName(algorithm)) + " produced a non-deterministic node count for " + path.string());
                        }
                        output << resolution << ',' << path.filename().string() << ',' << image->getNumRows() << ',' << image->getNumColumns() << ','
                               << image->getSize() << ',' << repetition << ',' << position << ',' << algorithmName(algorithm) << ',' << measurement.milliseconds
                               << ',' << measurement.nodes << '\n';
                    }
                }
                output.flush();
                std::cerr << "completed " << resolution << '/' << path.filename().string() << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "tree construction comparison benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
