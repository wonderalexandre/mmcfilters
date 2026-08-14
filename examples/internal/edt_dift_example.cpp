/**
 * Internal visualization example for the exact distance-transform support used
 * by the `MAX_DIST` attribute family.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_internal_edt_dift tmp/edt_dift_example`.
 *
 * The program writes binary and distance-transform PNGs. It includes a
 * `detail/` header and is not a public API example.
 */
#include "mmcfilters/attributes/computers/detail/maxdist/EdtDIFT.hpp"
#include "stb_image_write.h"

#include <filesystem>
#include <iostream>
#include <string>

namespace {

void writePng(const std::filesystem::path& outputPath, const mmcfilters::ImageUInt8Ptr& image) {
    if (!stbi_write_png(outputPath.string().c_str(), image->getNumColumns(), image->getNumRows(), 1, image->rawData(), 0)) {
        throw std::runtime_error("Failed to write PNG file: " + outputPath.string());
    }
}

bool isContourPixel(int row, int column, int top, int left, int bottom, int right) {
    const bool inside = row >= top && row <= bottom && column >= left && column <= right;
    if (!inside) {
        return false;
    }
    return row == top || row == bottom || column == left || column == right;
}

} // namespace

int main(int argc, char** argv) {
    const std::filesystem::path outputDir = argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("tmp/edt_dift_example");

    std::filesystem::create_directories(outputDir);

    constexpr int rows = 15;
    constexpr int columns = 15;
    constexpr int top = 3;
    constexpr int left = 4;
    constexpr int bottom = 11;
    constexpr int right = 10;

    mmcfilters::attributes::computers::detail::maxdist::EdtDIFT edt(rows, columns);

    for (int row = top; row <= bottom; ++row) {
        for (int column = left; column <= right; ++column) {
            const int pidx = row * columns + column;
            edt.addPixelToBinaryImage(pidx);

            if (isContourPixel(row, column, top, left, bottom, right)) {
                edt.seed(pidx);
            } else {
                edt.open(pidx);
                edt.insertNeighborsPQueue(pidx);
            }
        }
    }

    edt.run();

    const auto binaryImage = edt.binaryImageForVisualisation();
    const auto distanceImage = edt.distanceTransformImage();

    const auto binaryOutput = outputDir / "binary.png";
    const auto distanceOutput = outputDir / "distance_transform.png";

    writePng(binaryOutput, binaryImage);
    writePng(distanceOutput, distanceImage);

    std::cout << "Saved EdtDIFT example outputs to:\n"
              << "  " << binaryOutput << "\n"
              << "  " << distanceOutput << "\n";

    return 0;
}
