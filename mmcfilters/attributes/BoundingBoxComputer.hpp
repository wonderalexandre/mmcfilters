#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"

namespace mmcfilters {

/**
 * @brief Calcula atributos baseados em caixa delimitadora (bounding box).
 */
class BoundingBoxComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {BOX_WIDTH, BOX_HEIGHT, RECTANGULARITY, RATIO_WH,BOX_COL_MIN, BOX_COL_MAX, BOX_ROW_MIN, BOX_ROW_MAX,DIAGONAL_LENGTH};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing BOUNDING_BOX group" << std::endl;

        auto indexOfWidth  = [&](NodeId idx) { return attrNames->linearIndex(idx, BOX_WIDTH); };
        auto indexOfHeight = [&](NodeId idx) { return attrNames->linearIndex(idx, BOX_HEIGHT); };
        auto indexOfRectangularity = [&](NodeId idx) { return attrNames->linearIndex(idx, RECTANGULARITY); };
        auto indexOfRatioWH = [&](NodeId idx) { return attrNames->linearIndex(idx, RATIO_WH); };
        auto indexOfColMin = [&](NodeId idx) { return attrNames->linearIndex(idx, BOX_COL_MIN); };
        auto indexOfColMax = [&](NodeId idx) { return attrNames->linearIndex(idx, BOX_COL_MAX); };
        auto indexOfRowMin = [&](NodeId idx) { return attrNames->linearIndex(idx, BOX_ROW_MIN); };
        auto indexOfRowMax = [&](NodeId idx) { return attrNames->linearIndex(idx, BOX_ROW_MAX); };
        auto indexOfDiagonalLength = [&](NodeId idx) { return attrNames->linearIndex(idx, DIAGONAL_LENGTH); };

        bool computeWidth  = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_WIDTH)  != requestedAttributes.end();
        bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_HEIGHT) != requestedAttributes.end();
        bool computeRectangularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), RECTANGULARITY) != requestedAttributes.end();
        bool computeRatioWH = std::find(requestedAttributes.begin(), requestedAttributes.end(), RATIO_WH) != requestedAttributes.end();
        bool computeColMin = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_COL_MIN) != requestedAttributes.end();
        bool computeColMax = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_COL_MAX) != requestedAttributes.end();
        bool computeRowMin = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_ROW_MIN) != requestedAttributes.end();
        bool computeRowMax = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_ROW_MAX) != requestedAttributes.end();
        bool computeDiagonalLength = std::find(requestedAttributes.begin(), requestedAttributes.end(), DIAGONAL_LENGTH) != requestedAttributes.end();

        int n = tree->getNumNodes();
        int numCols = tree->getNumColsOfImage();
        int numRows = tree->getNumRowsOfImage();

        std::vector<int> xmin(n, numCols);
        std::vector<int> xmax(n, 0);
        std::vector<int> ymin(n, numRows);
        std::vector<int> ymax(n, 0);

        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId idx) {
                xmin[idx] = numCols;
                xmax[idx] = 0;
                ymin[idx] = numRows;
                ymax[idx] = 0;

                for (int p : tree->getCNPsById(idx)) {
                    auto [y, x] = ImageUtils::to2D(p, numCols);
                    xmin[idx] = std::min(xmin[idx], x);
                    xmax[idx] = std::max(xmax[idx], x);
                    ymin[idx] = std::min(ymin[idx], y);
                    ymax[idx] = std::max(ymax[idx], y);
                }
            },
            [&](NodeId pid, NodeId cid) {
                xmin[pid] = std::min(xmin[pid], xmin[cid]);
                xmax[pid] = std::max(xmax[pid], xmax[cid]);
                ymin[pid] = std::min(ymin[pid], ymin[cid]);
                ymax[pid] = std::max(ymax[pid], ymax[cid]);
            },
            [&](NodeId idx) {
                if(computeWidth)
                    buffer[indexOfWidth(idx)]  = xmax[idx] - xmin[idx] + 1;
                if(computeHeight)
                    buffer[indexOfHeight(idx)] = ymax[idx] - ymin[idx] + 1;

                if(computeRectangularity) {
                    float area = tree->getAreaById(idx);
                    float width = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    float denom = width * height;
                    buffer[indexOfRectangularity(idx)] = (denom > 0.0f) ? (area / denom) : 0.0f;
                }
                if(computeRatioWH) {
                    float width  = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    if (width > 0 && height > 0) {
                        buffer[indexOfRatioWH(idx)] = std::max(width, height) / std::min(width, height);
                    } else {
                        buffer[indexOfRatioWH(idx)] = 0.0f;
                    }
                }
                if(computeColMin)
                    buffer[indexOfColMin(idx)]  = xmin[idx];
                if(computeColMax)
                    buffer[indexOfColMax(idx)]  = xmax[idx];
                if(computeRowMin)
                    buffer[indexOfRowMin(idx)]  = ymin[idx];
                if(computeRowMax)
                    buffer[indexOfRowMax(idx)]  = ymax[idx];
                if(computeDiagonalLength) {
                    float width  = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    buffer[indexOfDiagonalLength(idx)] = std::sqrt(width*width + height*height);
                }
            }
        );
    }
};

} // namespace mmcfilters

