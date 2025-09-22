#pragma once

#include "../include/Common.hpp"
#include "../include/AttributeNames.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerAttributeBasedBitQuads.hpp"
#include "../include/AttributeComputedIncrementally.hpp"

#define PI 3.14159265358979323846

/**
 * @brief Interface base para computadores de atributos associados a uma árvore morfológica.
 *
 * Define o contrato para classes que preenchem buffers de atributos, incluindo
 * metadados sobre atributos produzidos e dependências necessárias para o
 * cálculo incremental.
 */
class AttributeComputer {
	public:
		virtual ~AttributeComputer() = default;
	
		/// Executa a computação dos atributos produzidos por essa classe
		virtual void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const {
			compute(tree, buffer, attrNames, this->attributes(), dependencySources);
		}

		/// Executa a computação somente dos atributos solicitados
		virtual void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const = 0;
	
		/// Atributos produzidos
		virtual std::vector<Attribute> attributes() const = 0;
	
		/// Atributos requeridos para o cálculo (apenas metadado)
		virtual std::vector<AttributeOrGroup> requiredAttributes() const { return {}; }
	
};



/**
 * @brief Computa a área (número de pixels) de cada nó da árvore.
 */
class AreaComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override { return {AREA}; }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>&, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing AREA" << std::endl;
        auto indexOf = [&](NodeId idx) { return attrNames->linearIndex(idx, AREA); };
        for(NodeId id: tree->getNodeIds()){
            buffer[indexOf(id)] = static_cast<float>(tree->getAreaById(id));
        }
    }
};

/**
 * @brief Calcula volume e volume relativo acumulando níveis cinza sobre a árvore.
 */
class VolumeComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override { return {VOLUME, RELATIVE_VOLUME}; }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing VOLUME" << std::endl;
        auto indexOfVol = [&](NodeId idx) { return attrNames->linearIndex(idx, VOLUME); };
        auto indexOfRel = [&](NodeId idx) { return attrNames->linearIndex(idx, RELATIVE_VOLUME); };

        bool computeVolume = std::find(requestedAttributes.begin(), requestedAttributes.end(), VOLUME) != requestedAttributes.end();
        bool computeRelative = std::find(requestedAttributes.begin(), requestedAttributes.end(), RELATIVE_VOLUME) != requestedAttributes.end();

        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                if (computeVolume)
                    buffer[indexOfVol(node)] = static_cast<float>(tree->getNumCNPsById(node) * tree->getLevelById(node));
                if (computeRelative)
                    buffer[indexOfRel(node)] = 0.0f;
            },
            [&](NodeId parent, NodeId child) {
                if (computeVolume)
                    buffer[indexOfVol(parent)] += buffer[indexOfVol(child)];
                if (computeRelative)
                    buffer[indexOfRel(parent)] += buffer[indexOfRel(child)] + static_cast<float>(tree->getAreaById(child) * std::abs(tree->getLevelById(child) - tree->getLevelById(parent)));
            },
            [&](NodeId node) {
            if (computeRelative)
                buffer[indexOfRel(node)] += static_cast<float>(tree->getAreaById(node));
        });
    }
};

/**
 * @brief Calcula estatísticas básicas de níveis de cinza (média, variância, altura).
 */
class GrayLevelStatsComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT};
    }

    std::vector<AttributeOrGroup> requiredAttributes() const override {
        return {VOLUME};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing GrayLevelStatsComputer " << std::endl;

        auto indexOfMean = [&](NodeId idx) { return attrNames->linearIndex(idx, MEAN_LEVEL); };
        auto indexOfLevel = [&](NodeId idx) { return attrNames->linearIndex(idx, LEVEL); };
        auto indexOfVariance = [&](NodeId idx) { return attrNames->linearIndex(idx, VARIANCE_LEVEL); };
        auto indexOfGrayHeight = [&](NodeId idx) { return attrNames->linearIndex(idx, GRAY_HEIGHT); };

        bool computeMeanLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), MEAN_LEVEL) != requestedAttributes.end();
        bool computeVarianceLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), VARIANCE_LEVEL) != requestedAttributes.end();
        bool computeLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEVEL) != requestedAttributes.end();
        bool computeGrayHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), GRAY_HEIGHT) != requestedAttributes.end();

        auto [dependencyAttrNamesVol, bufferVol] = dependencySources[0];
        auto indexOfVol = [&](NodeId idx) { return dependencyAttrNamesVol->linearIndex(idx, VOLUME); };

        std::shared_ptr<long[]> sumGrayLevelSquare = nullptr;
        if (computeVarianceLevel) {
            sumGrayLevelSquare = std::shared_ptr<long[]>(new long[tree->getNumNodes()]);
        }

        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                if (computeVarianceLevel)
                    sumGrayLevelSquare[node] = static_cast<long>(tree->getNumCNPsById(node) * std::pow(tree->getLevelById(node), 2));
                if (computeLevel)
                    buffer[indexOfLevel(node)] = static_cast<float>(tree->getLevelById(node));
                if (computeGrayHeight)
                    buffer[indexOfGrayHeight(node)] = static_cast<float>(tree->getLevelById(node));
            },
            [&](NodeId parent, NodeId child) {
                if (computeVarianceLevel)
                    sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
                if (computeGrayHeight) {
                    float childValue = buffer[indexOfGrayHeight(child)];
                    float& parentValue = buffer[indexOfGrayHeight(parent)];
                    if (tree->isMaxtree() || tree->isMaxtreeNodeById(parent))
                        parentValue = std::max(parentValue, childValue);
                    else
                        parentValue = std::min(parentValue, childValue);
                }
            },
            [&](NodeId node) {
                float area = static_cast<float>(tree->getAreaById(node));
                if (computeMeanLevel)
                    buffer[indexOfMean(node)] = bufferVol[indexOfVol(node)] / area;
                if (computeVarianceLevel) {
                    float meanGrayLevel = bufferVol[indexOfVol(node)] / area; //mean graylevel = E(f)
                    double meanGrayLevelSquare = sumGrayLevelSquare[node] / area;  // E(f^2)
                    float var = static_cast<float>(meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel)); //variance: E(f^2) - E(f)^2
                    buffer[indexOfVariance(node)] = var > 0.0f ? var : 0.0f; //variance

                }
            }
        );

        if (computeGrayHeight) {
            for (NodeId node : tree->getIteratorPostOrderTraversalById()) {
                if (tree->isLeafById(node))
                    buffer[indexOfGrayHeight(node)] = 0.0f;
                else
                    buffer[indexOfGrayHeight(node)] = std::abs(tree->getLevelById(node) - buffer[indexOfGrayHeight(node)]) + 1.0f;
            }
        }
    }
};

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

/**
 * @brief Calcula momentos centrais geométricos até a terceira ordem.
 */
class CentralMomentsComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30, CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requested, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing CENTRAL_MOMENT group" << std::endl;

        int numCols = tree->getNumColsOfImage();
        int n = tree->getNumNodes();
        std::vector<long long> sumX(n, 0);
        std::vector<long long> sumY(n, 0);

        auto indexOf = [&](NodeId idx, Attribute attr) { return attrNames->linearIndex(idx, attr); };

        bool computeMu20 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_20) != requested.end();
        bool computeMu02 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_02) != requested.end();
        bool computeMu11 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_11) != requested.end();
        bool computeMu30 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_30) != requested.end();
        bool computeMu03 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_03) != requested.end();
        bool computeMu21 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_21) != requested.end();
        bool computeMu12 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_12) != requested.end();

        // Computa sumX e sumY para calcular os centroides
        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                sumX[node] = 0;
                sumY[node] = 0;
                for (int p : tree->getCNPsById(node)) {
                    auto [py, px] = ImageUtils::to2D(p, numCols);
                    sumX[node] += px;
                    sumY[node] += py;
                }
            },
            [&](NodeId parent, NodeId child) {
                sumX[parent] += sumX[child];
                sumY[parent] += sumY[child];
            },
            [](NodeId) {}
        );

        // Computação dos momentos centrais
        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                if (computeMu20) buffer[indexOf(node, CENTRAL_MOMENT_20)] = 0.0f;
                if (computeMu02) buffer[indexOf(node, CENTRAL_MOMENT_02)] = 0.0f;
                if (computeMu11) buffer[indexOf(node, CENTRAL_MOMENT_11)] = 0.0f;
                if (computeMu30) buffer[indexOf(node, CENTRAL_MOMENT_30)] = 0.0f;
                if (computeMu03) buffer[indexOf(node, CENTRAL_MOMENT_03)] = 0.0f;
                if (computeMu21) buffer[indexOf(node, CENTRAL_MOMENT_21)] = 0.0f;
                if (computeMu12) buffer[indexOf(node, CENTRAL_MOMENT_12)] = 0.0f;

                // Cálculo do centroide
                float area = static_cast<float>(tree->getAreaById(node));
                if (area <= 0.0f) return;
                float xCentroid = static_cast<float>(sumX[node]) / area;
                float yCentroid = static_cast<float>(sumY[node]) / area;

                for (int p : tree->getCNPsById(node)) {
                    auto [py, px] = ImageUtils::to2D(p, numCols);
                    float dx = px - xCentroid;
                    float dy = py - yCentroid;
                    if (computeMu20) buffer[indexOf(node, CENTRAL_MOMENT_20)] += dx * dx;
                    if (computeMu02) buffer[indexOf(node, CENTRAL_MOMENT_02)] += dy * dy;
                    if (computeMu11) buffer[indexOf(node, CENTRAL_MOMENT_11)] += dx * dy;
                    if (computeMu30) buffer[indexOf(node, CENTRAL_MOMENT_30)] += dx * dx * dx;
                    if (computeMu03) buffer[indexOf(node, CENTRAL_MOMENT_03)] += dy * dy * dy;
                    if (computeMu21) buffer[indexOf(node, CENTRAL_MOMENT_21)] += dx * dx * dy;
                    if (computeMu12) buffer[indexOf(node, CENTRAL_MOMENT_12)] += dx * dy * dy;
                }
            },
            [&](NodeId parent, NodeId child) {
                if (computeMu20) buffer[indexOf(parent, CENTRAL_MOMENT_20)] += buffer[indexOf(child, CENTRAL_MOMENT_20)];
                if (computeMu02) buffer[indexOf(parent, CENTRAL_MOMENT_02)] += buffer[indexOf(child, CENTRAL_MOMENT_02)];
                if (computeMu11) buffer[indexOf(parent, CENTRAL_MOMENT_11)] += buffer[indexOf(child, CENTRAL_MOMENT_11)];
                if (computeMu30) buffer[indexOf(parent, CENTRAL_MOMENT_30)] += buffer[indexOf(child, CENTRAL_MOMENT_30)];
                if (computeMu03) buffer[indexOf(parent, CENTRAL_MOMENT_03)] += buffer[indexOf(child, CENTRAL_MOMENT_03)];
                if (computeMu21) buffer[indexOf(parent, CENTRAL_MOMENT_21)] += buffer[indexOf(child, CENTRAL_MOMENT_21)];
                if (computeMu12) buffer[indexOf(parent, CENTRAL_MOMENT_12)] += buffer[indexOf(child, CENTRAL_MOMENT_12)];
            },
            [](NodeId) {}
        );
    }
};

/**
 * @brief Calcula atributos baseados em momentos (eixos principais, excentricidade, etc.).
 */
class MomentBasedAttributeComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {COMPACTNESS, ECCENTRICITY, LENGTH_MAJOR_AXIS, LENGTH_MINOR_AXIS,AXIS_ORIENTATION, INERTIA, CIRCULARITY};
    }

    std::vector<AttributeOrGroup> requiredAttributes() const override {
        return {AttributeGroup::CENTRAL_MOMENTS};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing MOMENT_BASED group" << std::endl;


        auto indexOfMajorAxis = [&](int idx) { return attrNames->linearIndex(idx, LENGTH_MAJOR_AXIS); };
        auto indexOfMinorAxis = [&](int idx) { return attrNames->linearIndex(idx, LENGTH_MINOR_AXIS); };
        auto indexOfEccentricity = [&](int idx) { return attrNames->linearIndex(idx, ECCENTRICITY); };
        auto indexOfCompactness = [&](int idx) { return attrNames->linearIndex(idx, COMPACTNESS); };
        auto indexOfAxisOrientation = [&](int idx) { return attrNames->linearIndex(idx, AXIS_ORIENTATION); };
        auto indexOfInertia = [&](int idx) { return attrNames->linearIndex(idx, INERTIA); };
        auto indexOfCircularity = [&](int idx) { return attrNames->linearIndex(idx, CIRCULARITY); };

        bool computeMajorAxis  = std::find(requestedAttributes.begin(), requestedAttributes.end(), LENGTH_MAJOR_AXIS)  != requestedAttributes.end();
        bool computeMinorAxis = std::find(requestedAttributes.begin(), requestedAttributes.end(), LENGTH_MINOR_AXIS) != requestedAttributes.end();
        bool computeEccentricity = std::find(requestedAttributes.begin(), requestedAttributes.end(), ECCENTRICITY) != requestedAttributes.end();
        bool computeCompactness = std::find(requestedAttributes.begin(), requestedAttributes.end(), COMPACTNESS) != requestedAttributes.end();
        bool computeAxisOrientation = std::find(requestedAttributes.begin(), requestedAttributes.end(), AXIS_ORIENTATION) != requestedAttributes.end();
        bool computeInertia = std::find(requestedAttributes.begin(), requestedAttributes.end(), INERTIA) != requestedAttributes.end();
        bool computeCircularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), CIRCULARITY) != requestedAttributes.end();

        auto [namesMom, bufMom] = dependencySources[0];
        auto indexMu20 = [&](int idx) { return namesMom->linearIndex(idx, CENTRAL_MOMENT_20); };
        auto indexMu02 = [&](int idx) { return namesMom->linearIndex(idx, CENTRAL_MOMENT_02); };
        auto indexMu11 = [&](int idx) { return namesMom->linearIndex(idx, CENTRAL_MOMENT_11); };
        
        
        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId) {},
            [&](NodeId, NodeId) {},
            [&](NodeId idx) {
                float mu20 = bufMom[indexMu20(idx)];
                float mu02 = bufMom[indexMu02(idx)];
                float mu11 = bufMom[indexMu11(idx)];
                float area = tree->getAreaById(idx);

                float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);
                discriminant = std::max(discriminant, 0.0f);
                float lambda1 = mu20 + mu02 + std::sqrt(discriminant);  // maior autovalor
                float lambda2 = mu20 + mu02 - std::sqrt(discriminant);  // menor autovalor

                if(computeMajorAxis){
                    if (area > 0.0f && lambda1 > 0.0f) {
                        buffer[indexOfMajorAxis(idx)] = std::sqrt((2.0f * lambda1) / area);
                    } else {
                        buffer[indexOfMajorAxis(idx)] = 0.0f;
                    }
                }
                if(computeMinorAxis){
                    if (area > 0.0f && lambda2 > 0.0f) {
                        buffer[indexOfMinorAxis(idx)] = std::sqrt((2.0f * lambda2) / area);
                    } else {
                        buffer[indexOfMinorAxis(idx)] = 0.0f;
                    }
                }
                if(computeEccentricity){	
                    if (std::abs(lambda2) > std::numeric_limits<float>::epsilon()) {
                        buffer[indexOfEccentricity(idx)] = lambda1 / lambda2;
                    } else {
                        buffer[indexOfEccentricity(idx)] = lambda1 / 0.1f; // fallback para evitar divisão por zero
                    }
                }
                if(computeCompactness){
                    float denom = mu20 + mu02;
                    if (denom > std::numeric_limits<float>::epsilon()) {
                        buffer[indexOfCompactness(idx)] = (1.0f / (2.0f * static_cast<float>(M_PI))) * (area / denom);
                    } else {
                        buffer[indexOfCompactness(idx)] = 0.0f;
                    }
                }
                if(computeAxisOrientation){
                    // Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
                    if (mu20 != mu02 || mu11 != 0) {
                        float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
                        float degrees = radians * (180.0 / M_PI); // Converter para graus
                        buffer[indexOfAxisOrientation(idx)] = std::fmod(std::abs(degrees), 360.0f); ; // Armazenar a orientação no intervalo [0, 360]
                    } else {
                        buffer[indexOfAxisOrientation(idx)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
                    }
                }
                if(computeInertia){
                    float normMu20 = mu20 / std::pow(area, 2.0f);
                    float normMu02 = mu02 / std::pow(area, 2.0f);
                    buffer[indexOfInertia(idx)] = normMu20 + normMu02;
                }
                if(computeCircularity){	
                    float circularity;
                    if (std::abs(lambda1) > std::numeric_limits<float>::epsilon()) {
                        buffer[indexOfCircularity(idx)] = lambda2 / lambda1;
                    } else {
                        buffer[indexOfCircularity(idx)] = 0.0f; // forma degenerada → circularidade indefinida
                    }
                }

                
            }
        );
    }
};


/**
 * @brief Calcula os sete momentos invariantes de Hu
 */
class HuMomentsComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {HU_MOMENT_1, HU_MOMENT_2, HU_MOMENT_3, HU_MOMENT_4, HU_MOMENT_5, HU_MOMENT_6, HU_MOMENT_7};
    }

    std::vector<AttributeOrGroup> requiredAttributes() const override {
        return {AttributeGroup::CENTRAL_MOMENTS};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requested, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing HU_MOMENT group" << std::endl;

        auto indexOf = [&](NodeId idx, Attribute attr) { return attrNames->linearIndex(idx, attr); };
        auto [dependencyAttrNamesMu, bufferMu] = dependencySources[0];
        auto indexOfMu = [&](NodeId idx, Attribute attr) { return dependencyAttrNamesMu->linearIndex(idx, attr); };

        auto normMoment = [](int area, float moment, int p, int q) {
            return moment / std::pow(static_cast<float>(area), (p + q + 2.0f) / 2.0f);
        };

        bool computeHu1 = std::find(requested.begin(), requested.end(), HU_MOMENT_1) != requested.end();
        bool computeHu2 = std::find(requested.begin(), requested.end(), HU_MOMENT_2) != requested.end();
        bool computeHu3 = std::find(requested.begin(), requested.end(), HU_MOMENT_3) != requested.end();
        bool computeHu4 = std::find(requested.begin(), requested.end(), HU_MOMENT_4) != requested.end();
        bool computeHu5 = std::find(requested.begin(), requested.end(), HU_MOMENT_5) != requested.end();
        bool computeHu6 = std::find(requested.begin(), requested.end(), HU_MOMENT_6) != requested.end();
        bool computeHu7 = std::find(requested.begin(), requested.end(), HU_MOMENT_7) != requested.end();

        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [](NodeId) {},
            [](NodeId, NodeId) {},
            [&](NodeId idx) {
                float mu20 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_20)];
                float mu02 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_02)];
                float mu11 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_11)];
                float mu30 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_30)];
                float mu03 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_03)];
                float mu21 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_21)];
                float mu12 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_12)];
                int area = tree->getAreaById(idx);

                // Calcular os momentos normalizados
                float eta20 = normMoment(area, mu20, 2, 0);
                float eta02 = normMoment(area, mu02, 0, 2);
                float eta11 = normMoment(area, mu11, 1, 1);
                float eta30 = normMoment(area, mu30, 3, 0);
                float eta03 = normMoment(area, mu03, 0, 3);
                float eta21 = normMoment(area, mu21, 2, 1);
                float eta12 = normMoment(area, mu12, 1, 2);

                // Cálculo dos momentos de Hu
                if(computeHu1)
                    buffer[indexOf(idx, HU_MOMENT_1)] = eta20 + eta02; // primeiro momento de Hu => inertia
                if(computeHu2)
                    buffer[indexOf(idx, HU_MOMENT_2)]  = std::pow(eta20 - eta02, 2) + 4 * std::pow(eta11, 2);
                if(computeHu3)
                    buffer[indexOf(idx, HU_MOMENT_3)]  = std::pow(eta30 - 3 * eta12, 2) + std::pow(3 * eta21 - eta03, 2);
                if(computeHu4)
                    buffer[indexOf(idx, HU_MOMENT_4)]  = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
                if(computeHu5)
                    buffer[indexOf(idx, HU_MOMENT_5)] = (eta30 - 3 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) +
                                                    (3 * eta21 - eta03) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
                if(computeHu6)
                    buffer[indexOf(idx, HU_MOMENT_6)] = (eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + 
                                                    4 * eta11 * (eta30 + eta12) * (eta21 + eta03);
                if(computeHu7)
                    buffer[indexOf(idx, HU_MOMENT_7)] = (3 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) -
                                                    (eta30 - 3 * eta12) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
            }
        );
    }
};

/**
 * @brief Calcula atributos estruturais relacionados à topologia da árvore.
 */
class TreeTopologyComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {HEIGHT_NODE, DEPTH_NODE, IS_LEAF_NODE, IS_ROOT_NODE, NUM_CHILDREN_NODE, NUM_SIBLINGS_NODE, NUM_DESCENDANTS_NODE, NUM_LEAF_DESCENDANTS_NODE, LEAF_RATIO_NODE, BALANCE_NODE, AVG_CHILD_HEIGHT_NODE};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing STRUCTURE_TREE group" << std::endl;

        bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), HEIGHT_NODE) != requestedAttributes.end();
        bool computeDepth = std::find(requestedAttributes.begin(), requestedAttributes.end(), DEPTH_NODE) != requestedAttributes.end();
        bool computeIsLeaf = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_LEAF_NODE) != requestedAttributes.end();
        bool computeIsRoot = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_ROOT_NODE) != requestedAttributes.end();
        bool computeNumChildren = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_CHILDREN_NODE) != requestedAttributes.end();
        bool computeNumSiblings = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_SIBLINGS_NODE) != requestedAttributes.end();
        bool computeNumDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeNumLeafDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_LEAF_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeLeafRatio = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEAF_RATIO_NODE) != requestedAttributes.end();
        bool computeBalance = std::find(requestedAttributes.begin(), requestedAttributes.end(), BALANCE_NODE) != requestedAttributes.end();
        bool computeAvgChildHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), AVG_CHILD_HEIGHT_NODE) != requestedAttributes.end();

        std::shared_ptr<float[]> bufferHeight = computeHeight ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfHeight = [&](NodeId idx) {
            return computeHeight ? attrNames->linearIndex(idx, HEIGHT_NODE) : idx;
        };

        std::shared_ptr<float[]> bufferNumDesc = computeNumDescendants ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfNumDescendants = [&](NodeId idx) {
            return computeNumDescendants ? attrNames->linearIndex(idx, NUM_DESCENDANTS_NODE) : idx;
        };

        std::shared_ptr<float[]> bufferNumLeafDesc = computeNumLeafDescendants ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfNumLeafDescendants = [&](NodeId idx) {
            return computeNumLeafDescendants ? attrNames->linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE) : idx;
        };

        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                NodeId parent = tree->getParentById(node);
                int parentDepth = parent != -1? static_cast<int>(bufferHeight[indexOfHeight(parent)]) : 0;
                bufferHeight[indexOfHeight(node)] = parent!=-1 ? parentDepth + 1.0f : 0.0f;

                bufferNumDesc[indexOfNumDescendants(node)] = 0.0f;
                bufferNumLeafDesc[indexOfNumLeafDescendants(node)] = tree->getNumChildrenById(node) == 0 ? 1.0f : 0.0f;

                if (computeHeight)
                    buffer[attrNames->linearIndex(node, HEIGHT_NODE)] = 0.0f;
                if (computeIsLeaf)
                    buffer[attrNames->linearIndex(node, IS_LEAF_NODE)] = tree->getNumChildrenById(node) == 0 ? 1.0f : 0.0f;
                if (computeIsRoot)
                    buffer[attrNames->linearIndex(node, IS_ROOT_NODE)] = parent!=-1 ? 0.0f : 1.0f;
                if (computeNumChildren)
                    buffer[attrNames->linearIndex(node, NUM_CHILDREN_NODE)] = static_cast<float>(tree->getNumChildrenById(node));
                if (computeNumSiblings)
                    buffer[attrNames->linearIndex(node, NUM_SIBLINGS_NODE)] = parent!=-1 ? static_cast<float>( tree->getNumChildrenById(parent) - 1) : 0.0f;
                if (computeLeafRatio)
                    buffer[attrNames->linearIndex(node, LEAF_RATIO_NODE)] = 0.0f;
                if (computeBalance)
                    buffer[attrNames->linearIndex(node, BALANCE_NODE)] = 0.0f;
                if (computeAvgChildHeight)
                    buffer[attrNames->linearIndex(node, AVG_CHILD_HEIGHT_NODE)] = 0.0f;
            },
            [&](NodeId parent, NodeId child) {
                
                bufferNumDesc[indexOfNumDescendants(parent)] += bufferNumDesc[indexOfNumDescendants(child)] + 1.0f;
                bufferNumLeafDesc[indexOfNumLeafDescendants(parent)] += bufferNumLeafDesc[indexOfNumLeafDescendants(child)];

                float childHeight = bufferHeight[indexOfHeight(child)];
                float& parentHeight = bufferHeight[indexOfHeight(parent)];
                parentHeight = std::max(parentHeight, childHeight + 1.0f);
                int numChildren = tree->getNumChildrenById(parent);

                if (computeBalance) {
                    float& minH = buffer[attrNames->linearIndex(parent, BALANCE_NODE)];
                    if (numChildren == 1)
                        minH = childHeight;
                    else
                        minH = std::min(minH, childHeight);
                }

                if (computeAvgChildHeight) {
                    float& sumH = buffer[attrNames->linearIndex(parent, AVG_CHILD_HEIGHT_NODE)];
                    if (numChildren == 1)
                        sumH = childHeight;
                    else
                        sumH += childHeight;
                }
            },
            [&](NodeId idx) {
                
                if (computeLeafRatio) {
                    float desc = bufferNumDesc[indexOfNumDescendants(idx)];
                    float folhas = bufferNumLeafDesc[indexOfNumLeafDescendants(idx)];
                    buffer[attrNames->linearIndex(idx, LEAF_RATIO_NODE)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
                }

                if (tree->getNumChildrenById(idx) > 0) {
                    if (computeBalance) {
                        float alturaMax = bufferHeight[indexOfHeight(idx)];
                        float alturaMin = buffer[attrNames->linearIndex(idx, BALANCE_NODE)];
                        buffer[attrNames->linearIndex(idx, BALANCE_NODE)] = alturaMax - alturaMin;
                    }

                    if (computeAvgChildHeight) {
                        buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] /
                                                                                     static_cast<float>(tree->getNumChildrenById(idx));
                    }
                }
            }
        );
    }
};




/**
 * @brief Calcula atributos derivados dos padrões Bit-Quads
 */
class BitquadsComputer : public AttributeComputer {
	public:
		
		std::vector<Attribute> attributes() const override {
			return {BITQUADS_AREA,
					BITQUADS_NUMBER_EULER,
					BITQUADS_NUMBER_HOLES,
					BITQUADS_PERIMETER,
					BITQUADS_PERIMETER_CONTINUOUS,
					BITQUADS_CIRCULARITY,
					BITQUADS_PERIMETER_AVERAGE,
					BITQUADS_LENGTH_AVERAGE,
					BITQUADS_WIDTH_AVERAGE};
		}

		void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing BITQUADS group" << std::endl;
			int numCols = tree->getNumColsOfImage();
			auto indexOf = [&](int idx, Attribute attr) {
				return attrNames->linearIndex(idx, attr);
			};
			bool computeArea = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_AREA) != requestedAttributes.end();
			bool computeNumberEuler = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_NUMBER_EULER) != requestedAttributes.end();
			bool computeNumberHoles = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_NUMBER_HOLES) != requestedAttributes.end();
			bool computePerimeter = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_PERIMETER) != requestedAttributes.end();
			bool computePerimeterCont = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_PERIMETER_CONTINUOUS) != requestedAttributes.end();
			bool computeCircularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_CIRCULARITY) != requestedAttributes.end();
			bool computePerimeterAverage = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_PERIMETER_AVERAGE) != requestedAttributes.end();
			bool computeLengthAverage = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_LENGTH_AVERAGE) != requestedAttributes.end();
			bool computeWithAverage = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_WIDTH_AVERAGE) != requestedAttributes.end();


			ComputerAttributeBasedBitQuads computerBitQuads(tree);
			std::vector<AttributeBasedBitQuads> attr = computerBitQuads.getAttributes();
			for(NodeId node: tree->getNodeIds()){
				if(computeArea)
					buffer[indexOf(node, BITQUADS_AREA)] = attr[node].getAreaDuda();
				if(computeNumberEuler)
					buffer[indexOf(node, BITQUADS_NUMBER_EULER)] = attr[node].getNumberEuler();
				if(computeNumberHoles)
					buffer[indexOf(node, BITQUADS_NUMBER_HOLES)] = attr[node].getNumberHoles();
				if(computePerimeter)
					buffer[indexOf(node, BITQUADS_PERIMETER)] = attr[node].getPerimeter();
				if(computePerimeterCont)
					buffer[indexOf(node, BITQUADS_PERIMETER_CONTINUOUS)] = attr[node].getPerimeterContinuous();
				if(computeCircularity)
					buffer[indexOf(node, BITQUADS_CIRCULARITY)] = attr[node].getCircularity();
				if(computePerimeterAverage)
					buffer[indexOf(node, BITQUADS_PERIMETER_AVERAGE)] = attr[node].getPerimeterAverage();
				if(computeLengthAverage)
					buffer[indexOf(node, BITQUADS_LENGTH_AVERAGE)] = attr[node].getLengthAverage();
				if(computeWithAverage)
					buffer[indexOf(node, BITQUADS_WIDTH_AVERAGE)] = attr[node].getWidthAverage();
			}
			
		}
};


/**
 * @brief Fábrica responsável por instanciar computadores de atributos sob demanda.
 */
class AttributeFactory {
	private:
		static std::shared_ptr<AttributeComputer> createImpl(Attribute attr) {
			switch (attr) {
				case AREA: return std::make_shared<AreaComputer>();
				
				case RELATIVE_VOLUME:
				case VOLUME: return std::make_shared<VolumeComputer>();
				
				case GRAY_HEIGHT: 
				case LEVEL: 
				case MEAN_LEVEL:
				case VARIANCE_LEVEL: return std::make_shared<GrayLevelStatsComputer>();

				case BOX_COL_MIN:
				case BOX_COL_MAX:
				case BOX_ROW_MIN:
				case BOX_ROW_MAX:
				case RATIO_WH: 
				case RECTANGULARITY: 
				case DIAGONAL_LENGTH:
				case BOX_HEIGHT:
				case BOX_WIDTH: 
					return std::make_shared<BoundingBoxComputer>();
				

				case AXIS_ORIENTATION: 
				case LENGTH_MAJOR_AXIS: 
				case LENGTH_MINOR_AXIS: 
				case ECCENTRICITY: 
				case INERTIA:
				case COMPACTNESS: 
				case CIRCULARITY:
					return std::make_shared<MomentBasedAttributeComputer>();


				case CENTRAL_MOMENT_20:
				case CENTRAL_MOMENT_02:
				case CENTRAL_MOMENT_11:
				case CENTRAL_MOMENT_30:
				case CENTRAL_MOMENT_03:
				case CENTRAL_MOMENT_21:
				case CENTRAL_MOMENT_12:
					return std::make_shared<CentralMomentsComputer>();

				
				case HU_MOMENT_1: 
				case HU_MOMENT_2:
				case HU_MOMENT_3:
				case HU_MOMENT_4:
				case HU_MOMENT_5:
				case HU_MOMENT_6:
				case HU_MOMENT_7:
					return std::make_shared<HuMomentsComputer>();


				case HEIGHT_NODE:
				case DEPTH_NODE:
				case IS_LEAF_NODE:
				case IS_ROOT_NODE:
				case NUM_CHILDREN_NODE:
				case NUM_SIBLINGS_NODE:
				case NUM_DESCENDANTS_NODE:
				case NUM_LEAF_DESCENDANTS_NODE:
				case LEAF_RATIO_NODE:
				case BALANCE_NODE:
				case AVG_CHILD_HEIGHT_NODE:
					return std::make_shared<TreeTopologyComputer>();
				
				case BITQUADS_AREA:
				case BITQUADS_NUMBER_EULER:
				case BITQUADS_NUMBER_HOLES:
				case BITQUADS_PERIMETER:
				case BITQUADS_PERIMETER_CONTINUOUS:
				case BITQUADS_CIRCULARITY:
				case BITQUADS_PERIMETER_AVERAGE:
				case BITQUADS_LENGTH_AVERAGE:
				case BITQUADS_WIDTH_AVERAGE:
					return std::make_shared<BitquadsComputer>();

				default:
					throw std::runtime_error("Attribute not supported.");
			}
		}

		static std::shared_ptr<AttributeComputer> createImpl(AttributeGroup group) {
			switch (group) {
				case AttributeGroup::BOUNDING_BOX:
					return std::make_shared<BoundingBoxComputer>();
				case AttributeGroup::CENTRAL_MOMENTS:
					return std::make_shared<CentralMomentsComputer>();
				case AttributeGroup::HU_MOMENTS:
					return std::make_shared<HuMomentsComputer>();
				case AttributeGroup::MOMENT_BASED:
					return std::make_shared<MomentBasedAttributeComputer>();
				case AttributeGroup::TREE_TOPOLOGY:
					return std::make_shared<TreeTopologyComputer>();
				case AttributeGroup::BITQUADS:
					return std::make_shared<BitquadsComputer>();
				default:
					throw std::runtime_error("Attribute group not supported.");
			}
		}

	public:
		static std::shared_ptr<AttributeComputer> create(const AttributeOrGroup& attr) {
			return std::visit([](auto&& actualAttr) -> std::shared_ptr<AttributeComputer> {
				return AttributeFactory::createImpl(actualAttr); // Correção aqui!
			}, attr);
		}



};

