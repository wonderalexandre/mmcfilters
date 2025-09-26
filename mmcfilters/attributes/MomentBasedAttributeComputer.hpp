#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {



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
                        buffer[indexOfCompactness(idx)] = (1.0f / (2.0f * std::numbers::pi) ) * (area / denom);
                    } else {
                        buffer[indexOfCompactness(idx)] = 0.0f;
                    }
                }
                if(computeAxisOrientation){
                    // Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
                    if (mu20 != mu02 || mu11 != 0) {
                        float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
                        float degrees = radians * (180.0 / std::numbers::pi); // Converter para graus
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





} // namespace mmcfilters

