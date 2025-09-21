#pragma once

#include <vector>
#include <set>
#include <memory>
#include <functional>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdint>
#include <cassert>
#include "../include/Common.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AttributeComputedIncrementallyCT.hpp"

//---------------------------------------------
// CLASSES QuadBit e padrões
//---------------------------------------------
using NonComparablePixels = std::vector<std::set<int>>;

enum class BitQuadType {
    StrictAncestor,
    Ancestor,
    StrictDescendant,
    Descendant
};

class BitQuadComparatorCT {
public:
    int rowOffset;
    int colOffset;
    std::function<bool(int, int, ComponentTree*, NonComparablePixels&)> comparator;
    BitQuadType type;

    bool isValid(int row, int col, ComponentTree* tree) const {
        return row + rowOffset >= 0 && row + rowOffset < tree->getNumRowsOfImage() && col + colOffset >= 0 && col + colOffset < tree->getNumColsOfImage();
    }

    BitQuadComparatorCT(int rowOffset, int colOffset, BitQuadType type) : rowOffset(rowOffset), colOffset(colOffset), type(type) {
        switch (type) {
            case BitQuadType::StrictAncestor:
                comparator = [=](int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA) {
                    
                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    
                    NodeId nodeP = tree->getSCById(idP);
                    NodeId nodeQ = tree->getSCById(idQ);
                    /*if(tree->isComparable(nodeP, nodeQ) == false) {
                        NodeId lca = tree->findLowestCommonAncestor(nodeP, nodeQ);
                        pixelsOfLCA[lca->getIndex()].push_back( ImageUtils::to1D(row, col, tree->getNumColsOfImage()) );
                        return false;
                    }*/

                    return tree->isStrictAncestor(nodeP, nodeQ);
                };
                break;
            case BitQuadType::Ancestor:
                comparator = [=](int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA) {
                    
                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    NodeId nodeP = tree->getSCById(idP);
                    NodeId nodeQ = tree->getSCById(idQ);
                    /*if(tree->isComparable(nodeP, nodeQ) == false) {
                        NodeId lca = tree->findLowestCommonAncestor(nodeP, nodeQ);
                        pixelsOfLCA[lca->getIndex()].push_back( ImageUtils::to1D(row, col, tree->getNumColsOfImage()) );
                        return false;
                    }*/
                    return tree->isAncestor(nodeP, nodeQ);
                };
                break;
            case BitQuadType::StrictDescendant:
                comparator = [=](int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA) {
                    
                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    NodeId nodeP = tree->getSCById(idP);
                    NodeId nodeQ = tree->getSCById(idQ);
                    /*if(tree->isComparable(nodeP, nodeQ) == false) {
                        NodeId lca = tree->findLowestCommonAncestor(nodeP, nodeQ);
                        pixelsOfLCA[lca->getIndex()].push_back( ImageUtils::to1D(row, col, tree->getNumColsOfImage()) );
                        return false;
                    }*/
                    return tree->isStrictDescendant(nodeP, nodeQ);
                };
                break;
            case BitQuadType::Descendant:
                comparator = [=](int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA) {
                    
                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    NodeId nodeP = tree->getSCById(idP);
                    NodeId nodeQ = tree->getSCById(idQ);
                    /*if(tree->isComparable(nodeP, nodeQ) == false) {
                        NodeId lca = tree->findLowestCommonAncestor(nodeP, nodeQ);
                        pixelsOfLCA[lca->getIndex()].push_back( ImageUtils::to1D(row, col, tree->getNumColsOfImage()) );
                        return false;
                    }*/
                    return tree->isDescendant(nodeP, nodeQ);
                };
                break;
        }
    }

    bool compare(int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA) const {
        if (!isValid(row, col, tree)){
            if(type == BitQuadType::StrictDescendant || type == BitQuadType::Descendant) {
                return true;
            }else{
                return false;
            }
        }
            
        return comparator(row, col, tree, pixelsOfLCA);
    }
};

//---------------------------------------------
// Padrão e grupo de padrões
//---------------------------------------------
class BitQuadCT {
    std::vector<BitQuadComparatorCT> quads;
public:
    BitQuadCT() = default;
    explicit BitQuadCT(size_t size) { quads.reserve(size); }
    
    BitQuadCT& add(BitQuadComparatorCT quad) {
        quads.push_back(quad);
        return *this;
    }
    // Função que retorna o símbolo Unicode para o tipo
    std::string symbolForType(BitQuadType type) {
        switch (type) {
            case BitQuadType::StrictAncestor:    return "A";
            case BitQuadType::Ancestor:          return "Ā";
            case BitQuadType::StrictDescendant:  return "D";
            case BitQuadType::Descendant:        return "Ḏ";
            default:                             return "?";
        }
    }

    // Imprime os padrões BitQuadCT em uma grade 3x3 com símbolos Unicode
    void print() {
        const int SIZE = 3;
        std::vector<std::vector<std::string>> matrix(SIZE, std::vector<std::string>(SIZE, " "));

        int center = SIZE / 2;
        matrix[center][center] = "o";

        for (const auto& quad : quads) {
            int row = center + quad.rowOffset;
            int col = center + quad.colOffset;

            if (row >= 0 && row < SIZE && col >= 0 && col < SIZE) {
                matrix[row][col] = symbolForType(quad.type);
            } else {
                std::cerr << "Aviso: posição fora do grid: (" << row << ", " << col << ")\n";
            }
        }

        // Impressão com grid (7x7 visual)
        for (int i = 0; i < SIZE; ++i) {
            std::cout << "+---+---+---+\n";
            std::cout << "| ";
            for (int j = 0; j < SIZE; ++j) {
                std::cout << matrix[i][j] << " | ";
            }
            std::cout << "\n";
        }
        std::cout << "+---+---+---+\n\n";
    }

    bool match(int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA) const {
        for (const auto& quad : quads) {
            if (!quad.compare(row, col, tree, pixelsOfLCA))
                return false;
        }
        return true;
    }
};

class BitQuadPatternCT {
    std::vector<BitQuadCT> patterns;
public:
    BitQuadPatternCT() = default;
    explicit BitQuadPatternCT(size_t size) { patterns.reserve(size); }
    
    BitQuadPatternCT& addBitQuad(const BitQuadCT& pattern) {
        patterns.push_back(pattern);
        return *this;
    }

    void print() {
        for (auto& bitquad: patterns) {
            bitquad.print();
        }
    }

    int count(int row, int col, ComponentTree* tree, NonComparablePixels& pixelsOfLCA)  {
        int c = 0;
        for (auto& pattern : patterns)
            if (pattern.match(row, col, tree, pixelsOfLCA)){
                ++c;
                //pattern.print(); 
            }
        return c;
    }
};

//---------------------------------------------
// AttributeBasedBitQuadsCT
//---------------------------------------------

struct AttributeBasedBitQuadsCT {
    int countPatternC1C4 = 0;
    int countPatternC1 = 0;
    int countPatternC2 = 0;
    int countPatternCD = 0;
    int countPatternC3 = 0;
    int countPatternC4 = 0;
    
    int countPatternCT1C4 = 0;
    int countPatternCT1 = 0;
    int countPatternCT2 = 0;
    int countPatternCTD = 0;
    int countPatternCT3 = 0;
    AdjacencyRelationPtr adj;

    std::string printPattern() const {
        std::ostringstream oss;
        oss << std::left
            << "Q1:" << std::setw(3) << countPatternC1
            << "  Q1C4:" << std::setw(3) << countPatternC1C4
            << "  Q2:" << std::setw(3) << countPatternC2
            << "  Q3:" << std::setw(3) << countPatternC3
            << "  QD:" << std::setw(3) << countPatternCD
            << "  Q4:" << std::setw(3) << countPatternC4
            << "  QT1C4:" << std::setw(3) << countPatternCT1C4
            << "  QT1:" << std::setw(3) << countPatternCT1
            << "  QT2:" << std::setw(3) << countPatternCT2
            << "  QT3:" << std::setw(3) << countPatternCT3
            << "  QTD:" << std::setw(3) << countPatternCTD;
        return oss.str();
    }

    AttributeBasedBitQuadsCT(AdjacencyRelationPtr adj) : adj(adj) {}

    int getNumberEuler() const {
        if (adj || adj->is4connectivity()) // ou use uma constante do seu projeto
            return (countPatternC1C4 - countPatternC3) / 4;
        else
            return (countPatternC1 - countPatternC3 - (2 * countPatternCD)) / 4;
    }

    int getNumberHoles() const {
        return 1 - getNumberEuler();
    }

    int getPerimeter() const {
        return countPatternC1 + countPatternC2 + countPatternC3 + (2 * countPatternCD);
    }

    int getArea() const {
        return (countPatternC1 + 2 * countPatternC2 + 3 * countPatternC3 + 4 * countPatternC4 + 2 * countPatternCD) / 4;
    }

    double getAreaDuda() const {
        return (1.0/4.0*countPatternC1 + 1.0/2.0*countPatternC2 + 7.0/8.0*countPatternC3 + countPatternC4 + 3.0/4.0*countPatternCD);
    }

    double getPerimeterContinuous() const {
        return countPatternC2 + ((countPatternC1 + countPatternC3) / 1.5);
    }

    double getCircularity() const {
        double area = getAreaDuda();
        double per = getPerimeterContinuous();
        return (4.0 * M_PI * area) / (per * per);
    }

    double getAreaAverage() const {
        double area = getAreaDuda();
        return area / static_cast<double>(getNumberEuler());
    }

    double getPerimeterAverage() const {
        return getPerimeterContinuous() / static_cast<double>(getNumberEuler());
    }

    double getLengthAverage() const {
        return getPerimeterAverage() / 2.0;
    }

    double getWidthAverage() const {
        return (2.0 * getAreaAverage()) / getPerimeterAverage();
    }

};



class ComputerAttributeBasedBitQuadsCT {
private:
    BitQuadPatternCT Q1;
    BitQuadPatternCT Q1C4;
    BitQuadPatternCT Q2;
    BitQuadPatternCT QD;
    BitQuadPatternCT Q3;
    BitQuadPatternCT Q4;

    BitQuadPatternCT Q1T;
    BitQuadPatternCT Q1C4T;
    BitQuadPatternCT Q2T;
    BitQuadPatternCT QDT;
    BitQuadPatternCT Q3T;

    ComponentTree* tree;
    AdjacencyRelationPtr adj;
    std::vector<AttributeBasedBitQuadsCT> attr;
    NonComparablePixels pixelsOfLCA;


    void initializePatterns() {
        if(!adj || adj->is4connectivity()) {
            createQ1C4Patterns();
            createQ1C4TPatterns();
        }else { // 8-connectivity
            createQ1Patterns();
            createQDPatterns();
            createQ1TPatterns();
            createQDTPattern();
        }
        createQ2Patterns();
        createQ3Patterns();
        createQ4Patterns();
        createQ2TPattern();
        createQ3TPattern();
    }

    void computerLocalPattern(NodeId nodeId, int p, std::vector<AttributeBasedBitQuadsCT>& attr) {
        auto [row,col] = ImageUtils::to2D(p, tree->getNumColsOfImage());
        

        if (!adj || adj->is4connectivity()) {
            attr[nodeId].countPatternC1C4  += Q1C4.count(row, col, tree, pixelsOfLCA);
            attr[nodeId].countPatternCT1C4 += Q1C4T.count(row, col, tree, pixelsOfLCA);
        } else { // 8-connectivity
            attr[nodeId].countPatternC1   += Q1.count(row, col, tree, pixelsOfLCA);
            attr[nodeId].countPatternCD   += QD.count(row, col, tree, pixelsOfLCA);
            attr[nodeId].countPatternCTD  += QDT.count(row, col, tree, pixelsOfLCA);
            attr[nodeId].countPatternCT1  += Q1T.count(row, col, tree, pixelsOfLCA);
        }

        attr[nodeId].countPatternC2  += Q2.count(row, col, tree, pixelsOfLCA);
        attr[nodeId].countPatternC3  += Q3.count(row, col, tree, pixelsOfLCA);
        attr[nodeId].countPatternC4  += Q4.count(row, col, tree, pixelsOfLCA);
        attr[nodeId].countPatternCT2 += Q2T.count(row, col, tree, pixelsOfLCA);
        attr[nodeId].countPatternCT3 += Q3T.count(row, col, tree, pixelsOfLCA);
    }

    void createQ1Patterns() {
        BitQuadCT Q1P1(3), Q1P2(3), Q1P3(3), Q1P4(3);

        Q1P1.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));

        Q1P2.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        Q1P3.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        Q1P4.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        Q1 = BitQuadPatternCT(4);
        Q1.addBitQuad(Q1P1).addBitQuad(Q1P2).addBitQuad(Q1P3).addBitQuad(Q1P4);
        if(PRINT_LOG){
            std::cout << "Q1 Patterns created:\n";
            Q1.print();
        }
    }

    void createQ1C4Patterns() {
        BitQuadCT Q1C4P1(2), Q1C4P2(2), Q1C4P3(2), Q1C4P4(2);

        Q1C4P1.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));
        

        Q1C4P2.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        Q1C4P3.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        Q1C4P4.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        Q1C4 = BitQuadPatternCT(4);
        Q1C4.addBitQuad(Q1C4P1).addBitQuad(Q1C4P2).addBitQuad(Q1C4P3).addBitQuad(Q1C4P4);
        if(PRINT_LOG){
            std::cout << "Q1C4 Patterns created:\n";
            Q1C4.print();
        }
    }


    void createQ2Patterns() {
        BitQuadCT Q2P1(3), Q2P2(3), Q2P3(3), Q2P4(3);
        BitQuadCT Q2P5(3), Q2P6(3), Q2P7(3), Q2P8(3);

        Q2P1.add(BitQuadComparatorCT(1, 0, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        Q2P2.add(BitQuadComparatorCT(0, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        Q2P3.add(BitQuadComparatorCT(-1, 0, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        Q2P4.add(BitQuadComparatorCT(0, -1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));

        Q2P5.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q2P6.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q2P7.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q2P8.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        Q2 = BitQuadPatternCT(8);
        Q2.addBitQuad(Q2P1).addBitQuad(Q2P2).addBitQuad(Q2P3).addBitQuad(Q2P4).addBitQuad(Q2P5).addBitQuad(Q2P6).addBitQuad(Q2P7).addBitQuad(Q2P8);
        if(PRINT_LOG){
            std::cout << "Q2 Patterns created:\n";
            Q2.print();
        }
    }

    void createQDPatterns() {
        BitQuadCT QDP1(3), QDP2(3), QDP3(3), QDP4(3);

        QDP1.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        QDP2.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));

        QDP3.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        QDP4.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        QD = BitQuadPatternCT(4);
        QD.addBitQuad(QDP1).addBitQuad(QDP2).addBitQuad(QDP3).addBitQuad(QDP4);
        if(PRINT_LOG){
            std::cout << "QD Patterns created:\n";
            QD.print();
        }
    }

    void createQ3Patterns() {
        BitQuadCT Q3P1(3), Q3P2(3), Q3P3(3), Q3P4(3);
        BitQuadCT Q3P5(3), Q3P6(3), Q3P7(3), Q3P8(3);
        BitQuadCT Q3P9(3), Q3P10(3), Q3P11(3), Q3P12(3);

        Q3P1.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q3P2.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        Q3P3.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q3P4.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q3P5.add(BitQuadComparatorCT(1, 0, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        Q3P6.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::Ancestor));

        Q3P7.add(BitQuadComparatorCT(-1, 0, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        Q3P8.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::Ancestor));

        Q3P9.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::Ancestor));

        Q3P10.add(BitQuadComparatorCT(0, -1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));

        Q3P11.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::Ancestor));

        Q3P12.add(BitQuadComparatorCT(0, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        Q3 = BitQuadPatternCT(12);
        Q3.addBitQuad(Q3P1).addBitQuad(Q3P2).addBitQuad(Q3P3).addBitQuad(Q3P4).addBitQuad(Q3P5).addBitQuad(Q3P6).addBitQuad(Q3P7).addBitQuad(Q3P8).addBitQuad(Q3P9).addBitQuad(Q3P10).addBitQuad(Q3P11).addBitQuad(Q3P12);
        if(PRINT_LOG){
            std::cout << "Q3 Patterns created:\n";
            Q3.print();
        }
    }

    void createQ4Patterns() {
        BitQuadCT Q4P1(3), Q4P2(3), Q4P3(3), Q4P4(3);

        Q4P1.add(BitQuadComparatorCT(1, 0, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::Ancestor));

        Q4P2.add(BitQuadComparatorCT(0, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q4P3.add(BitQuadComparatorCT(-1, 0, BitQuadType::Ancestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q4P4.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q4 = BitQuadPatternCT(4);
        Q4.addBitQuad(Q4P1).addBitQuad(Q4P2).addBitQuad(Q4P3).addBitQuad(Q4P4);
        
        if(PRINT_LOG){
            std::cout << "Q4 Patterns created:\n";
            Q4.print();
        }
    }

    void createQ1C4TPatterns() {
        BitQuadCT Q1C4TP1(2), Q1C4TP2(2), Q1C4TP3(2), Q1C4TP4(2);
        BitQuadCT Q1C4TP5(2), Q1C4TP6(2), Q1C4TP7(2), Q1C4TP8(2);

        Q1C4TP1.add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
                .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        Q1C4TP2.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor))
                .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictDescendant));

        Q1C4TP3.add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictDescendant))
                .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q1C4TP4.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
                .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant));

        Q1C4TP5.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor))
                .add(BitQuadComparatorCT(-1, -1, BitQuadType::Descendant));

        Q1C4TP6.add(BitQuadComparatorCT(1, -1, BitQuadType::Descendant))
                .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q1C4TP7.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor))
                .add(BitQuadComparatorCT(1, 1, BitQuadType::Descendant));

        Q1C4TP8.add(BitQuadComparatorCT(-1, 1, BitQuadType::Descendant))
                .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q1C4T = BitQuadPatternCT(8);
        Q1C4T.addBitQuad(Q1C4TP1).addBitQuad(Q1C4TP2).addBitQuad(Q1C4TP3).addBitQuad(Q1C4TP4).addBitQuad(Q1C4TP5).addBitQuad(Q1C4TP6).addBitQuad(Q1C4TP7).addBitQuad(Q1C4TP8);
        if(PRINT_LOG){
            std::cout << "Q1C4T Patterns created:\n";
            Q1C4T.print();
        }
    }

    void createQ1TPatterns() {
        BitQuadCT Q1TP1(3), Q1TP2(3), Q1TP3(3), Q1TP4(3);
        BitQuadCT Q1TP5(3), Q1TP6(3), Q1TP7(3), Q1TP8(3);
        BitQuadCT Q1TP9(3), Q1TP10(3), Q1TP11(3), Q1TP12(3);

        Q1TP1.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        Q1TP2.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        Q1TP3.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));

        Q1TP4.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        Q1TP5.add(BitQuadComparatorCT(1, 0, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        Q1TP6.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::Descendant));

        Q1TP7.add(BitQuadComparatorCT(-1, 0, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q1TP8.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::Descendant));

        Q1TP9.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::Descendant));

        Q1TP10.add(BitQuadComparatorCT(0, -1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q1TP11.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::Descendant));

        Q1TP12.add(BitQuadComparatorCT(0, 1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q1T = BitQuadPatternCT(12);
        Q1T.addBitQuad(Q1TP1).addBitQuad(Q1TP2).addBitQuad(Q1TP3).addBitQuad(Q1TP4).addBitQuad(Q1TP5).addBitQuad(Q1TP6).addBitQuad(Q1TP7).addBitQuad(Q1TP8).addBitQuad(Q1TP9).addBitQuad(Q1TP10).addBitQuad(Q1TP11).addBitQuad(Q1TP12);
        if(PRINT_LOG){
            std::cout << "Q1T Patterns created:\n";
            Q1T.print();
        }
    }

    void createQ2TPattern() {
        BitQuadCT Q2TP1(3), Q2TP2(3), Q2TP3(3), Q2TP4(3);
        BitQuadCT Q2TP5(3), Q2TP6(3), Q2TP7(3), Q2TP8(3);

        Q2TP1.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictDescendant));

        Q2TP2.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictDescendant));

        Q2TP3.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictDescendant));

        Q2TP4.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictDescendant));

        Q2TP5.add(BitQuadComparatorCT(-1, 0, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q2TP6.add(BitQuadComparatorCT(0, -1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q2TP7.add(BitQuadComparatorCT(1, 0, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        Q2TP8.add(BitQuadComparatorCT(0, 1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q2T = BitQuadPatternCT(8);
        Q2T.addBitQuad(Q2TP1).addBitQuad(Q2TP2).addBitQuad(Q2TP3).addBitQuad(Q2TP4).addBitQuad(Q2TP5).addBitQuad(Q2TP6).addBitQuad(Q2TP7).addBitQuad(Q2TP8);
        if(PRINT_LOG){
            std::cout << "Q2T Patterns created:\n";
            Q2T.print();
        }
    }

    void createQDTPattern() {
        BitQuadCT QDTP1(3), QDTP2(3), QDTP3(3), QDTP4(3);

        QDTP1.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        QDTP2.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictDescendant))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        QDTP3.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        QDTP4.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::Descendant))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        QDT = BitQuadPatternCT(4);
        QDT.addBitQuad(QDTP1).addBitQuad(QDTP2).addBitQuad(QDTP3).addBitQuad(QDTP4);
        if(PRINT_LOG){
            std::cout << "QDT Patterns created:\n";
            QDT.print();
        }
    }

    void createQ3TPattern() {
        BitQuadCT Q3TP1(3), Q3TP2(3), Q3TP3(3), Q3TP4(3);

        Q3TP1.add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor));

        Q3TP2.add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, -1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, -1, BitQuadType::StrictAncestor));

        Q3TP3.add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(-1, 0, BitQuadType::StrictAncestor));

        Q3TP4.add(BitQuadComparatorCT(1, 0, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(1, 1, BitQuadType::StrictAncestor))
            .add(BitQuadComparatorCT(0, 1, BitQuadType::StrictAncestor));

        Q3T = BitQuadPatternCT(4);
        Q3T.addBitQuad(Q3TP1).addBitQuad(Q3TP2).addBitQuad(Q3TP3).addBitQuad(Q3TP4);
        if(PRINT_LOG){
            std::cout << "Q3T Patterns created:\n";
            Q3T.print();
        }
    }

public:


    // Construtor principal
    ComputerAttributeBasedBitQuadsCT(ComponentTreePtr tree): ComputerAttributeBasedBitQuadsCT(tree.get()) {}
    ComputerAttributeBasedBitQuadsCT(ComponentTree* tree) : tree(tree), adj(tree->getAdjacencyRelation()), attr(tree->getNumNodes(), AttributeBasedBitQuadsCT(adj)) {
        
        assert(tree->getTreeType() != ComponentTree::TREE_OF_SHAPES && "Não está implementado para tree of shapes!");
        
        initializePatterns();
        AttributeComputedIncrementallyCT::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                for(int p: tree->getCNPsById(node)){
                    computerLocalPattern(node, p, attr) ;
                }	
            },
            [&](NodeId parent, NodeId child) {
                if (!adj || adj->is4connectivity()) // 4-connectivity
                    attr[parent].countPatternC1C4 += attr[child].countPatternC1C4;
                else {
                    attr[parent].countPatternC1 += attr[child].countPatternC1;
                    attr[parent].countPatternCD += attr[child].countPatternCD;
                }
                attr[parent].countPatternC2 += attr[child].countPatternC2;
                attr[parent].countPatternC3 += attr[child].countPatternC3;
                attr[parent].countPatternC4 += attr[child].countPatternC4;
            },
            [&](NodeId node) {

                /*std::vector<int>& pixelsNonComparable = pixelsOfLCA[node];
                std::cout << "Node: " << node << " - Non-comparable pixels: " << pixelsNonComparable.size() << std::endl;
                for (int p : pixelsNonComparable) {
                    computerLocalPattern(node, p, attr);
                }*/

                if (!adj || adj->is4connectivity()) 
                    attr[node].countPatternC1C4 = attr[node].countPatternC1C4 - attr[node].countPatternCT1C4;
                else {
                    attr[node].countPatternC1 = attr[node].countPatternC1 - attr[node].countPatternCT1;
                    attr[node].countPatternCD = attr[node].countPatternCD - attr[node].countPatternCTD;
                }
                attr[node].countPatternC2 = attr[node].countPatternC2 - attr[node].countPatternCT2;
                attr[node].countPatternC3 = attr[node].countPatternC3 - attr[node].countPatternCT3;
            }
        );

    }

    std::vector<AttributeBasedBitQuadsCT> getAttributes() const {
        return attr;
    }

};


