
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

#ifndef COMPUTER_ATTRIBUTE_BASED_BIT_QUADS_CT_HPP
#define COMPUTER_ATTRIBUTE_BASED_BIT_QUADS_CT_HPP

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

    void initializePatterns();

    void createQ1Patterns();
    void createQ1C4Patterns();
    void createQ2Patterns();
    void createQDPatterns();
    void createQ3Patterns();
    void createQ4Patterns();

    void createQ1C4TPatterns();
    void createQ1TPatterns();
    void createQ2TPattern();
    void createQDTPattern();
    void createQ3TPattern();

    void computerLocalPattern(NodeId node, int p, std::vector<AttributeBasedBitQuadsCT>& attr);

public:


    // Construtor principal
    ComputerAttributeBasedBitQuadsCT(ComponentTree* tree);

    std::vector<AttributeBasedBitQuadsCT> getAttributes() const;

    


};


#endif // COMPUTER_ATTRIBUTE_BASED_BIT_QUADS_CT_HPP