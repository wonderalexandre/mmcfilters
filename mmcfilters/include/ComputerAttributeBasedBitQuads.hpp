
#include <vector>
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
#include "../include/MorphologicalTree.hpp"
#include "../include/NodeMT.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/ImageUtils.hpp"

#ifndef COMPUTER_ATTRIBUTE_BASED_BIT_QUADS_HPP
#define COMPUTER_ATTRIBUTE_BASED_BIT_QUADS_HPP


//---------------------------------------------
// CLASSES QuadBit e padrões
//---------------------------------------------
enum class BitQuadType {
    StrictAncestor,
    Ancestor,
    StrictDescendant,
    Descendant
};

class BitQuadComparator {
public:
    int rowOffset;
    int colOffset;
    std::function<bool(int, int, MorphologicalTreePtr, bool)> comparator;



    BitQuadComparator(int rowOffset, int colOffset, BitQuadType type) : rowOffset(rowOffset), colOffset(colOffset){


        switch (type) {
            case BitQuadType::StrictAncestor:
                comparator = [=](int row, int col, MorphologicalTreePtr tree, bool isMaxtree) {

                    bool isValid = (row + rowOffset >= 0 && row + rowOffset < tree->getNumRowsOfImage() && col + colOffset >= 0 && col + colOffset < tree->getNumColsOfImage());
                    if (!isValid) return false;

                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    return tree->isStrictAncestor(tree->getSC(idP), tree->getSC(idQ));
                    
                };
                break;
            case BitQuadType::Ancestor:
                comparator = [=](int row, int col, MorphologicalTreePtr tree, bool isMaxtree) {
                    bool isValid = (row + rowOffset >= 0 && row + rowOffset < tree->getNumRowsOfImage() && col + colOffset >= 0 && col + colOffset < tree->getNumColsOfImage());
                    if (!isValid) return false;

                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    return tree->isAncestor(tree->getSC(idP), tree->getSC(idQ));
                };
                break;
            case BitQuadType::StrictDescendant:
                comparator = [=](int row, int col, MorphologicalTreePtr tree, bool isMaxtree) {
                    bool isValid = (row + rowOffset >= 0 && row + rowOffset < tree->getNumRowsOfImage() && col + colOffset >= 0 && col + colOffset < tree->getNumColsOfImage());
                    if (!isValid) return true;

                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    return tree->isStrictDescendant(tree->getSC(idP), tree->getSC(idQ));
                };
                break;
            case BitQuadType::Descendant:
                comparator = [=](int row, int col, MorphologicalTreePtr tree, bool isMaxtree) {
                    bool isValid = (row + rowOffset >= 0 && row + rowOffset < tree->getNumRowsOfImage() && col + colOffset >= 0 && col + colOffset < tree->getNumColsOfImage());
                    if (!isValid) return true;
                    
                    auto idP = ImageUtils::to1D(row, col, tree->getNumColsOfImage());
                    auto idQ = ImageUtils::to1D(row + rowOffset, col + colOffset, tree->getNumColsOfImage());
                    return tree->isDescendant(tree->getSC(idP), tree->getSC(idQ));
                };
                break;
        }
    }

    bool compare(int row, int col, MorphologicalTreePtr tree, bool isMaxtree) const {
        return comparator(row, col, tree, isMaxtree);
    }
};

//---------------------------------------------
// Padrão e grupo de padrões
//---------------------------------------------

class BitQuad {
    std::vector<std::shared_ptr<BitQuadComparator>> quads;
public:
    BitQuad() = default;
    explicit BitQuad(size_t size) { quads.reserve(size); }
    
    BitQuad& add(std::shared_ptr<BitQuadComparator> quad) {
        quads.push_back(quad);
        return *this;
    }

    bool match(int row, int col, MorphologicalTreePtr tree, bool isMaxtree) const {
        for (const auto& q : quads) {
            if (!q->compare(row, col, tree, isMaxtree))
                return false;
        }
        return true;
    }
};

class BitQuadPattern {
    std::vector<BitQuad> patterns;
public:
    BitQuadPattern() = default;
    explicit BitQuadPattern(size_t size) { patterns.reserve(size); }
    
    BitQuadPattern& addBitQuad(const BitQuad& pattern) {
        patterns.push_back(pattern);
        return *this;
    }

    int count(int row, int col, MorphologicalTreePtr tree, bool isMaxtree) const {
        int c = 0;
        for (const auto& pattern : patterns)
            if (pattern.match(row, col, tree, isMaxtree))
                ++c;
        return c;
    }
};

//---------------------------------------------
// AttributeBasedBitQuads
//---------------------------------------------

struct AttributeBasedBitQuads {
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
};
using AttributeBasedBitQuadsPtr = std::shared_ptr<AttributeBasedBitQuads>;



class ComputerAttributeBasedBitQuads {
private:
    // Patterns para conectividade
    BitQuadPattern Q1;
    BitQuadPattern Q1C4;
    BitQuadPattern Q2;
    BitQuadPattern QD;
    BitQuadPattern Q3;
    BitQuadPattern Q4;

    BitQuadPattern Q1T;
    BitQuadPattern Q1C4T;
    BitQuadPattern Q2T;
    BitQuadPattern QDT;
    BitQuadPattern Q3T;

    // Imagem e adjacência
    MorphologicalTreePtr tree;
    AdjacencyRelationPtr adj;
    bool isMaxtree;

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

    void computerLocalPattern(NodeMTPtr node, int p, std::vector<AttributeBasedBitQuadsPtr>& attr);

    std::vector<AttributeBasedBitQuadsPtr> attr;

public:


    // Construtor principal
    ComputerAttributeBasedBitQuads(MorphologicalTreePtr tree) : tree(tree), adj(tree->getAdjacencyRelation()), attr(tree->getNumNodes()) {
        initializePatterns();
        AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&](NodeMTPtr node) {
			    attr[node->getIndex()] = std::make_shared<AttributeBasedBitQuads>();
		        for(int p: node->getCNPs()){
			        computerLocalPattern(node, p, attr) ;
		        }	

			},
			[&](NodeMTPtr parent, NodeMTPtr child) {
                if (adj->is4connectivity()) // 4-connectivity
                    attr[parent->getIndex()]->countPatternC1C4 += attr[child->getIndex()]->countPatternC1C4;
                else {
                    attr[parent->getIndex()]->countPatternC1 += attr[child->getIndex()]->countPatternC1;
                    attr[parent->getIndex()]->countPatternCD += attr[child->getIndex()]->countPatternCD;
                }
                attr[parent->getIndex()]->countPatternC2 += attr[child->getIndex()]->countPatternC2;
                attr[parent->getIndex()]->countPatternC3 += attr[child->getIndex()]->countPatternC3;
                attr[parent->getIndex()]->countPatternC4 += attr[child->getIndex()]->countPatternC4;
            },
			[&](NodeMTPtr node) {
                if (adj->is4connectivity()) 
                    attr[node->getIndex()]->countPatternC1C4 = attr[node->getIndex()]->countPatternC1C4 - attr[node->getIndex()]->countPatternCT1C4;
                else {
                    attr[node->getIndex()]->countPatternC1 = attr[node->getIndex()]->countPatternC1 - attr[node->getIndex()]->countPatternCT1;
                    attr[node->getIndex()]->countPatternCD = attr[node->getIndex()]->countPatternCD - attr[node->getIndex()]->countPatternCTD;
                }
                
                attr[node->getIndex()]->countPatternC2 = attr[node->getIndex()]->countPatternC2 - attr[node->getIndex()]->countPatternCT2;
                attr[node->getIndex()]->countPatternC3 = attr[node->getIndex()]->countPatternC3 - attr[node->getIndex()]->countPatternCT3;
                
            }
		);

    }

    std::vector<AttributeBasedBitQuadsPtr> getAttributes() const {
        return attr;
    }


};


#endif // COMPUTER_ATTRIBUTE_BASED_BIT_QUADS_HPP