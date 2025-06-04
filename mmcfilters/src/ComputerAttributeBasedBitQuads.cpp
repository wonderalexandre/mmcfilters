#include "../include/ComputerAttributeBasedBitQuads.hpp"



    void ComputerAttributeBasedBitQuads::initializePatterns() {
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

    void ComputerAttributeBasedBitQuads::computerLocalPattern(NodeMTPtr node, int p, std::vector<AttributeBasedBitQuads>& attr) {
        auto [row,col] = ImageUtils::to2D(p, tree->getNumColsOfImage());
        int nodeId = node->getIndex();

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

    void ComputerAttributeBasedBitQuads::createQ1Patterns() {
        BitQuad Q1P1(3), Q1P2(3), Q1P3(3), Q1P4(3);

        Q1P1.add(BitQuadComparator(0, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));

        Q1P2.add(BitQuadComparator(-1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        Q1P3.add(BitQuadComparator(0, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        Q1P4.add(BitQuadComparator(1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        Q1 = BitQuadPattern(4);
        Q1.addBitQuad(Q1P1).addBitQuad(Q1P2).addBitQuad(Q1P3).addBitQuad(Q1P4);
        if(PRINT_LOG){
            std::cout << "Q1 Patterns created:\n";
            Q1.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ1C4Patterns() {
        BitQuad Q1C4P1(2), Q1C4P2(2), Q1C4P3(2), Q1C4P4(2);

        Q1C4P1.add(BitQuadComparator(0, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));
        

        Q1C4P2.add(BitQuadComparator(-1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        Q1C4P3.add(BitQuadComparator(0, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        Q1C4P4.add(BitQuadComparator(1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        Q1C4 = BitQuadPattern(4);
        Q1C4.addBitQuad(Q1C4P1).addBitQuad(Q1C4P2).addBitQuad(Q1C4P3).addBitQuad(Q1C4P4);
        if(PRINT_LOG){
            std::cout << "Q1C4 Patterns created:\n";
            Q1C4.print();
        }
    }


    void ComputerAttributeBasedBitQuads::createQ2Patterns() {
        BitQuad Q2P1(3), Q2P2(3), Q2P3(3), Q2P4(3);
        BitQuad Q2P5(3), Q2P6(3), Q2P7(3), Q2P8(3);

        Q2P1.add(BitQuadComparator(1, 0, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        Q2P2.add(BitQuadComparator(0, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        Q2P3.add(BitQuadComparator(-1, 0, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        Q2P4.add(BitQuadComparator(0, -1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));

        Q2P5.add(BitQuadComparator(0, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q2P6.add(BitQuadComparator(-1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q2P7.add(BitQuadComparator(0, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q2P8.add(BitQuadComparator(1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        Q2 = BitQuadPattern(8);
        Q2.addBitQuad(Q2P1).addBitQuad(Q2P2).addBitQuad(Q2P3).addBitQuad(Q2P4).addBitQuad(Q2P5).addBitQuad(Q2P6).addBitQuad(Q2P7).addBitQuad(Q2P8);
        if(PRINT_LOG){
            std::cout << "Q2 Patterns created:\n";
            Q2.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQDPatterns() {
        BitQuad QDP1(3), QDP2(3), QDP3(3), QDP4(3);

        QDP1.add(BitQuadComparator(1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(1, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        QDP2.add(BitQuadComparator(0, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, -1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));

        QDP3.add(BitQuadComparator(-1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        QDP4.add(BitQuadComparator(0, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        QD = BitQuadPattern(4);
        QD.addBitQuad(QDP1).addBitQuad(QDP2).addBitQuad(QDP3).addBitQuad(QDP4);
        if(PRINT_LOG){
            std::cout << "QD Patterns created:\n";
            QD.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ3Patterns() {
        BitQuad Q3P1(3), Q3P2(3), Q3P3(3), Q3P4(3);
        BitQuad Q3P5(3), Q3P6(3), Q3P7(3), Q3P8(3);
        BitQuad Q3P9(3), Q3P10(3), Q3P11(3), Q3P12(3);

        Q3P1.add(BitQuadComparator(0, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q3P2.add(BitQuadComparator(1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        Q3P3.add(BitQuadComparator(0, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q3P4.add(BitQuadComparator(-1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q3P5.add(BitQuadComparator(1, 0, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        Q3P6.add(BitQuadComparator(0, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 0, BitQuadType::SubsetEq));

        Q3P7.add(BitQuadComparator(-1, 0, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        Q3P8.add(BitQuadComparator(0, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, 0, BitQuadType::SubsetEq));

        Q3P9.add(BitQuadComparator(-1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(-1, -1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(0, -1, BitQuadType::SubsetEq));

        Q3P10.add(BitQuadComparator(0, -1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, -1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));

        Q3P11.add(BitQuadComparator(1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(1, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(0, 1, BitQuadType::SubsetEq));

        Q3P12.add(BitQuadComparator(0, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        Q3 = BitQuadPattern(12);
        Q3.addBitQuad(Q3P1).addBitQuad(Q3P2).addBitQuad(Q3P3).addBitQuad(Q3P4).addBitQuad(Q3P5).addBitQuad(Q3P6).addBitQuad(Q3P7).addBitQuad(Q3P8).addBitQuad(Q3P9).addBitQuad(Q3P10).addBitQuad(Q3P11).addBitQuad(Q3P12);
        if(PRINT_LOG){
            std::cout << "Q3 Patterns created:\n";
            Q3.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ4Patterns() {
        BitQuad Q4P1(3), Q4P2(3), Q4P3(3), Q4P4(3);

        Q4P1.add(BitQuadComparator(1, 0, BitQuadType::SubsetEq))
            .add(BitQuadComparator(1, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(0, 1, BitQuadType::SubsetEq));

        Q4P2.add(BitQuadComparator(0, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, 1, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q4P3.add(BitQuadComparator(-1, 0, BitQuadType::SubsetEq))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q4P4.add(BitQuadComparator(0, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q4 = BitQuadPattern(4);
        Q4.addBitQuad(Q4P1).addBitQuad(Q4P2).addBitQuad(Q4P3).addBitQuad(Q4P4);
        
        if(PRINT_LOG){
            std::cout << "Q4 Patterns created:\n";
            Q4.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ1C4TPatterns() {
        BitQuad Q1C4TP1(2), Q1C4TP2(2), Q1C4TP3(2), Q1C4TP4(2);
        BitQuad Q1C4TP5(2), Q1C4TP6(2), Q1C4TP7(2), Q1C4TP8(2);

        Q1C4TP1.add(BitQuadComparator(1, 1, BitQuadType::Superset))
               .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        Q1C4TP2.add(BitQuadComparator(0, 1, BitQuadType::Subset))
               .add(BitQuadComparator(-1, 1, BitQuadType::Subset));

        Q1C4TP3.add(BitQuadComparator(-1, -1, BitQuadType::Superset))
               .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q1C4TP4.add(BitQuadComparator(0, -1, BitQuadType::Subset))
               .add(BitQuadComparator(1, -1, BitQuadType::Superset));

        Q1C4TP5.add(BitQuadComparator(-1, 0, BitQuadType::Subset))
               .add(BitQuadComparator(-1, -1, BitQuadType::SupersetEq));

        Q1C4TP6.add(BitQuadComparator(1, -1, BitQuadType::SupersetEq))
               .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q1C4TP7.add(BitQuadComparator(1, 0, BitQuadType::Subset))
               .add(BitQuadComparator(1, 1, BitQuadType::SupersetEq));

        Q1C4TP8.add(BitQuadComparator(-1, 1, BitQuadType::SupersetEq))
               .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q1C4T = BitQuadPattern(8);
        Q1C4T.addBitQuad(Q1C4TP1).addBitQuad(Q1C4TP2).addBitQuad(Q1C4TP3).addBitQuad(Q1C4TP4).addBitQuad(Q1C4TP5).addBitQuad(Q1C4TP6).addBitQuad(Q1C4TP7).addBitQuad(Q1C4TP8);
        if(PRINT_LOG){
            std::cout << "Q1C4T Patterns created:\n";
            Q1C4T.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ1TPatterns() {
        BitQuad Q1TP1(3), Q1TP2(3), Q1TP3(3), Q1TP4(3);
        BitQuad Q1TP5(3), Q1TP6(3), Q1TP7(3), Q1TP8(3);
        BitQuad Q1TP9(3), Q1TP10(3), Q1TP11(3), Q1TP12(3);

        Q1TP1.add(BitQuadComparator(0, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        Q1TP2.add(BitQuadComparator(1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        Q1TP3.add(BitQuadComparator(0, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));

        Q1TP4.add(BitQuadComparator(-1, 0, BitQuadType::Superset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        Q1TP5.add(BitQuadComparator(1, 0, BitQuadType::SupersetEq))
            .add(BitQuadComparator(1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        Q1TP6.add(BitQuadComparator(0, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(-1, 0, BitQuadType::SupersetEq));

        Q1TP7.add(BitQuadComparator(-1, 0, BitQuadType::SupersetEq))
            .add(BitQuadComparator(-1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q1TP8.add(BitQuadComparator(0, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::SupersetEq));

        Q1TP9.add(BitQuadComparator(-1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(-1, -1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(0, -1, BitQuadType::SupersetEq));

        Q1TP10.add(BitQuadComparator(0, -1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(1, -1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q1TP11.add(BitQuadComparator(1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(1, 1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(0, 1, BitQuadType::SupersetEq));

        Q1TP12.add(BitQuadComparator(0, 1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(-1, 1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q1T = BitQuadPattern(12);
        Q1T.addBitQuad(Q1TP1).addBitQuad(Q1TP2).addBitQuad(Q1TP3).addBitQuad(Q1TP4).addBitQuad(Q1TP5).addBitQuad(Q1TP6).addBitQuad(Q1TP7).addBitQuad(Q1TP8).addBitQuad(Q1TP9).addBitQuad(Q1TP10).addBitQuad(Q1TP11).addBitQuad(Q1TP12);
        if(PRINT_LOG){
            std::cout << "Q1T Patterns created:\n";
            Q1T.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ2TPattern() {
        BitQuad Q2TP1(3), Q2TP2(3), Q2TP3(3), Q2TP4(3);
        BitQuad Q2TP5(3), Q2TP6(3), Q2TP7(3), Q2TP8(3);

        Q2TP1.add(BitQuadComparator(0, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, 0, BitQuadType::Superset));

        Q2TP2.add(BitQuadComparator(1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(0, 1, BitQuadType::Superset));

        Q2TP3.add(BitQuadComparator(0, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Superset));

        Q2TP4.add(BitQuadComparator(-1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Superset));

        Q2TP5.add(BitQuadComparator(-1, 0, BitQuadType::SupersetEq))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q2TP6.add(BitQuadComparator(0, -1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q2TP7.add(BitQuadComparator(1, 0, BitQuadType::SupersetEq))
            .add(BitQuadComparator(1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        Q2TP8.add(BitQuadComparator(0, 1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(-1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q2T = BitQuadPattern(8);
        Q2T.addBitQuad(Q2TP1).addBitQuad(Q2TP2).addBitQuad(Q2TP3).addBitQuad(Q2TP4).addBitQuad(Q2TP5).addBitQuad(Q2TP6).addBitQuad(Q2TP7).addBitQuad(Q2TP8);
        if(PRINT_LOG){
            std::cout << "Q2T Patterns created:\n";
            Q2T.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQDTPattern() {
        BitQuad QDTP1(3), QDTP2(3), QDTP3(3), QDTP4(3);

        QDTP1.add(BitQuadComparator(0, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, -1, BitQuadType::Superset))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        QDTP2.add(BitQuadComparator(1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(1, 1, BitQuadType::Superset))
            .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        QDTP3.add(BitQuadComparator(0, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        QDTP4.add(BitQuadComparator(-1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(-1, -1, BitQuadType::SupersetEq))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        QDT = BitQuadPattern(4);
        QDT.addBitQuad(QDTP1).addBitQuad(QDTP2).addBitQuad(QDTP3).addBitQuad(QDTP4);
        if(PRINT_LOG){
            std::cout << "QDT Patterns created:\n";
            QDT.print();
        }
    }

    void ComputerAttributeBasedBitQuads::createQ3TPattern() {
        BitQuad Q3TP1(3), Q3TP2(3), Q3TP3(3), Q3TP4(3);

        Q3TP1.add(BitQuadComparator(0, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(1, 0, BitQuadType::Subset));

        Q3TP2.add(BitQuadComparator(-1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(-1, -1, BitQuadType::Subset))
            .add(BitQuadComparator(0, -1, BitQuadType::Subset));

        Q3TP3.add(BitQuadComparator(0, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(-1, 0, BitQuadType::Subset));

        Q3TP4.add(BitQuadComparator(1, 0, BitQuadType::Subset))
            .add(BitQuadComparator(1, 1, BitQuadType::Subset))
            .add(BitQuadComparator(0, 1, BitQuadType::Subset));

        Q3T = BitQuadPattern(4);
        Q3T.addBitQuad(Q3TP1).addBitQuad(Q3TP2).addBitQuad(Q3TP3).addBitQuad(Q3TP4);
        if(PRINT_LOG){
            std::cout << "Q3T Patterns created:\n";
            Q3T.print();
        }
    }