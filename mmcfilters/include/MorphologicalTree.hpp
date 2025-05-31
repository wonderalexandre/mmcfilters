#include <list>
#include <vector>

#include "../include/NodeMT.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/BuilderComponentTreeByUnionFind.hpp"
#include "../include/Common.hpp"
#include <iostream>
#include <cmath>
#include <type_traits>

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H

class MorphologicalTree;
using MorphologicalTreePtr = std::shared_ptr<MorphologicalTree>; 

class MorphologicalTree {

protected:
	int numRows;
    int numCols;
	int treeType; //0-mintree, 1-maxtree, 2-tree of shapes
	NodeMTPtr root;
	int numNodes;
	std::vector<NodeMTPtr> indexToNode;
	std::vector<NodeMTPtr> nodes;
	AdjacencyRelationPtr adj;
	int depth;
	
	void reconstruction(NodeMTPtr node, uint8_t* data);
	void computerTreeAttribute();
	MorphologicalTree(){ }
    
    
public:
   	static const int MAX_TREE = 0;
	static const int MIN_TREE = 1;
	static const int TREE_OF_SHAPES = 2;

	MorphologicalTree(ImageUInt8Ptr img, std::string ToSInperpolation="self-dual");
	explicit MorphologicalTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5);
    MorphologicalTree(ImageUInt8Ptr img, const char* ToSInterpolation) : MorphologicalTree(img, std::string(ToSInterpolation)) {}

    ~MorphologicalTree();

    static MorphologicalTreePtr create(int rows, int cols, bool isMaxtree, AdjacencyRelationPtr adj) {
        struct Enabler : public MorphologicalTree {
            Enabler(int r, int c, bool m, AdjacencyRelationPtr a) {
                this->numRows = r;
                this->numCols = c;
                this->treeType = m ? MAX_TREE : MIN_TREE;
                this->adj = a;
            }
        };
        return std::make_shared<Enabler>(rows, cols, isMaxtree, adj);
    }
    

    template <typename PixelType>
    static MorphologicalTreePtr createFromAttributeMapping(ImagePtr<PixelType> attrMappingPtr, ImageUInt8Ptr imgPtr, bool isMaxtree, double radius) {
        AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(imgPtr->getNumRows(), imgPtr->getNumCols(), radius);	
        MorphologicalTreePtr tree = MorphologicalTree::create(imgPtr->getNumRows(), imgPtr->getNumCols(), isMaxtree, adj);
        BuilderComponentTreeByUnionFind<PixelType> builder(attrMappingPtr, isMaxtree, adj);
        int n = imgPtr->getNumRows() * imgPtr->getNumCols();
        auto img = imgPtr->rawData();
        auto attrMapping = attrMappingPtr->rawData();

        auto orderedPixels = builder.getOrderedPixels();
        auto parent = builder.getParent();
        float epsilon = 1e-5f; // Tolerância para comparação de pixels flutuantes

        tree->nodes.resize(n, nullptr);

        tree->numNodes = 0;
        for (int i = 0; i < n; i++) {
            int p = orderedPixels[i];
            if (p == parent[p]) { //representante do node raiz
                tree->root = tree->nodes[p] = std::make_shared<NodeMT>(tree->numNodes++, nullptr, img[p]);
                tree->nodes[p]->addCNPs(p);
            }
            if constexpr (std::is_floating_point<PixelType>::value) {            
                if (std::fabs(attrMapping[p] - attrMapping[parent[p]]) > epsilon) { //representante de um node
                    tree->nodes[p] = std::make_shared<NodeMT>(tree->numNodes++, tree->nodes[parent[p]], img[p]);
                    tree->nodes[p]->addCNPs(p);
                    tree->nodes[parent[p]]->addChild(tree->nodes[p]);
                }
                else {
                    tree->nodes[parent[p]]->addCNPs(p);
                    tree->nodes[p] = tree->nodes[parent[p]];
                }
            }
            else{
                if (attrMapping[p] != attrMapping[parent[p]]) { //representante de um node
                    tree->nodes[p] = std::make_shared<NodeMT>(tree->numNodes++, tree->nodes[parent[p]], img[p]);
                    tree->nodes[p]->addCNPs(p);
                    tree->nodes[parent[p]]->addChild(tree->nodes[p]);
                }
                else {
                    tree->nodes[parent[p]]->addCNPs(p);
                    tree->nodes[p] = tree->nodes[parent[p]];
                }
            }
        }
        
        tree->computerTreeAttribute();

        //ajustar o level de cada nó
        for(NodeMTPtr node: tree->getIndexNode()){
            int media = 0;
            for(int p: node->getCNPs()){
                media += img[p];
            }
            media /= node->getCNPs().size();
            node->setLevel(media);
        }


        return tree;

    }


	NodeMTPtr getRoot();

	bool isMaxtree();

	int getTreeType();

	NodeMTPtr getSC(int pixel);

    NodeMTPtr getNodeByIndex(int index);

	std::vector<NodeMTPtr>& getIndexNode();

	int getNumNodes();

	int getNumRowsOfImage();

	int getNumColsOfImage();

	ImageUInt8Ptr reconstructionImage();

	ImageUInt8Ptr getImageAferPruning(NodeMTPtr node);

	void pruning(NodeMTPtr node);

	bool isAncestor(NodeMTPtr u, NodeMTPtr v);
	
	bool isDescendant(NodeMTPtr u, NodeMTPtr v);

	bool isComparable(NodeMTPtr u, NodeMTPtr v);

	bool isStrictAncestor(NodeMTPtr u, NodeMTPtr v);
	
	bool isStrictDescendant(NodeMTPtr u, NodeMTPtr v);
	
	bool isStrictComparable(NodeMTPtr u, NodeMTPtr v);
	
	NodeMTPtr findLowestCommonAncestor(NodeMTPtr u, NodeMTPtr v);

	int getDepth();

    std::list<NodeMTPtr> getLeaves();

	std::vector<std::vector<NodeMTPtr>> getNodesByDepth();

	static void extractDepthMap(NodeMTPtr node, int depth, std::vector<std::vector<NodeMTPtr>>& nodesByDepth){
		nodesByDepth[depth].push_back(node);
		for (NodeMTPtr child : node->getChildren()) {
		  extractDepthMap(child, depth + 1, nodesByDepth);
		}
	}
    
};


/**
 * Método Euler Tour + RMQ
  
 Etapa 1: Euler Tour
    Realiza um DFS na árvore e registra:
	1.	A ordem dos nós visitados → euler[]
	2.	A profundidade de cada nó na árvore durante o percurso → depth[]
	3.	A índice da primeira ocorrência de cada nó no vetor euler → firstOccurrence[]

  Etapa 2: RMQ na profundidade
    Para responder LCA(u, v):
	1.	Pegue i = firstOccurrence[u], 
              j = firstOccurrence[v]
	2.	Realize um RMQ (Range Minimum Query) sobre o vetor depth[] entre as posições min(i, j) e max(i, j) no vetor euler[].
	3.	O resultado do RMQ será o índice do nó com menor profundidade entre u e v no caminho — ou seja, o LCA!

    Exemplo:
      0
     / \
    1   2
   /
  3
  Índices:         0  1  2  3  4  5  6
  euler =         [0, 1, 3, 1, 0, 2, 0]
  depth =         [0, 1, 2, 1, 0, 1, 0]
  firstOccurrence=[0, 1, 5, 2         ]
    
  LCA(3, 2) = 0
    i = firstOccurrence[3] = 2
    j = firstOccurrence[2] = 5
    RMQ: 
      1. Descobrir o intervalo no vetor depth: depth[2..5] = [2, 1, 0, 1]
      2. Encontrar a posição do menor valor: O mínimo é 0, que ocorre em depth[4]
      3. O correspondente em no vetor euler: euler[4] = 0 que é o indice do LCA
	
 */

class LCAEulerRMQ {
private:
    std::vector<int> euler;            // timePreOrder dos nós na ordem de visita
    std::vector<int> depth;            // profundidade associada a cada posição em euler
    std::vector<int> firstOccurrence;  // [timePreOrder] = posição no vetor euler
    std::vector<std::vector<int>> st;  // Sparse Table para RMQ
    std::vector<NodeMTPtr> indexToNode;  // acesso direto aos nós via timePreOrder

public:
    LCAEulerRMQ(MorphologicalTreePtr tree) {
        indexToNode = tree->getIndexNode(); // indexado por timePreOrder
        int n = indexToNode.size();
        firstOccurrence.resize(n, -1);

        depthFirstTraversal(tree->getRoot(), 0);
        buildSparseTable();
    }

	NodeMTPtr findLowestCommonAncestor(NodeMTPtr u, NodeMTPtr v) {
        int uTime = u->getIndex();
        int vTime = v->getIndex();
        int i = firstOccurrence[uTime];
        int j = firstOccurrence[vTime];
        if (i > j) std::swap(i, j);
        int idx = rmq(i, j);
        return indexToNode[euler[idx]];
    }


private:
    void depthFirstTraversal(NodeMTPtr node, int d) {
        int time = node->getIndex();
        if (firstOccurrence[time] == -1)
            firstOccurrence[time] = euler.size();

        euler.push_back(time);
        depth.push_back(d);

        for (NodeMTPtr child : node->getChildren()) {
            depthFirstTraversal(child, d + 1);
            euler.push_back(time);
            depth.push_back(d);
        }
    }

    void buildSparseTable() {
        int n = depth.size();
        int logn = std::log2(n) + 1;
        st.assign(n, std::vector<int>(logn));

        for (int i = 0; i < n; ++i)
            st[i][0] = i;

        for (int j = 1; (1 << j) <= n; ++j) {
            for (int i = 0; i + (1 << j) <= n; ++i) {
                int l = st[i][j - 1];
                int r = st[i + (1 << (j - 1))][j - 1];
                st[i][j] = (depth[l] < depth[r]) ? l : r;
            }
        }
    }

    int rmq(int l, int r) {
        int len = r - l + 1;
        int k = std::log2(len);
        int a = st[l][k];
        int b = st[r - (1 << k) + 1][k];
        return (depth[a] < depth[b]) ? a : b;
    }
};


#endif