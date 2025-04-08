#include <list>
#include <vector>

#include "../include/NodeMT.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"
#include <iostream>

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H


class MorphologicalTree {

protected:
	int numCols;
	int numRows;
	int treeType; //0-mintree, 1-maxtree, 2-tree of shapes
	NodeMTPtr root;
	int numNodes;
	std::vector<NodeMTPtr> indexToNode;
	std::vector<NodeMTPtr> nodes;
	AdjacencyRelationPtr adj;
	int depth;
	
	void reconstruction(NodeMTPtr node, int* imgOut);
	void computerTreeAttribute();

public:
   	static const int MAX_TREE = 0;
	static const int MIN_TREE = 1;
	static const int TREE_OF_SHAPES = 2;

	MorphologicalTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	MorphologicalTree(int* img, int numRows, int numCols, bool isMaxtree);

	MorphologicalTree(int* img, int numRows, int numCols);

    ~MorphologicalTree();

	int* getInputImage();
	
	NodeMTPtr getRoot();

	bool isMaxtree();

	int getTreeType();

	NodeMTPtr getSC(int pixel);

    NodeMTPtr getNodeByIndex(int index);

	std::vector<NodeMTPtr>& getIndexNode();

	int getNumNodes();

	int getNumRowsOfImage();

	int getNumColsOfImage();

	int* reconstructionImage();

	int* getImageAferPruning(NodeMTPtr node);

	void pruning(NodeMTPtr node);

	bool isAncestor(NodeMTPtr u, NodeMTPtr v);
	
	bool isDescendant(NodeMTPtr u, NodeMTPtr v);

	bool isComparable(NodeMTPtr u, NodeMTPtr v);

	bool isStrictAncestor(NodeMTPtr u, NodeMTPtr v);
	
	bool isStrictDescendant(NodeMTPtr u, NodeMTPtr v);
	
	bool isStrictComparable(NodeMTPtr u, NodeMTPtr v);
	
	NodeMTPtr findLowestCommonAncestor(NodeMTPtr u, NodeMTPtr v);

	int getDepth();

	std::vector<std::vector<NodeMTPtr>> getNodesByDepth();

	static void extractDepthMap(NodeMTPtr node, int depth, std::vector<std::vector<NodeMTPtr>>& nodesByDepth){
		nodesByDepth[depth].push_back(node);
		for (NodeMTPtr child : node->getChildren()) {
		  extractDepthMap(child, depth + 1, nodesByDepth);
		}
	}
};




class LCAEulerRMQ {
private:
    std::vector<int> euler;            // timePreOrder dos nós na ordem de visita
    std::vector<int> depth;            // profundidade associada a cada posição em euler
    std::vector<int> firstOccurrence;  // [timePreOrder] = posição no vetor euler
    std::vector<std::vector<int>> st;  // Sparse Table para RMQ
    std::vector<NodeMTPtr> indexToNode;  // acesso direto aos nós via timePreOrder

public:
    LCAEulerRMQ(const MorphologicalTreePtr& tree) {
        indexToNode = tree->getIndexNode(); // indexado por timePreOrder
        int n = indexToNode.size();
        firstOccurrence.resize(n, -1);

        dfs(tree->getRoot(), 0);
        buildSparseTable();
    }

	NodeMTPtr findLowestCommonAncestor(const NodeMTPtr& u, const NodeMTPtr& v) {
        int uTime = u->getIndex();
        int vTime = v->getIndex();
        int i = firstOccurrence[uTime];
        int j = firstOccurrence[vTime];
        if (i > j) std::swap(i, j);
        int idx = rmq(i, j);
        return indexToNode[euler[idx]];
    }


private:
    void dfs(const NodeMTPtr& node, int d) {
        int time = node->getIndex();
        if (firstOccurrence[time] == -1)
            firstOccurrence[time] = euler.size();

        euler.push_back(time);
        depth.push_back(d);

        for (const auto& child : node->getChildren()) {
            dfs(child, d + 1);
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