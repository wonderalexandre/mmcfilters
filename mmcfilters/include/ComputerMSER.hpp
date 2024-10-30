
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include <vector>

#ifndef COMPUTER_MSER_H
#define COMPUTER_MSER_H


#define UNDEF -999999999999


class ComputerMSER {
private:
	
	ComponentTree* tree;
	double maxVariation;
	int minArea;
	int maxArea;
	int num;
	std::vector<NodeCT*> ascendants;
	std::vector<NodeCT*> descendants;
	std::vector<double> stability;
	
	double* attr_mser;
	double* attr_area;

	NodeCT* getNodeAscendant(NodeCT* node, int h);

	void maxAreaDescendants(NodeCT* nodeAsc, NodeCT* nodeDes);
	
public:
	ComputerMSER(ComponentTree* tree);
	ComputerMSER(ComponentTree* tree, double* attr_increasing);
	
	~ComputerMSER();

	std::vector<bool> computerMSER(int delta);

	int getNumNodes();

	NodeCT* descendantWithMaxStability(NodeCT* node);
	
	NodeCT* ascendantWithMaxStability(NodeCT* node);

	std::vector<double> getStabilities();

	double getStability(NodeCT* node);

	std::vector<NodeCT*> getAscendants();

	std::vector<NodeCT*> getDescendants();

	NodeCT* getAscendant(NodeCT* node);
	
	NodeCT* getDescendant(NodeCT* node);

	NodeCT* getNodeInPathWithMaxStability(NodeCT* node, std::vector<bool> isMSER);

	void setMaxVariation(double maxVariation);
	void setMinArea(int minArea);
	void setMaxArea(int maxArea);
	
};

#endif