
#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include <vector>

#ifndef COMPUTER_MSER_H
#define COMPUTER_MSER_H


#define UNDEF -999999999999


class ComputerMSER {
private:
	
	MorphologicalTreePtr tree;
	double maxVariation;
	int minArea;
	int maxArea;
	int num;
	std::vector<NodeMTPtr> ascendants;
	std::vector<NodeMTPtr> descendants;
	std::vector<double> stability;
	
	double* attr_mser;
	double* attr_area;

	NodeMTPtr getNodeAscendant(NodeMTPtr node, int h);

	void maxAreaDescendants(NodeMTPtr nodeAsc, NodeMTPtr nodeDes);
	
public:
	ComputerMSER(MorphologicalTreePtr tree);
	ComputerMSER(MorphologicalTreePtr tree, double* attr_increasing);
	
	~ComputerMSER();

	std::vector<bool> computerMSER(int delta);

	int getNumNodes();

	NodeMTPtr descendantWithMaxStability(NodeMTPtr node);
	
	NodeMTPtr ascendantWithMaxStability(NodeMTPtr node);

	std::vector<double> getStabilities();

	double getStability(NodeMTPtr node);

	std::vector<NodeMTPtr> getAscendants();

	std::vector<NodeMTPtr> getDescendants();

	NodeMTPtr getAscendant(NodeMTPtr node);
	
	NodeMTPtr getDescendant(NodeMTPtr node);

	NodeMTPtr getNodeInPathWithMaxStability(NodeMTPtr node, std::vector<bool> isMSER);

	void setMaxVariation(double maxVariation);
	void setMinArea(int minArea);
	void setMaxArea(int maxArea);
	
};

#endif