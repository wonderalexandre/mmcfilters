#ifndef ATTRIBUTE_COMPUTED_INCREMENTALLY_PYBIND_H
#define ATTRIBUTE_COMPUTED_INCREMENTALLY_PYBIND_H


#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/NodeCT.hpp"

#include "../pybind/ComponentTreePybind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>

class AttributeComputedIncrementallyPybind : public AttributeComputedIncrementally{

    public:
    using AttributeComputedIncrementally::AttributeComputedIncrementally;



	static py::array_t<float> computerArea(ComponentTreePybind *tree){
		const int n = tree->getNumNodes();
		float *attrs = new float[n];
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&attrs](NodeCT* node) -> void { //pre-processing
				attrs[node->getIndex()] = node->getCNPs().size(); //area
			},
			[&attrs](NodeCT* parent, NodeCT* child) -> void { //merge-processing
				attrs[parent->getIndex()] += attrs[child->getIndex()];
			},
			[](NodeCT* node) -> void { //post-processing			
		});
	    py::array_t<float> numpy = py::array(py::buffer_info(
			attrs,            
			sizeof(float),     
			py::format_descriptor<float>::value, 
			1,         
			{  n }, 
			{ sizeof(float) }
	    ));
		return numpy;
	}

    static std::pair<py::dict, py::array_t<float>> computerBasicAttributes(ComponentTreePybind *tree){
        const int numAttribute = 19;
		const int n = tree->getNumNodes();
        float* attrs = AttributeComputedIncrementally::computerBasicAttributes(tree);
		
		py::dict dict;
        dict[py::str("AREA")] = 0;
		dict[py::str("VOLUME")] = 1;
		dict[py::str("LEVEL")] = 2;
		dict[py::str("MEAN_LEVEL")] = 3;
		dict[py::str("VARIANCE_LEVEL")] = 4;
		dict[py::str("WIDTH")] = 5;
		dict[py::str("HEIGHT")] = 6;
		dict[py::str("RETANGULARITY")] = 7;
		dict[py::str("RATIO_WH")] = 8;
		dict[py::str("MOMENT_02")] = 9;
		dict[py::str("MOMENT_20")] = 10;
		dict[py::str("MOMENT_11")] = 11;
		dict[py::str("INERTIA")] = 12;
		dict[py::str("ORIENTATION")] = 13;
		dict[py::str("LEN_MAJOR_AXIS")] = 14;
		dict[py::str("LEN_MINOR_AXIS")] = 15;
		dict[py::str("ECCENTRICITY")] = 16;
		dict[py::str("COMPACTNESS")] = 17;
		dict[py::str("STD_LEVEL")] = 18;
		
	    py::array_t<float> numpy = py::array(py::buffer_info(
			attrs,            
			sizeof(float),     
			py::format_descriptor<float>::value, 
			2,         
			{  n,  numAttribute }, 
			{ sizeof(float), sizeof(float) * n }
	    ));
		
		return std::make_pair(dict, numpy);
	}
};

#endif 