#ifndef ATTRIBUTE_COMPUTED_INCREMENTALLY_PYBIND_H
#define ATTRIBUTE_COMPUTED_INCREMENTALLY_PYBIND_H


#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/NodeMT.hpp"

#include "../pybind/MorphologicalTreePybind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>

class AttributeComputedIncrementallyPybind : public AttributeComputedIncrementally{

    public:
    using AttributeComputedIncrementally::AttributeComputedIncrementally;


	static py::array_t<float> computerArea(MorphologicalTreePybindPtr tree){
		const int n = tree->getNumNodes();
		float *attrs = new float[n];
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&attrs](NodeMTPtr node) -> void { //pre-processing
				attrs[node->getIndex()] = node->getCNPs().size(); //area
			},
			[&attrs](NodeMTPtr parent, NodeMTPtr child) -> void { //merge-processing
				attrs[parent->getIndex()] += attrs[child->getIndex()];
			},
			[](NodeMTPtr node) -> void { //post-processing			
		});
	    
		py::capsule free_when_done(attrs, [](void *f) {
			delete[] reinterpret_cast<float*>(f);
		});
		
		py::array_t<float> numpy = py::array(py::buffer_info(
			attrs,
			sizeof(float),
			py::format_descriptor<float>::value,
			1,
			{ n },
			{ sizeof(float) }
		), free_when_done);

		return numpy;
	}

	static py::dict extractCountors(MorphologicalTreePybindPtr tree) {
		auto contours = AttributeComputedIncrementally::extractCountors(tree);  // chama o método original
	
		py::dict pyContours;
		for (size_t nodeIdx = 0; nodeIdx < contours.size(); ++nodeIdx) {
			py::set pySet;
			for (int pixel : contours[nodeIdx]) {
				pySet.add(pixel);
			}
			pyContours[py::int_(nodeIdx)] = pySet;
		}
	
		return pyContours;
	}
	

    static std::pair<py::dict, py::array_t<float>> computerBasicAttributes(MorphologicalTreePybindPtr tree){
        
		auto [attributeNames, ptrValues] = AttributeComputedIncrementally::computerBasicAttributes(tree);
		const int numAttribute = attributeNames.NUM_ATTRIBUTES;
		const int n = tree->getNumNodes();
        
		
		std::vector<std::string> keys;
		std::vector<int> values;

		// 1. Copiar chaves e valores para vetores separados
		for (const auto& pair : attributeNames.mapIndexes) {
			keys.push_back(pair.first);
			values.push_back(pair.second);
		}

		// 2. Criar um vetor de índices para ordenar os valores
		std::vector<size_t> indices(values.size());
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&values](size_t i1, size_t i2) { return values[i1] < values[i2]; });

		py::dict dict;
		for (size_t i = 0; i < indices.size(); ++i) {
			dict[py::str( keys[indices[i]] )] = values[indices[i]] / n; 
		}

		py::capsule free_when_done(ptrValues, [](void *f) {
			delete[] reinterpret_cast<float*>(f);
		});
		
		py::array_t<float> numpy = py::array(py::buffer_info(
			ptrValues,
			sizeof(float),
			py::format_descriptor<float>::value,
			2,
			{  n,  numAttribute }, 
			{ sizeof(float), sizeof(float) * n }
		), free_when_done);

		
		return std::make_pair(dict, numpy);
	}
};

#endif 