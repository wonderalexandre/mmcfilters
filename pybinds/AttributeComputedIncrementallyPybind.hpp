#pragma once

#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"

#include "MorphologicalTreePybind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>
namespace mmcfilters {

class AttributeComputedIncrementallyPybind;
using AttributeComputedIncrementallyPybindPtr = std::shared_ptr<AttributeComputedIncrementallyPybind>;

/**
 * @brief Camada Pybind que expõe utilitários incrementais de atributos.
 */
class AttributeComputedIncrementallyPybind : public AttributeComputedIncrementally{

    public:
    using AttributeComputedIncrementally::AttributeComputedIncrementally;

	

	static std::string describeAttribute(Attribute attribute) {
		return AttributeNames::describe(attribute);
	}

	/*
	static std::vector<std::tuple<NodeMT, NodeMT, float>> extractionExtinctionValues(MorphologicalTreePybindPtr tree, py::array_t<float>& attr){

			std::vector<std::tuple<NodeMT, NodeMT, float>> extinctionValues;
			std::shared_ptr<float[]> attribute = PybindUtils::toShared_ptr(attr);
			auto extValuesPtr = AttributeComputedIncrementally::getExtinctionValue(tree, attribute);
			extinctionValues.reserve(extValuesPtr.size());
			for (const auto& extValue : extValuesPtr) {
					extinctionValues.push_back(std::make_tuple(tree->proxy(extValue->leaf), tree->proxy(extValue->cutoffNode), extValue->extinction));
			}
			return extinctionValues;
	}*/

	static py::array_t<float> computeSingleAttribute(MorphologicalTreePybindPtr tree, Attribute attribute){
		auto [attributeNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(tree, attribute);
		return PybindUtils::toNumpyShared_ptr(buffer, tree->getNumNodes());
	}

	static std::pair<py::dict, py::array_t<float>> computeSingleAttributeWithDelta(MorphologicalTreePybindPtr tree, Attribute attribute, int delta, std::string padding = "null-padding") {
		auto [attributeNames, buffer] = AttributeComputedIncrementally::computeSingleAttributeWithDelta(tree, attribute, delta, padding, /*deps*/{});

		const int numAttribute = attributeNames->NUM_ATTRIBUTES;
		const int n = tree->getNumNodes();

		// Extrai nomes e offsets das chaves compostas
		std::vector<std::string> keys;
		std::vector<int> values;

		for (const auto& pair : attributeNames->indexMap) {
			const AttributeKey& attrKey = pair.first;
			int offset = pair.second;
			keys.push_back(AttributeNamesWithDelta::toString(attrKey.attr, attrKey.delta));
			values.push_back(offset);
		}

		// Ordena pelo offset, garantindo alinhamento com as colunas do array
		std::vector<size_t> indices(values.size());
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&values](size_t i1, size_t i2) { return values[i1] < values[i2]; });

		py::dict dict;
		for (size_t i = 0; i < indices.size(); ++i) {
			dict[py::str(keys[indices[i]])] = values[indices[i]];
		}

		// Prepara o capsule para o numpy array
		py::capsule free_when_done(new std::shared_ptr<float[]>(buffer), [](void* ptr) {
			delete reinterpret_cast<std::shared_ptr<float[]>*>(ptr);
		});

		// Atenção para strides: (linha, coluna) = (num_attrs*sizeof(float), sizeof(float))
		py::array_t<float> numpy = py::array(py::buffer_info(
			buffer.get(),
			sizeof(float),
			py::format_descriptor<float>::value,
			2,
			{ n, numAttribute },
			{ sizeof(float) * numAttribute, sizeof(float) }
		), free_when_done);

		return std::make_pair(dict, numpy);
	}

	static py::array_t<float> computerAttributeMapping(MorphologicalTreePybindPtr tree, Attribute attribute) {
		auto imgFloatPtr = AttributeComputedIncrementally::computerAttributeMapping(tree, attribute);
		return PybindUtils::toNumpy(imgFloatPtr);
	} 

	static std::pair<py::dict, py::array_t<float>> computeAttributesFromList(MorphologicalTreePybindPtr tree, const std::vector<AttributeOrGroup>& attributes) {
		auto [attributeNames, buffer] = AttributeComputedIncrementally::computeAttributes(tree, attributes);
		// Cria uma cópia do shared_ptr para manter ownership no Python
		std::shared_ptr<float[]> bufferCopy = buffer;
		
		const int numAttribute = attributeNames->NUM_ATTRIBUTES;
		const int n = tree->getNumNodes();
        
		std::vector<std::string> keys;
		std::vector<int> values;

		// 1. Copiar chaves e valores para vetores separados
		for (const auto& pair : attributeNames->indexMap) {
			Attribute attribute = pair.first;
			int offset = pair.second;

			keys.push_back( attributeNames->toString(attribute) );
			values.push_back(offset);
		}

		// 2. Criar um vetor de índices para ordenar os valores
		std::vector<size_t> indices(values.size());
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&values](size_t i1, size_t i2) { return values[i1] < values[i2]; });

		py::dict dict;
		for (size_t i = 0; i < indices.size(); ++i) {
			dict[py::str( keys[indices[i]] )] = values[indices[i]];
		}

		
		py::capsule free_when_done(new std::shared_ptr<float[]>(bufferCopy), [](void* ptr) {
			// Converte de volta e destrói corretamente
			delete reinterpret_cast<std::shared_ptr<float[]>*>(ptr);
		});
		
		py::array_t<float> numpy = py::array(py::buffer_info(
			buffer.get(),
			sizeof(float),
			py::format_descriptor<float>::value,
			2,
			{  n,  numAttribute }, 
			{ sizeof(float)*numAttribute, sizeof(float) }
		), free_when_done);
		
		
		return std::make_pair(dict, numpy);
		
	}

};

} // namespace mmcfilters
