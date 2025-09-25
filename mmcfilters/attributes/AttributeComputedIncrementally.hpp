#pragma once

#include "../utils/Common.hpp"
#include "../attributes/AttributeNames.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {



//Forward declaration
class AttributeComputer; 

using DependencyMap = std::unordered_map<Attribute, std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>>>;

/**
 * @brief Funções utilitárias para computar atributos em árvores de forma incremental.
 *
 * A classe expõe algoritmos genéricos de travessia pós-ordem permitindo
 * compor etapas de pré-processamento, mesclagem e pós-processamento sem criar
 * estruturas auxiliares temporárias. 
 */
class AttributeComputedIncrementally {
public:
    
    template<class PreProcessing, class MergeProcessing, class PostProcessing>
    static void computerAttribute(MorphologicalTree* tree, NodeId root, 
                                    PreProcessing&& preProcessing, 
                                    MergeProcessing&& mergeProcessing, 
                                    PostProcessing&& postProcessing) {
        preProcessing(root);
        for (NodeId child : tree->getChildrenById(root)) {
            AttributeComputedIncrementally::computerAttribute(tree, child, preProcessing, mergeProcessing, postProcessing); // passar por ref (sem cópia)
            mergeProcessing(root, child);
        }
        postProcessing(root);
    }

    static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributesByComputer(MorphologicalTree* tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& available = {});
    static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributesByComputer(MorphologicalTreePtr tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& available = {}){ return computeAttributesByComputer(tree.get(), comp, available); }
	
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTree* tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {});
    static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {}){ return computeSingleAttribute(tree.get(), attr, availableDeps); }
	
	static std::pair<std::shared_ptr<AttributeNamesWithDelta>, std::shared_ptr<float[]>> computeSingleAttributeWithDelta(MorphologicalTree* tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={});
    static std::pair<std::shared_ptr<AttributeNamesWithDelta>, std::shared_ptr<float[]>> computeSingleAttributeWithDelta(MorphologicalTreePtr tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={}){ return computeSingleAttributeWithDelta(tree.get(), attribute, delta, padding, availableDeps); }

	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTree* tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={});
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={}){ return computeAttributes(tree.get(), attributes, providedDependencies); }

	static ImageFloatPtr computerAttributeMapping(MorphologicalTree* tree, Attribute attribute);
    static ImageFloatPtr computerAttributeMapping(MorphologicalTreePtr tree, Attribute attribute){ return computerAttributeMapping(tree.get(), attribute); }
    
};

} // namespace mmcfilters

