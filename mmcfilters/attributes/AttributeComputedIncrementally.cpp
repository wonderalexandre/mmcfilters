
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "AttributeFactory.hpp"

namespace mmcfilters {

ImageFloatPtr AttributeComputedIncrementally::computerAttributeMapping(MorphologicalTree* tree, Attribute attribute) {
    auto [attrNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(tree, attribute);
    ImageFloatPtr imgPtr = std::make_shared<ImageFloat>(tree->getNumRowsOfImage(), tree->getNumColsOfImage());
    float* img = imgPtr->rawData();
    for(int p=0; p < imgPtr->getSize(); ++p){
        int index = tree->getSCById(p);
        img[p] = buffer[attrNames->linearIndex(index, attribute)];
    }
    return imgPtr;
}


std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> AttributeComputedIncrementally::computeAttributesByComputer(MorphologicalTree* tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& availableDeps) {
    // Lambda para obter os atributos de um grupo
    auto attributesOf = [](AttributeGroup group) -> std::vector<Attribute> {
        const auto it = ATTRIBUTE_GROUPS.find(group);
        if (it != ATTRIBUTE_GROUPS.end()) return it->second;
        throw std::runtime_error("Unknown AttributeGroup in attributesOf.");
    };

    DependencyMap available = availableDeps;

    // Resolver dependências automaticamente
    for (const AttributeOrGroup& dep : comp->requiredAttributes()) {
        bool needsCompute = false;

        if (std::holds_alternative<Attribute>(dep)) {
            needsCompute = !available.count(std::get<Attribute>(dep));
        } else {
            const auto& groupAttrs = attributesOf(std::get<AttributeGroup>(dep));
            for (const auto& attr : groupAttrs) {
                if (!available.count(attr)) {
                    needsCompute = true;
                    break;
                }
            }
        }

        if (needsCompute) {
            auto [depNames, depBuf] = computeSingleAttribute(tree, dep, available);
            for (const auto& [a, _] : depNames->indexMap) {
                available[a] = {depNames, depBuf};
            }
        }
    }

    // Coletar dependências em ordem
    std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>> dependencySources;
    for (const auto& dep : comp->requiredAttributes()) {
        std::vector<Attribute> requiredAttrs;
        if (std::holds_alternative<Attribute>(dep)) {
            requiredAttrs.push_back(std::get<Attribute>(dep));
        } else {
            auto groupAttrs = attributesOf(std::get<AttributeGroup>(dep));
            requiredAttrs.insert(requiredAttrs.end(), groupAttrs.begin(), groupAttrs.end());
        }

        for (const auto& attr : requiredAttrs) {
            dependencySources.push_back(available.at(attr));
        }
    }

    // Todos os atributos produzidos por esse computador
    const auto& computedAttrs = comp->attributes();

    // Construir o objeto AttributeNames
    std::unordered_map<Attribute, int> attrOffsets;
    for (int i = 0; i < static_cast<int>(computedAttrs.size()); ++i) {
        attrOffsets[computedAttrs[i]] = i;
    }
    auto attrNames = std::make_shared<AttributeNames>(std::move(attrOffsets));

    // Alocar buffer
    int n = tree->getNumNodes();
    std::shared_ptr<float[]> buffer(new float[n * attrNames->NUM_ATTRIBUTES]());

    // Computar atributos (modo completo)
    comp->compute(tree, buffer, attrNames, dependencySources);

    return {attrNames, buffer};
}

std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> AttributeComputedIncrementally::computeSingleAttribute(MorphologicalTree* tree, AttributeOrGroup attrOrGroup, const DependencyMap& availableDeps) {

    // Lambda para obter os atributos de um grupo
    auto attributesOf = [](AttributeGroup group) -> std::vector<Attribute> {
        const auto it = ATTRIBUTE_GROUPS.find(group);
        if (it != ATTRIBUTE_GROUPS.end()) return it->second;
        throw std::runtime_error("Unknown AttributeGroup in attributesOf.");
    };

    auto comp = AttributeFactory::create(attrOrGroup);
    DependencyMap available = availableDeps;

    // Descobrir quais atributos estão sendo solicitados
    std::vector<Attribute> requestedAttrs;
    if (std::holds_alternative<Attribute>(attrOrGroup)) {
        requestedAttrs.push_back(std::get<Attribute>(attrOrGroup));
    } else {
        requestedAttrs = attributesOf(std::get<AttributeGroup>(attrOrGroup));
    }

    // Resolver dependências automaticamente
    for (const AttributeOrGroup& dep : comp->requiredAttributes()) {
        // Verifica se todos os atributos do grupo (ou atributo individual) já estão disponíveis
        bool needsCompute = false;

        if (std::holds_alternative<Attribute>(dep)) {
            needsCompute = !available.count(std::get<Attribute>(dep));
        } else {
            const auto& groupAttrs = attributesOf(std::get<AttributeGroup>(dep));
            for (const auto& attr : groupAttrs) {
                if (!available.count(attr)) {
                    needsCompute = true;
                    break;
                }
            }
        }

        if (needsCompute) {
            auto [depNames, depBuf] = computeSingleAttribute(tree, dep, available);
            for (const auto& [a, _] : depNames->indexMap) {
                available[a] = {depNames, depBuf};
            }
        }
    }

    // Coletar dependências
    std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>> dependencySources;
    for (const AttributeOrGroup& dep : comp->requiredAttributes()) {
        std::vector<Attribute> requiredList;
        if (std::holds_alternative<Attribute>(dep)) {
            requiredList.push_back(std::get<Attribute>(dep));
        } else {
            const auto& groupAttrs = attributesOf(std::get<AttributeGroup>(dep));
            requiredList.insert(requiredList.end(), groupAttrs.begin(), groupAttrs.end());
        }

        for (const Attribute& attr : requiredList) {
            dependencySources.push_back(available.at(attr));
        }
    }

    // Construir AttributeNames com apenas os atributos solicitados
    std::unordered_map<Attribute, int> attrOffsets;
    for (int i = 0; i < static_cast<int>(requestedAttrs.size()); ++i) {
        attrOffsets[requestedAttrs[i]] = i;
    }
    auto attrNames = std::make_shared<AttributeNames>(std::move(attrOffsets));

    // Alocar buffer de saída
    int n = tree->getNumNodes();
    std::shared_ptr<float[]> buffer(new float[n * attrNames->NUM_ATTRIBUTES]());

    // Computar atributos
    comp->compute(tree, buffer, attrNames, requestedAttrs, dependencySources);

    return {attrNames, buffer};
}



std::pair<std::shared_ptr<AttributeNamesWithDelta>, std::shared_ptr<float[]>> AttributeComputedIncrementally::computeSingleAttributeWithDelta(MorphologicalTree* tree, Attribute attribute, int delta, std::string padding, const DependencyMap& availableDeps) {
	/*
		Valores de padding:
		 - zero-padding: preenchimento com zero
		 - nan-padding: preenchimento com o valor NaN 
		 - last-padding: preenchimento com o ultimo valor valido
		 - null-padding: preenchimento com 0 todo os nos do caminho
	*/

    // 1. Computa atributo base (delta = 0)
    auto [attributeNamesBase, attrsBase] = AttributeComputedIncrementally::computeSingleAttribute(tree, attribute, availableDeps);

    int n = tree->getNumNodes();
    std::vector<Attribute> attrVec = {attribute};
    auto attributeNamesDelta = std::make_shared<AttributeNamesWithDelta>(
        AttributeNamesWithDelta::create(delta, attrVec)
    );

    // Aloca buffer do novo atributo delta (inicializado com zero)
    std::shared_ptr<float[]> attrsDelta(new float[n * attributeNamesDelta->NUM_ATTRIBUTES]());

    // Para d = 0, copia valor do próprio nó SEMPRE
    for (NodeId nodeIndex : tree->getNodeIds()) {
        int outIdx = attributeNamesDelta->linearIndex(nodeIndex, attribute, 0);
        int baseIdx = attributeNamesBase->linearIndex(nodeIndex, attribute);
        attrsDelta[outIdx] = attrsBase[baseIdx];
    }

    // Para d > 0, só copia se realmente existe ascendente/descendente (zero-padding padrão)
    for (int d = 1; d <= delta; ++d) {
        auto [ascendants, descendants] = tree->computerAscendantsAndDescendants(d);
        
        for (NodeId nodeIndex : tree->getNodeIds()) {
            
            // Ascendente (-d)
            int ascIndex = (ascendants[nodeIndex] != InvalidNode ? ascendants[nodeIndex] : nodeIndex);
            if (ascIndex != nodeIndex) {
                int outIdxAsc = attributeNamesDelta->linearIndex(nodeIndex, attribute, -d);
                int baseIdxAsc = attributeNamesBase->linearIndex(ascIndex, attribute);
                attrsDelta[outIdxAsc] = attrsBase[baseIdxAsc];
            }
            // Descendente (+d)
            int descIndex = (descendants[nodeIndex] != InvalidNode ? descendants[nodeIndex] : nodeIndex);
            if (descIndex != nodeIndex) {
                int outIdxDesc = attributeNamesDelta->linearIndex(nodeIndex, attribute, +d);
                int baseIdxDesc = attributeNamesBase->linearIndex(descIndex, attribute);
                attrsDelta[outIdxDesc] = attrsBase[baseIdxDesc];
            }
        }
    }

    // 4. Preenche delta > 0 (só copia se realmente existe ascendente/descendente)
    for (int d = 1; d <= delta; ++d) {
        auto [ascendants, descendants] = tree->computerAscendantsAndDescendants(d);
        
        for (NodeId nodeIndex : tree->getNodeIds()) {

            // Ascendente (-d)
            int ascIndex = (ascendants[nodeIndex]!=InvalidNode ? ascendants[nodeIndex] : nodeIndex);
            if (ascIndex != nodeIndex) {
                int outIdxAsc = attributeNamesDelta->linearIndex(nodeIndex, attribute, -d);
                int baseIdxAsc = attributeNamesBase->linearIndex(ascIndex, attribute);
                attrsDelta[outIdxAsc] = attrsBase[baseIdxAsc];
            }

            // Descendente (+d)
            int descIndex = (descendants[nodeIndex] !=InvalidNode ? descendants[nodeIndex] : nodeIndex);
            if (descIndex != nodeIndex) {
                int outIdxDesc = attributeNamesDelta->linearIndex(nodeIndex, attribute, +d);
                int baseIdxDesc = attributeNamesBase->linearIndex(descIndex, attribute);
                attrsDelta[outIdxDesc] = attrsBase[baseIdxDesc];
            }
        }
    }

    // 5. Aplica padding para modos diferentes de zero-padding
    if (padding != "zero-padding") {
        for (NodeId nodeIndex : tree->getNodeIds()) {

            // Ascendente (-d)
            for (int d = 1; d <= delta; ++d) {
                int outIdx = attributeNamesDelta->linearIndex(nodeIndex, attribute, -d);
                int refIdx = attributeNamesDelta->linearIndex(nodeIndex, attribute, -(d - 1));

                if (attrsDelta[outIdx] == 0) {
                    if (padding == "last-padding") {
                        attrsDelta[outIdx] = attrsDelta[refIdx];
                    } else if (padding == "nan-padding") {
                        attrsDelta[outIdx] = std::numeric_limits<float>::quiet_NaN();
                    } else if (padding == "null-padding") {
                        for (int k = 0; k <= d; ++k) {
                            // Marca todos os deltas negativos e positivos até d como NaN
                            attrsDelta[attributeNamesDelta->linearIndex(nodeIndex, attribute, -k)] = std::numeric_limits<float>::quiet_NaN();
                            attrsDelta[attributeNamesDelta->linearIndex(nodeIndex, attribute, +k)] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                }
            }

            // Descendente (+d)
            for (int d = 1; d <= delta; ++d) {
                int outIdx = attributeNamesDelta->linearIndex(nodeIndex, attribute, +d);
                int refIdx = attributeNamesDelta->linearIndex(nodeIndex, attribute, +(d - 1));

                if (tree->isLeafById(nodeIndex) || attrsDelta[outIdx] == 0) {
                    if (padding == "last-padding") {
                        attrsDelta[outIdx] = attrsDelta[refIdx];
                    } else if (padding == "nan-padding") {
                        attrsDelta[outIdx] = std::numeric_limits<float>::quiet_NaN();
                    }
                }
                if (attrsDelta[outIdx] == 0 && padding == "null-padding") {
                    for (int k = 0; k <= d; ++k) {
                        attrsDelta[attributeNamesDelta->linearIndex(nodeIndex, attribute, +k)] = std::numeric_limits<float>::quiet_NaN();
                        attrsDelta[attributeNamesDelta->linearIndex(nodeIndex, attribute, -k)] = std::numeric_limits<float>::quiet_NaN();
                    }
                }
            }
        }
    }

    return {attributeNamesDelta, attrsDelta};
}






static std::vector<std::shared_ptr<AttributeComputer>> getOrderedComputers(const std::vector<AttributeOrGroup>& attrOrGroups) {
    using ACptr = std::shared_ptr<AttributeComputer>;
    using TIndex = std::type_index;

    auto attributesOf = [](AttributeGroup group) -> std::vector<Attribute> {
        const auto it = ATTRIBUTE_GROUPS.find(group);
        if (it != ATTRIBUTE_GROUPS.end()) return it->second;
        throw std::runtime_error("Unknown AttributeGroup in attributesOf.");
    };

    std::map<TIndex, ACptr> computerMap;
    std::map<TIndex, std::set<TIndex>> dependencyGraph;
    std::set<Attribute> visitedAttrs;

    std::function<void(Attribute)> collect = [&](Attribute attr) {
        if (visitedAttrs.count(attr)) return;
        visitedAttrs.insert(attr);

        auto comp = AttributeFactory::create(attr);
        const auto& ref = *comp;
        TIndex id(typeid(ref));

        computerMap[id] = comp;
        dependencyGraph[id];

        for (const auto& depOrGroup : comp->requiredAttributes()) {
            std::vector<Attribute> deps;
            if (std::holds_alternative<Attribute>(depOrGroup)) {
                deps.push_back(std::get<Attribute>(depOrGroup));
            } else {
                auto groupAttrs = attributesOf(std::get<AttributeGroup>(depOrGroup));
                deps.insert(deps.end(), groupAttrs.begin(), groupAttrs.end());
            }

            for (const auto& depAttr : deps) {
                collect(depAttr);
                auto depComp = AttributeFactory::create(depAttr);
                const auto& depRef = *depComp;
                std::type_index depId(typeid(depRef));
                dependencyGraph[id].insert(depId);
            }
        }
    };

    for (const auto& item : attrOrGroups) {
        std::vector<Attribute> attrs;
        if (std::holds_alternative<Attribute>(item)) {
            attrs.push_back(std::get<Attribute>(item));
        } else {
            auto groupAttrs = attributesOf(std::get<AttributeGroup>(item));
            attrs.insert(attrs.end(), groupAttrs.begin(), groupAttrs.end());
        }

        for (const auto& attr : attrs) {
            collect(attr);
        }
    }

    std::vector<ACptr> ordered;
    std::set<TIndex> visited;

    std::function<void(TIndex)> visit = [&](TIndex id) {
        if (visited.count(id)) return;
        visited.insert(id);
        for (const auto& depId : dependencyGraph[id]) {
            visit(depId);
        }
        ordered.push_back(computerMap.at(id));
    };

    for (const auto& [id, _] : computerMap) {
        visit(id);
    }

    return ordered;
}

std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> AttributeComputedIncrementally::computeAttributes(MorphologicalTree* tree, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& providedDependencies) {
    DependencyMap available = providedDependencies;

    auto attributesOf = [](AttributeGroup group) -> std::vector<Attribute> {
        const auto it = ATTRIBUTE_GROUPS.find(group);
        if (it != ATTRIBUTE_GROUPS.end()) return it->second;
        throw std::runtime_error("Unknown AttributeGroup.");
    };

    // Etapa 1: Expansão dos atributos
    std::set<Attribute> uniqueExpandedAttrs;
    std::unordered_map<std::type_index, std::vector<Attribute>> attributesPerComputer;

    for (const auto& item : attributes) {
        std::vector<Attribute> attrs;
        if (std::holds_alternative<Attribute>(item)) {
            attrs.push_back(std::get<Attribute>(item));
        } else {
            auto groupAttrs = attributesOf(std::get<AttributeGroup>(item));
            attrs.insert(attrs.end(), groupAttrs.begin(), groupAttrs.end());
        }

        for (const auto& attr : attrs) {
            uniqueExpandedAttrs.insert(attr);
            auto comp = AttributeFactory::create(attr);
            const auto& ref = *comp;
            std::type_index id(typeid(ref));
            attributesPerComputer[id].push_back(attr);
        }
    }

    // Etapa 2: Construção do buffer principal (apenas atributos solicitados)
    std::unordered_map<Attribute, int> attrOffsets;
    int offset = 0;
    for (const auto& attr : uniqueExpandedAttrs) {
        attrOffsets[attr] = offset++;
    }

    auto attrNames = std::make_shared<AttributeNames>(std::move(attrOffsets));
    int n = tree->getNumNodes();
    std::shared_ptr<float[]> buffer(new float[n * attrNames->NUM_ATTRIBUTES]());

    // Etapa 3: Obter ordem de execução
    auto orderedComputers = getOrderedComputers(attributes);

    // Etapa 4: Executar os computadores em ordem
    for (const auto& comp : orderedComputers) {
        const auto& ref = *comp;
        std::type_index id(typeid(ref));
        
        if (!attributesPerComputer.count(id)) { //Se o computador não tem atributos solicitados
            auto [tempNames, tempBuffer] = computeAttributesByComputer(tree, comp, available);
            for (const auto& attrTemp : comp->attributes()) {
                available[attrTemp] = {tempNames, tempBuffer};
            }
            continue;
        }

        // Atributos solicitados para esse computador        
        std::vector<Attribute> userRequestedAttrs  = attributesPerComputer.at(id); 
        
        // Verifica se os atributos solicitados já foram computados anteirmente. Se todos os atributos já estão disponíveis, não precisa computar novamente
        bool alreadyAvailable = std::all_of(userRequestedAttrs.begin(), userRequestedAttrs.end(), [&](const Attribute& a) { return available.count(a); });
        if (alreadyAvailable) continue;

        // Resolver dependências: 
        std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>> depsForThis;
        for (const auto& req : comp->requiredAttributes()) {
            std::vector<Attribute> deps;
            if (std::holds_alternative<Attribute>(req)) {
                deps.push_back(std::get<Attribute>(req));
            } else {
                auto groupAttrs = attributesOf(std::get<AttributeGroup>(req));
                deps.insert(deps.end(), groupAttrs.begin(), groupAttrs.end());
            }

            for (const auto& dep : deps) {
                if (!available.count(dep)) {
                    auto [depNames, depBuf] = computeSingleAttribute(tree, dep, available);
                    for (const auto& [a, _] : depNames->indexMap) {
                        available[a] = {depNames, depBuf};
                    }
                }
                depsForThis.push_back(available.at(dep));
            }
        }

        // Computa os atributos solicitados no buffer principal
        comp->compute(tree, buffer, attrNames, userRequestedAttrs, depsForThis);

        // Marca os atributos solicitados como disponíveis
        for (const auto& attr : userRequestedAttrs) {
            available[attr] = {attrNames, buffer};
        }
    }

    return {attrNames, buffer};
}




} // namespace mmcfilters
