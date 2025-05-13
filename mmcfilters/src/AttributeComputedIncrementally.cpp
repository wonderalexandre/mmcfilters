#include <set>
#include <map>
#include <typeindex>
#include "../include/AttributeComputedIncrementally.hpp"


void AttributeComputedIncrementally::preProcessing(NodeMTPtr v){}

void AttributeComputedIncrementally::mergeChildren(NodeMTPtr parent, NodeMTPtr child){}

void AttributeComputedIncrementally::postProcessing(NodeMTPtr parent){}

void AttributeComputedIncrementally::computerAttribute(NodeMTPtr root) {
        preProcessing(root);
        for (NodeMTPtr child : root->getChildren())
        {
            computerAttribute(child);
            mergeChildren(root, child);
        }
        postProcessing(root);
}


std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> AttributeComputedIncrementally::computeSingleAttribute(MorphologicalTreePtr tree, Attribute attr, const DependencyMap& availableDeps){
    // Cria instância do computador para o atributo solicitado
    auto comp = AttributeFactory::create(attr);

    // Cria um mapa local copiando availableDeps (para permitir recursão e resolução de dependências)
    DependencyMap available = availableDeps;

    // Etapa 1: resolver dependências automaticamente
    for (const Attribute& dep : comp->requiredAttributes()) {
        if (!available.count(dep)) {
            auto [depNames, depBuf] = computeSingleAttribute(tree, dep, available);
            for (const auto& [a, _] : depNames->indexMap) {
                available[a] = {depNames, depBuf};
            }
        }
    }

    // Etapa 2: extrair somente as dependências requeridas (na ordem esperada pela função compute)
    std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>> dependencySources;
    for (const Attribute& dep : comp->requiredAttributes()) {
        dependencySources.push_back(available.at(dep));
    }

    // Etapa 3: construir o AttributeNames com offsets
    std::unordered_map<Attribute, int> attrOffsets;
    const auto& attrs = comp->attributes();
    for (int i = 0; i < attrs.size(); ++i) {
        attrOffsets[attrs[i]] = i;
    }
    auto attrNames = std::make_shared<AttributeNames>(std::move(attrOffsets));

    // Etapa 4: alocar buffer de saída
    int n = tree->getNumNodes();
    std::shared_ptr<float[]> buffer(new float[n * attrNames->NUM_ATTRIBUTES]());

    // Etapa 5: computar os atributos
    comp->compute(tree, buffer, attrNames, dependencySources);

    return {attrNames, buffer};
}

std::vector<std::shared_ptr<AttributeComputer>> getOrderedComputers(const std::vector<Attribute>& attributes) {
    using ACptr = std::shared_ptr<AttributeComputer>;
    using TIndex = std::type_index;

    // Mapa: tipo do computador -> instância compartilhada
    std::map<TIndex, ACptr> computerMap;

    // Grafo de dependência: tipo -> tipos dos requisitos
    std::map<TIndex, std::set<TIndex>> dependencyGraph;

    // Fila de atributos pendentes (inclui transitivos)
    std::set<Attribute> allAttributes(attributes.begin(), attributes.end());

    // Expande os requisitos de forma transitiva
    std::function<void(Attribute)> collect = [&](Attribute attr) {
        auto comp = AttributeFactory::create(attr);
        const auto& ref = *comp;
        TIndex id(typeid(ref));

        if (computerMap.count(id)) return; // já visitado

        computerMap[id] = comp;
        dependencyGraph[id] = {};

        for (const auto& dep : comp->requiredAttributes()) {
            auto depComp = AttributeFactory::create(dep);
            const auto& depRef = *depComp;
            TIndex depId(typeid(depRef));
            dependencyGraph[id].insert(depId);
            allAttributes.insert(dep);
            collect(dep); // recursivo
        }
    };

    for (const auto& attr : attributes) {
        collect(attr);
    }

    // Ordenação topológica
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

std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> AttributeComputedIncrementally::computeAttributes(MorphologicalTreePtr tree, const std::vector<Attribute>& attributes, const DependencyMap& providedDependencies){
    DependencyMap available = providedDependencies;
    std::set<std::type_index> addedComputers;
    std::vector<std::shared_ptr<AttributeComputer>> computers;

    // Etapa 1: criar lista de computadores únicos necessários
    for (const auto& attr : attributes) {
        auto comp = AttributeFactory::create(attr);
        const auto& ref = *comp;
        std::type_index id(typeid(ref));
        if (!addedComputers.count(id)) {
            addedComputers.insert(id);
            computers.push_back(comp);
        }
    }

    // Etapa 2: definir offsets no buffer final
    std::unordered_map<Attribute, int> attrOffsets;
    int offset = 0;
    for (const auto& comp : computers) {
        for (const auto& attr : comp->attributes()) {
            attrOffsets[attr] = offset++;
        }
    }

    auto attrNames = std::make_shared<AttributeNames>(std::move(attrOffsets));
    int n = tree->getNumNodes();
    std::shared_ptr<float[]> buffer(new float[n * attrNames->NUM_ATTRIBUTES]());

    // Etapa 3: computar cada grupo de atributos com dependências resolvidas
    auto orderedComputers = getOrderedComputers(attributes);
    for (const auto& comp : orderedComputers) {
        // Verifica se já estão disponíveis
        bool alreadyAvailable = true;
        for (const auto& attr : comp->attributes()) {
            if (!available.count(attr)) {
                alreadyAvailable = false;
                break;
            }
        }
        if (alreadyAvailable) continue;

        // Resolve dependências recursivamente
        for (const auto& req : comp->requiredAttributes()) {
            if (!available.count(req)) {
                auto [depNames, depBuf] = computeSingleAttribute(tree, req, available);
                for (const auto& [a, _] : depNames->indexMap) {
                    available[a] = {depNames, depBuf};
                }
            }
        }

        // Coleta dependências em ordem
        std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>> depsForThis;
        for (const auto& req : comp->requiredAttributes()) {
            depsForThis.push_back(available.at(req));
        }

        // Computa os atributos no buffer final
        comp->compute(tree, buffer, attrNames, depsForThis);

        // Marca atributos como disponíveis
        for (const auto& attr : comp->attributes()) {
            available[attr] = {attrNames, buffer};
        }
    }

    return {attrNames, buffer};
}