#pragma once

#include "../utils/Common.hpp"
#include "AttributeNames.hpp"


namespace mmcfilters {

class MorphologicalTree;
class AttributeNames;

class AttributeComputer {
public:
    virtual ~AttributeComputer() = default;

    virtual void compute(
        MorphologicalTree* tree,
        std::shared_ptr<float[]> buffer,
        std::shared_ptr<AttributeNames> attrNames,
        const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const
    {
        compute(tree, buffer, attrNames, this->attributes(), dependencySources);
    }

    virtual void compute(
        MorphologicalTree* tree,
        std::shared_ptr<float[]> buffer,
        std::shared_ptr<AttributeNames> attrNames,
        const std::vector<Attribute>& requestedAttributes,
        const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const = 0;

    virtual std::vector<Attribute> attributes() const = 0;

    virtual std::vector<AttributeOrGroup> requiredAttributes() const { return {}; }
};

} // namespace mmcfilters

