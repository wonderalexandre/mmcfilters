
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
