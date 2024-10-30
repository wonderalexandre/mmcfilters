
#include "../include/AttributeComputedIncrementally.hpp"


void AttributeComputedIncrementally::preProcessing(NodeCT *v){}

void AttributeComputedIncrementally::mergeChildren(NodeCT *parent, NodeCT *child){}

void AttributeComputedIncrementally::postProcessing(NodeCT *parent){}

void AttributeComputedIncrementally::computerAttribute(NodeCT *root) {
        preProcessing(root);
        for (NodeCT *child : root->getChildren())
        {
            computerAttribute(child);
            mergeChildren(root, child);
        }
        postProcessing(root);
}
