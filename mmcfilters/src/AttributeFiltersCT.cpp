
#include "../include/AttributeFiltersCT.hpp"


    AttributeFiltersCT::AttributeFiltersCT(ComponentTree* tree){
        this->tree = tree;
    }

    AttributeFiltersCT::AttributeFiltersCT(ComponentTreePtr tree){
        this->tree = tree.get();
    }

    AttributeFiltersCT::~AttributeFiltersCT(){
        
    }
                           
    ImageUInt8Ptr AttributeFiltersCT::filteringByPruningMin(std::shared_ptr<float[]> attribute, float threshold){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringByPruningMin(this->tree, attribute, threshold, imgOutput);
        return imgOutput;
    }

    /*
    ImageUInt8Ptr AttributeFiltersCT::filteringByExtinctionValue(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, int k){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringByExtinctionValue(this->tree, attribute, k, imgOutput);
        return imgOutput;
    }*/

    ImageUInt8Ptr AttributeFiltersCT::filteringByPruningMax(std::shared_ptr<float[]> attribute, float threshold){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringByPruningMax(this->tree, attribute, threshold, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFiltersCT::filteringByPruningMin(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringByPruningMin(this->tree, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFiltersCT::filteringByDirectRule(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringByDirectRule(this->tree, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFiltersCT::filteringByPruningMax(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());

        AttributeFiltersCT::filteringByPruningMax(this->tree, criterion, imgOutput);

        return imgOutput;

    }

    ImageUInt8Ptr AttributeFiltersCT::filteringBySubtractiveRule(std::vector<bool>& criterion){

        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringBySubtractiveRule(this->tree, criterion, imgOutput);

        return imgOutput;

    }

    ImageFloatPtr AttributeFiltersCT::filteringBySubtractiveScoreRule(std::vector<float>& prob){
        ImageFloatPtr imgOutput = ImageFloat::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFiltersCT::filteringBySubtractiveScoreRule(this->tree, prob, imgOutput);
        return imgOutput;

    }

    std::vector<bool> AttributeFiltersCT::getAdaptativeCriterion(std::vector<bool>& criterion, int delta){
        return AttributeFiltersCT::getAdaptativeCriterion(this->tree, criterion, delta);
    }
