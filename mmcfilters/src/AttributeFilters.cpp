
#include "../include/AttributeFilters.hpp"


    AttributeFilters::AttributeFilters(MorphologicalTreePtr tree){
        this->tree = tree;
    }

    AttributeFilters::~AttributeFilters(){
        
    }
                           
    ImageUInt8Ptr AttributeFilters::filteringByPruningMin(std::shared_ptr<float[]> attribute, float threshold){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMin(this->tree, attribute, threshold, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFilters::filteringByExtinctionValue(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, int k){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByExtinctionValue(this->tree, attribute, k, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFilters::filteringByPruningMax(std::shared_ptr<float[]> attribute, float threshold){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMax(this->tree, attribute, threshold, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFilters::filteringByPruningMin(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMin(this->tree, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFilters::filteringByDirectRule(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByDirectRule(this->tree, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr AttributeFilters::filteringByPruningMax(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());

        AttributeFilters::filteringByPruningMax(this->tree, criterion, imgOutput);

        return imgOutput;

    }

    ImageUInt8Ptr AttributeFilters::filteringBySubtractiveRule(std::vector<bool>& criterion){

        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringBySubtractiveRule(this->tree, criterion, imgOutput);

        return imgOutput;

    }

    ImageFloatPtr AttributeFilters::filteringBySubtractiveScoreRule(std::vector<float>& prob){
        ImageFloatPtr imgOutput = ImageFloat::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringBySubtractiveScoreRule(this->tree, prob, imgOutput);
        return imgOutput;

    }

    std::vector<bool> AttributeFilters::getAdaptativeCriterion(std::vector<bool>& criterion, int delta){
        return AttributeFilters::getAdaptativeCriterion(this->tree, criterion, delta);
    }
