
#include "../include/AttributeFilters.hpp"


    AttributeFilters::AttributeFilters(MorphologicalTreePtr tree){
        this->tree = tree;
    }

    AttributeFilters::~AttributeFilters(){
        
    }
                           
    int* AttributeFilters::filteringByPruningMin(float* attribute, float threshold){
        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        int* imgOutput = new int[n];
        AttributeFilters::filteringByPruningMin(this->tree, attribute, threshold, imgOutput);
        return imgOutput;
    }

    int* AttributeFilters::filteringByPruningMax(float* attribute, float threshold){
        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        int* imgOutput = new int[n];
        AttributeFilters::filteringByPruningMax(this->tree, attribute, threshold, imgOutput);
        return imgOutput;
    }

    int* AttributeFilters::filteringByPruningMin(std::vector<bool>& criterion){
        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        int* imgOutput = new int[n];
        AttributeFilters::filteringByPruningMin(this->tree, criterion, imgOutput);
        return imgOutput;
    }

    int* AttributeFilters::filteringByDirectRule(std::vector<bool>& criterion){
        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        int* imgOutput = new int[n];
        AttributeFilters::filteringByDirectRule(this->tree, criterion, imgOutput);
        return imgOutput;
    }

    int* AttributeFilters::filteringByPruningMax(std::vector<bool>& criterion){
        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        int* imgOutput = new int[n];

        AttributeFilters::filteringByPruningMax(this->tree, criterion, imgOutput);

        return imgOutput;

    }

    int* AttributeFilters::filteringBySubtractiveRule(std::vector<bool>& criterion){

        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        int* imgOutput = new int[n];

        AttributeFilters::filteringBySubtractiveRule(this->tree, criterion, imgOutput);

        return imgOutput;

    }

    float* AttributeFilters::filteringBySubtractiveScoreRule(std::vector<float>& prob){
        int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
        float* imgOutput = new float[n];

        AttributeFilters::filteringBySubtractiveScoreRule(this->tree, prob, imgOutput);

        return imgOutput;

    }

    std::vector<bool> AttributeFilters::getAdaptativeCriterion(std::vector<bool>& criterion, int delta){
        return AttributeFilters::getAdaptativeCriterion(this->tree, criterion, delta);
    }
