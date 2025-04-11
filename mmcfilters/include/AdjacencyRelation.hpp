
#include "../include/Common.hpp"
#include "../include/ImageUtils.hpp"

#ifndef ADJACENCY_H
#define ADJACENCY_H

class AdjacencyRelation;  // forward declaration
using AdjacencyRelationPtr = std::shared_ptr<AdjacencyRelation>;

class AdjacencyRelation {
private:
    int id;
    
    int row;
    int col;    
    int numCols;
    int numRows;
    int n;
    bool validatePixels;
 
    std::unique_ptr<int[]> offsetRow;
    std::unique_ptr<int[]> offsetCol;
    
      

public:

    AdjacencyRelation(int numRows, int numCols, double radius);
    ~AdjacencyRelation();
    int nextValid();
    int getSize();
    AdjacencyRelation& getAdjPixels(int row, int col);
    AdjacencyRelation& getAdjPixels(int index);
    AdjacencyRelation& getNeighboringPixels(int row, int col);
    AdjacencyRelation& getNeighboringPixels(int index);
    bool isValid(int index);
    bool isValid(int row, int col);
    bool isBorderDomainImage(int index);
    bool isBorderDomainImage(int row, int col);
    class IteratorAdjacency{ 
        private:
    	    int index;
            AdjacencyRelation&  instance;
        public:
        	using iterator_category = std::input_iterator_tag;
            using value_type = int; 
            
            IteratorAdjacency(AdjacencyRelation& obj, int id): instance(obj), index(id)  { }

            IteratorAdjacency& operator++() { 
                this->index = instance.nextValid(); return *this; 
            }
            bool operator==(const IteratorAdjacency& other) const {
                return index == other.index;
            }
            bool operator!=(const IteratorAdjacency& other) const {
                return index != other.index;
            }

            int operator*() const { 
                return ImageUtils::to1D(instance.row + instance.offsetRow[index], instance.col + instance.offsetCol[index], instance.numCols);
                //return (instance.row + instance.offsetRow[index]) * instance.numCols + (instance.col + instance.offsetCol[index]); 
            }    
    };
    IteratorAdjacency begin();
    IteratorAdjacency end();	 
};

#endif
