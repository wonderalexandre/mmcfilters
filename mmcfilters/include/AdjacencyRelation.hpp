
#include "../include/Common.hpp"

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
 
    std::unique_ptr<int[]> offsetRow;
    std::unique_ptr<int[]> offsetCol;
    
      

public:

    AdjacencyRelation(int numRows, int numCols, double radius);
    ~AdjacencyRelation();
    int nextValid();
    int getSize();
    AdjacencyRelation& getAdjPixels(int row, int col);
    AdjacencyRelation& getAdjPixels(int index);

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
                return (instance.row + instance.offsetRow[index]) * instance.numCols + (instance.col + instance.offsetCol[index]); 
            }    
    };
    IteratorAdjacency begin();
    IteratorAdjacency end();	 
};

#endif
