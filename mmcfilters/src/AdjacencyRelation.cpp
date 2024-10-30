#include "../include/AdjacencyRelation.hpp"
#include <math.h>
#define PI 3.14159265358979323846

AdjacencyRelation::~AdjacencyRelation(){
    delete[] this->offsetCol;
    delete[] this->offsetRow;
}

AdjacencyRelation::AdjacencyRelation(int numRows, int numCols, double radius){
    this->numRows = numRows;
    this->numCols = numCols;
 
    int i, j, k, dx, dy, r0, r2, i0 = 0;
    this->n = 0;
    r0 = (int) radius;
    r2 = (int) (radius * radius);
	for (dy = -r0; dy <= r0; dy++)
		for (dx = -r0; dx <= r0; dx++)
			if (((dx * dx) + (dy * dy)) <= r2)
				this->n++;
	
	i = 0;
    this->offsetCol = new int[this->n];
    this->offsetRow = new int[this->n];

	for (dy = -r0; dy <= r0; dy++) {
		for (dx = -r0; dx <= r0; dx++) {
			if (((dx * dx) + (dy * dy)) <= r2) {
				this->offsetCol[i] =dx;
                this->offsetRow[i] =dy;
				if ((dx == 0) && (dy == 0))
					i0 = i;
				i++;
			}
		}
	}
		
	double aux;
	double da[this->n];
	double dr[this->n];

	/* Set clockwise */
	for (i = 0; i < n; i++) {
		dx = this->offsetCol[i];
		dy = this->offsetRow[i];
		dr[i] = sqrt((dx * dx) + (dy * dy));
		if (i != i0) {
			da[i] = (atan2(-dy, -dx) * 180.0 / PI);
			if (da[i] < 0.0)
				da[i] += 360.0;
		}
	}
	da[i0] = 0.0;
	dr[i0] = 0.0;

	/* place central pixel at first */
	aux = da[i0];
	da[i0] = da[0];
	da[0] = aux;

	aux = dr[i0];
	dr[i0] = dr[0];
	dr[0] = aux;

	int auxX, auxY;
	auxX = this->offsetCol[i0];
	auxY = this->offsetRow[i0];
	this->offsetCol[i0] = this->offsetCol[0];
	this->offsetRow[i0] = this->offsetRow[0];
		
	this->offsetCol[0] = auxX;
	this->offsetRow[0] = auxY;
		

	/* sort by angle */
	for (i = 1; i < n - 1; i++) {
		k = i;
		for (j = i + 1; j < n; j++)
			if (da[j] < da[k]) {
				k = j;
			}
		aux = da[i];
		da[i] = da[k];
		da[k] = aux;
		aux = dr[i];
		dr[i] = dr[k];
		dr[k] = aux;

		auxX = this->offsetCol[i];
		auxY = this->offsetRow[i];
		this->offsetCol[i] = this->offsetCol[k];
		this->offsetRow[i] = this->offsetRow[k];
			
		this->offsetCol[k] = auxX;
		this->offsetRow[k] = auxY;
	}

	/* sort by radius for each angle */
	for (i = 1; i < n - 1; i++) {
		k = i;
		for (j = i + 1; j < n; j++)
			if ((dr[j] < dr[k]) && (da[j] == da[k])) {
				k = j;
			}
		aux = dr[i];
		dr[i] = dr[k];
		dr[k] = aux;

		auxX = this->offsetCol[i];
		auxY = this->offsetRow[i];
		this->offsetCol[i] = this->offsetCol[k];
		this->offsetRow[i] = this->offsetRow[k];
			
		this->offsetCol[k] = auxX;
		this->offsetRow[k] = auxY;
			
	}

    
}

int AdjacencyRelation::getSize(){
	return this->n;
}

int AdjacencyRelation::nextValid(){
    this->id += 1;
    while (this->id < this->n){
        if (0 <= this->row + this->offsetRow[this->id] && this->row + this->offsetRow[this->id] < this->numRows && 0 <= this->col + this->offsetCol[this->id] && this->col + this->offsetCol[this->id] < this->numCols)
            break;
        this->id += 1;
    }
    return this->id;
} 

AdjacencyRelation::IteratorAdjacency AdjacencyRelation::begin() { 
    return IteratorAdjacency(*this, nextValid()); 
}

AdjacencyRelation::IteratorAdjacency AdjacencyRelation::end() { 
    return IteratorAdjacency(*this, this->n); 
}

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int row, int col){
    this->row = row;
    this->col = col;
    this->id = -1;
    return *this;
}

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int indexVector){
    return getAdjPixels(indexVector / this->numCols, indexVector % this->numCols);
}
