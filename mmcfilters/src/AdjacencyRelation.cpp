#include "../include/AdjacencyRelation.hpp"
#include "../include/ImageUtils.hpp"

#include <math.h>
#define PI 3.14159265358979323846

AdjacencyRelation::~AdjacencyRelation(){

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
    this->offsetCol = std::unique_ptr<int[]>(new int[n]);
    this->offsetRow = std::unique_ptr<int[]>(new int[n]);

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
	std::unique_ptr<double[]> da(new double[this->n]);
	std::unique_ptr<double[]> dr(new double[this->n]);

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
		if (0 <= this->row + this->offsetRow[this->id] && this->row + this->offsetRow[this->id] < this->numRows && 
			0 <= this->col + this->offsetCol[this->id] && this->col + this->offsetCol[this->id] < this->numCols){
			break;
		}
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

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int index){
	auto [row, col] = ImageUtils::to2D(index, this->numCols);
    return getAdjPixels(row, col);
}

AdjacencyRelation& AdjacencyRelation::getNeighboringPixels(int row, int col){
    this->row = row;
    this->col = col;
    this->id = 0;
    return *this;
}

AdjacencyRelation& AdjacencyRelation::getNeighboringPixels(int index){
	auto [row, col] = ImageUtils::to2D(index, this->numCols);
    return getNeighboringPixels(row, col);
}

bool AdjacencyRelation::isBorderDomainImage(int index){
	auto[row, col] = ImageUtils::to2D(index, this->numCols);
	return isBorderDomainImage(row, col);
}
bool AdjacencyRelation::isBorderDomainImage(int row, int col){
	return row == 0 || col == 0 || row == this->numRows - 1 || col == this->numCols - 1;
}

bool AdjacencyRelation::isValid(int index){
	auto [row, col] = ImageUtils::to2D(index, this->numCols);
	return isValid(row, col);
}

bool AdjacencyRelation::isValid(int row, int col){
	return row >= 0 && col >= 0 && row < this->numRows && col < this->numCols;
		
}
