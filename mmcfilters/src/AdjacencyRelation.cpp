#include "../include/AdjacencyRelation.hpp"
#include <math.h>
#include <cmath> 
#include <stdexcept>
#define PI 3.14159265358979323846


AdjacencyRelation::AdjacencyRelation(int numRows, int numCols, double radius){
    this->numRows = numRows;
    this->numCols = numCols;
	this->radius = radius;

    int i, j, k, dx, dy, r0, r2, i0 = 0;
    this->n = 0;
    r0 = (int) radius;
    r2 = (int) (radius * radius);
	for (dy = -r0; dy <= r0; dy++)
		for (dx = -r0; dx <= r0; dx++)
			if (((dx * dx) + (dy * dy)) <= r2)
				this->n++;
	
	i = 0;
    this->offsetCol.resize(this->n);
    this->offsetRow.resize(this->n);
	
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

double AdjacencyRelation::getRadius(){
	return this->radius;
}

bool AdjacencyRelation::isAdjacent(int px, int py, int qx, int qy) {
	double distance = std::sqrt(std::pow(px - qx, 2) + std::pow(py - qy, 2));
    return (distance <= radius);
}

bool AdjacencyRelation::isAdjacent(int p, int q) {
    int py = p / numCols, px = p % numCols;
    int qy = q / numCols, qx = q % numCols;

    return isAdjacent(px, py, qx, qy);
}

int AdjacencyRelation::nextValid() {
    this->id += 1;
    while (this->id < this->n) {
        int newRow = this->row + this->offsetRow[this->id];
        int newCol = this->col + this->offsetCol[this->id];

        if (newRow >= 0 && newRow < this->numRows && newCol >= 0 && newCol < this->numCols) {
            return this->id;
        }
        this->id += 1;
    }
    return this->n;
}

AdjacencyRelation::IteratorAdjacency AdjacencyRelation::begin() { 
    return IteratorAdjacency(this, nextValid()); 
}

AdjacencyRelation::IteratorAdjacency AdjacencyRelation::end() { 
    return IteratorAdjacency(this, this->n); 
}

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int row, int col){
	if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Índice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = -1;
    return *this;
}

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int indexVector){
    return getAdjPixels(indexVector / this->numCols, indexVector % this->numCols);
}
